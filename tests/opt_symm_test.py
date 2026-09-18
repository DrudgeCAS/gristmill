"""Tests for optimizations using the symmetry of the left-hand sides."""

import pytest
from drudge import PartHoleDrudge, Perm, NEG
from sympy import IndexedBase, Rational

from gristmill import (
    optimize,
    verify_eval_seq,
    get_flop_cost,
    EinsumPrinter,
)


@pytest.fixture(scope="module")
def parthole_drudge(spark_ctx):
    """The particle-hole drudge."""
    dr = PartHoleDrudge(spark_ctx)
    return dr


@pytest.fixture(scope="module")
def ccd_doubles(parthole_drudge):
    """The CCD doubles amplitude equation, derived from scratch.

    The residual is declared to have the double-bar symmetry of the amplitudes.
    """

    dr = parthole_drudge
    p = dr.names
    a, b = p.V_dumms[:2]
    i, j = p.O_dumms[:2]
    c_ = p.c_
    c_dag = p.c_dag

    t = IndexedBase("t")
    dr.set_dbbar_base(t, 2)
    t2 = dr.einst(
        Rational(1, 4) * t[a, b, i, j] * c_dag[a] * c_dag[b] * c_[j] * c_[i]
    )

    curr = dr.ham
    h_bar = dr.ham
    for order in range(2):
        curr = (curr | t2).simplify() / (order + 1)
        h_bar += curr
        continue
    h_bar = h_bar.simplify()

    proj = c_dag[i] * c_dag[j] * c_[b] * c_[a]
    r = IndexedBase("r")
    dr.set_symm(r, Perm([1, 0, 2, 3], NEG), Perm([0, 1, 3, 2], NEG), valence=4)
    return dr.define(r[a, b, i, j], (proj * h_bar).eval_fermi_vev().simplify())


def _leading_coeff(eval_seq, no, nv, o_power, v_power):
    """Get the coefficient of a scaling term in the leading FLOP cost."""
    cost = get_flop_cost(eval_seq, leading=True).expand()
    return cost.coeff(no**o_power * nv**v_power)


def test_ccd_doubles_lhs_symm(ccd_doubles):
    """Test the CCD doubles equation with the symmetry of its left-hand side.

    Without using the symmetry, three separate contractions scaling as
    :math:`o^3 v^3` are needed for the ring terms, while with the symmetry the
    linear and quadratic ring terms can share a common intermediate, which is
    the scheme in the literature.
    """

    tdef = ccd_doubles
    dr = tdef.drudge
    p = dr.names
    no, nv = p.no, p.nv
    substs = {nv: 10 * no}

    plain_seq = optimize([tdef], substs=substs)
    assert verify_eval_seq(plain_seq, [tdef], simplify=True)
    assert _leading_coeff(plain_seq, no, nv, 3, 3) == 6

    symm_seq = optimize([tdef], substs=substs, lhs_symm=True)
    assert verify_eval_seq(symm_seq, [tdef], simplify=True)
    assert _leading_coeff(symm_seq, no, nv, 3, 3) == 4
    assert _leading_coeff(symm_seq, no, nv, 4, 2) == 4
    assert _leading_coeff(symm_seq, no, nv, 2, 4) == 2

    # The seed is an intermediate, and the result is assembled at the end by
    # four permuted copies of the seed.
    final = symm_seq[-1]
    assert final.base == tdef.base
    assert not final.if_interm
    assert final.n_terms == 4
    assert all(len(i.sums) == 0 for i in final.local_terms)

    seed = symm_seq[-2]
    assert str(seed.base) == "r_s"
    assert seed.if_interm

    # The printers should be able to handle the sequence.
    assert "r_s" in EinsumPrinter().doprint(symm_seq)


def test_independent_outputs_share_intermediates(parthole_drudge):
    """Test independent symmetric outputs optimized together.

    Two residual-like outputs with the same doubles ladder contraction should
    share its intermediate when they are optimized together.
    """

    dr = parthole_drudge
    p = dr.names
    a, b, c, d = p.V_dumms[:4]
    i, j, k, l = p.O_dumms[:4]
    no, nv = p.no, p.nv
    u = dr.two_body
    t = IndexedBase("t")
    dr.set_dbbar_base(t, 2)

    ladder = t[a, b, k, l] * t[c, d, i, j] * u[k, l, c, d]
    ring = t[a, c, i, k] * u[k, b, c, j]

    x1 = IndexedBase("x1")
    x2 = IndexedBase("x2")
    for base in [x1, x2]:
        dr.set_symm(
            base, Perm([1, 0, 2, 3], NEG), Perm([0, 1, 3, 2], NEG), valence=4
        )
        continue

    def antisymm(expr):
        swap_ab = expr.xreplace({a: b, b: a})
        return (
            expr
            - swap_ab
            - expr.xreplace({i: j, j: i})
            + swap_ab.xreplace({i: j, j: i})
        )

    def1 = dr.define(x1[a, b, i, j], dr.einst(ladder / 4 + antisymm(ring)))
    def2 = dr.define(x2[a, b, i, j], dr.einst(ladder / 2 + u[a, b, i, j]))
    targets = [def1, def2]

    seq = optimize(targets, substs={nv: 10 * no}, lhs_symm=True)
    assert verify_eval_seq(seq, targets, simplify=True)

    # Both results are assembled, both seeds are intermediates.
    finals = [i for i in seq if not i.if_interm]
    assert [i.base for i in finals] == [x1, x2]
    interms = [i for i in seq if i.if_interm]
    assert {str(i.base) for i in interms} >= {"x1_s", "x2_s"}

    # The ladder contraction is computed only once when the two outputs are
    # optimized together, while separate optimizations need it twice.
    assert _leading_coeff(seq, no, nv, 4, 2) == 4
    sep_seqs = [
        optimize([i], substs={nv: 10 * no}, lhs_symm=True) for i in targets
    ]
    assert sum(_leading_coeff(i, no, nv, 4, 2) for i in sep_seqs) == 8


@pytest.mark.parametrize("res_at_end", [True, False])
def test_dependent_outputs_are_assembled_before_use(
    parthole_drudge, res_at_end
):
    """Test symmetric outputs consumed by other outputs.

    The assembly of the symmetric result must come before the computation of
    a result using it, and the seed must not be freed before its assembly.
    """

    dr = parthole_drudge
    p = dr.names
    a, b, c = p.V_dumms[:3]
    i, j, k = p.O_dumms[:3]
    no, nv = p.no, p.nv
    u = dr.two_body
    f = dr.fock
    t = IndexedBase("t")
    dr.set_dbbar_base(t, 2)

    w = IndexedBase("w")
    dr.set_symm(w, Perm([1, 0, 2, 3], NEG), Perm([0, 1, 3, 2], NEG), valence=4)
    ring = t[a, c, i, k] * u[k, b, c, j]
    swap_ab = ring.xreplace({a: b, b: a})
    w_def = dr.define(
        w[a, b, i, j],
        dr.einst(
            ring
            - swap_ab
            - ring.xreplace({i: j, j: i})
            + swap_ab.xreplace({i: j, j: i})
        ),
    )

    # A plain output consuming the symmetric one, together with a plain output
    # independent of it, to test the pass-through.
    y = IndexedBase("y")
    y_def = dr.define(y[a, b, i, j], dr.einst(w[a, b, i, k] * f[k, j]))
    z = IndexedBase("z")
    z_def = dr.define(z[a, i], dr.einst(f[a, k] * t[k, i]))

    targets = [y_def, w_def, z_def]
    seq = optimize(
        targets, substs={nv: 10 * no}, res_at_end=res_at_end, lhs_symm=True
    )
    assert verify_eval_seq(seq, targets, simplify=True)

    bases = [i.base for i in seq]
    assert bases.index(w) < bases.index(y)
    assert bases.index(IndexedBase("w_s")) < bases.index(w)
    if res_at_end:
        assert bases[-3:] == [w, y, z]

    # Pass-through of the plain outputs.
    for def_ in seq:
        if def_.base in (y, z):
            assert not def_.if_interm
        continue
    z_final = next(i for i in seq if i.base == z)
    assert z_final == z_def.simplify()

    # The printer verifies the dependency order and the lifetime of the seed.
    code = EinsumPrinter().doprint(seq)
    assert code.index("del w_s") > code.index("w_s", code.index("w = zeros"))
