"""Tests for the block storage of input tensors in generated code."""

import numpy as np
import pytest
import os
import pickle
import subprocess
import tempfile
from drudge import PartHoleDrudge, Range, Perm, NEG
from sympy import IndexedBase, Rational, Symbol

from gristmill import (
    optimize,
    verify_eval_seq,
    EinsumPrinter,
    BlockSpec,
    split_input_blocks,
)


@pytest.fixture(scope="module")
def offset_parthole(spark_ctx):
    """The particle-hole drudge with the virtual range after the occupied.

    In this way, the indices in both ranges can be used directly to slice full
    arrays of the whole orbital space.
    """

    no = Symbol("no")
    nv = Symbol("nv")
    dr = PartHoleDrudge(
        spark_ctx,
        part_orb=(Range("V", no, no + nv), PartHoleDrudge.DEFAULT_PART_DUMMS),
    )
    return dr


@pytest.fixture(scope="module")
def ccd_doubles(offset_parthole):
    """The CCD doubles amplitude equation with symmetric left-hand side."""

    dr = offset_parthole
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


def _dbbar_random(rng, shape):
    """Form a random array with the double-bar symmetry.

    The symbolic manipulations rely on the declared antisymmetry of the
    tensors, so the numerical arrays need to carry it as well.
    """

    arr = rng.random(shape)
    arr = arr - arr.transpose(1, 0, 2, 3)
    arr = arr - arr.transpose(0, 1, 3, 2)
    if shape[0] == shape[2]:
        arr = arr + arr.transpose(2, 3, 0, 1)
    return arr


def _run_einsum(code, **arrays):
    """Run code, optionally checking PySCF in a separate Python environment."""
    env = {"einsum": np.einsum}
    env.update(arrays)
    exec(code, env)
    python = os.environ.get("PYSCF_PYTHON")
    if python:
        with tempfile.TemporaryDirectory() as directory:
            inp = os.path.join(directory, "input.pickle")
            out = os.path.join(directory, "output.pickle")
            with open(inp, "wb") as stream:
                pickle.dump((code, arrays), stream)
            subprocess.run(
                [
                    python,
                    "-c",
                    """
import pickle, sys
from pyscf.lib import einsum
with open(sys.argv[1], 'rb') as stream:
    code, arrays = pickle.load(stream)
env = dict(arrays, einsum=einsum)
exec(code, env)
result = {k: v for k, v in env.items()
          if not k.startswith('__') and not callable(v)}
with open(sys.argv[2], 'wb') as stream:
    pickle.dump(result, stream)
""",
                    inp,
                    out,
                ],
                check=True,
            )
            with open(out, "rb") as stream:
                actual = pickle.load(stream)
        for key, value in actual.items():
            assert np.allclose(value, env[key]), key
        return actual
    return env


def test_input_blocks_are_formed_and_released(offset_parthole):
    """Test the formation and release of blocks in generated code.

    Blocks of the two-body interaction should be formed right before their
    first use, named by the ranges of the indices, and released after their
    last use.
    """

    dr = offset_parthole
    p = dr.names
    a, b, c = p.V_dumms[:3]
    i, j = p.O_dumms[:2]
    u = dr.two_body
    t1 = IndexedBase("t1")

    tdef = dr.define_einst(
        IndexedBase("w")[a, b, i, j],
        u[a, b, i, c] * t1[c, j] - u[a, b, j, c] * t1[c, i],
    )
    eval_seq = optimize([tdef], substs={p.nv: 10 * p.no}, interm_fmt="tau{}")
    assert verify_eval_seq(eval_seq, [tdef])

    # The sequence transform alone.
    split_seq = split_input_blocks(eval_seq, BlockSpec([u]))
    blocks = [i for i in split_seq if hasattr(i, "block")]
    assert len(blocks) == 1
    block = blocks[0]
    assert str(block.base) == "u_vvov"
    assert block.if_interm
    assert block.block.base == "u"
    assert [i.label for i in block.block.ranges] == ["V", "V", "O", "V"]
    assert [str(i[0]) for i in block.exts] == ["a", "b", "i", "c"]
    # The block comes right before its first use, which no longer refers to u.
    users = [i for i in split_seq if i.free_vars & {Symbol("u_vvov")}]
    assert len(users) == 1
    assert split_seq.index(block) == split_seq.index(users[0]) - 1
    assert not any(
        i.free_vars & {Symbol("u")}
        for i in split_seq
        if not hasattr(i, "block")
    )

    # The generated code.
    printer = EinsumPrinter(blocks=[u], base_indent=0)
    code = printer.doprint(eval_seq)
    lines = code.splitlines()
    form = "u_vvov = u[no:no + nv, no:no + nv, 0:no, no:no + nv]"
    assert form in lines
    first_use = min(
        idx
        for idx, line in enumerate(lines)
        if "u_vvov" in line and line.strip() != form
    )
    last_use = max(
        idx
        for idx, line in enumerate(lines)
        if "u_vvov" in line and line.strip() != "del u_vvov"
    )
    assert lines.index(form) < first_use
    assert lines.index("del u_vvov") > last_use
    assert "einsum(" not in "".join(
        line for line in lines if "u_vvov = " in line
    )

    # Numerical check against direct slicing.
    no, nv = 3, 4
    n = no + nv
    rng = np.random.default_rng(1)
    u_arr = _dbbar_random(rng, (n, n, n, n))
    t1_arr = rng.random((nv, no))
    env = _run_einsum(code, no=no, nv=nv, u=u_arr, t1=t1_arr)
    u_vvov = u_arr[no:, no:, :no, no:]
    w_ref = np.einsum("abic,cj->abij", u_vvov, t1_arr)
    w_ref -= np.einsum("abjc,ci->abij", u_vvov, t1_arr)
    assert np.allclose(env["w"], w_ref)
    assert "u_vvov" not in env  # Released.

    # Blocks given as inputs: no formation, no release.
    printer = EinsumPrinter(blocks=BlockSpec([u], src_fmt=None), base_indent=0)
    code = printer.doprint(eval_seq)
    assert "u_vvov" in code
    assert "u_vvov = " not in code
    assert "del u_vvov" not in code
    env = _run_einsum(code, no=no, nv=nv, u_vvov=u_vvov, t1=t1_arr)
    assert np.allclose(env["w"], w_ref)

    # Customized naming and formation.
    spec = BlockSpec(
        [u],
        name_fmt="{base}{blocks}",
        range_labels={"O": "h", "V": "p"},
        src_fmt="get_block('{base}', '{blocks}')",
    )
    code = EinsumPrinter(blocks=spec, base_indent=0).doprint(eval_seq)
    assert "uppHp = get_block('u', 'ppHp')".replace("H", "h") in code


def test_unchanged_definition_and_permutation(offset_parthole):
    """Untouched definitions retain metadata; permutations use safe views."""
    dr = offset_parthole
    a, b = dr.names.V_dumms[:2]
    t, r = IndexedBase("t"), IndexedBase("r")
    def_ = dr.define_einst(r[a, b], t[a, b] - t[b, a])
    def_.if_interm = False
    assert split_input_blocks([def_], BlockSpec([dr.two_body]))[0] is def_
    code = EinsumPrinter(base_indent=0).doprint([def_])
    assert "einsum(" not in code
    assert ".transpose((1, 0))" in code
    # Exercise a noncontiguous input view and verify it is not modified.
    arr = np.arange(32.0).reshape(4, 8)[:, ::2]
    original = arr.copy()
    env = _run_einsum(code, t=arr, no=2, nv=4)
    assert np.array_equal(env["r"], original - original.T)
    assert np.array_equal(arr, original)


def test_slicing_needs_distinct_lower_bounds(spark_ctx):
    """Test the rejection of ambiguous slicing with zero-based ranges."""

    dr = PartHoleDrudge(
        spark_ctx,
        part_orb=(
            Range("V", 0, Symbol("nv")),
            PartHoleDrudge.DEFAULT_PART_DUMMS,
        ),
    )
    p = dr.names
    a, b = p.V_dumms[:2]
    i, j = p.O_dumms[:2]
    u = dr.two_body
    # Occupied and virtual indices share the axes of u across the terms.
    tdef = dr.define_einst(
        IndexedBase("w")[a, b, i, j], u[a, b, i, j] + u[i, j, a, b]
    )
    with pytest.raises(ValueError):
        EinsumPrinter(blocks=[u]).doprint([tdef])
    # Blocks given as inputs do not need slicing.
    code = EinsumPrinter(blocks=BlockSpec([u], src_fmt=None)).doprint([tdef])
    assert "u_vvoo" in code


def test_ccd_doubles_with_blocks(ccd_doubles):
    """Test the block storage on the optimized CCD doubles equation.

    The generated code with blocks is compared numerically with the code for
    the unoptimized equation, and each block should be formed once and
    released after use.
    """

    tdef = ccd_doubles
    dr = tdef.drudge
    p = dr.names
    u = dr.two_body
    f = dr.fock

    eval_seq = optimize(
        [tdef], substs={p.nv: 10 * p.no}, interm_fmt="tau{}", lhs_symm=True
    )
    assert verify_eval_seq(eval_seq, [tdef], simplify=True)

    printer = EinsumPrinter(blocks=[u, f], base_indent=0)
    code = printer.doprint(eval_seq)
    ref_code = printer.doprint([tdef])

    # Every block is formed exactly once, released exactly once, and u itself
    # is never contracted directly.
    lines = [i.strip() for i in code.splitlines()]
    formed = [i.split(" = ")[0] for i in lines if " = u[" in i or " = f[" in i]
    assert len(formed) == len(set(formed))
    assert {"u_oovv", "u_oooo", "u_vvvv", "u_ovov"} <= set(formed)
    for name in formed:
        assert lines.count("del {}".format(name)) == 1
        assert lines.index("del {}".format(name)) > max(
            idx
            for idx, line in enumerate(lines)
            if name in line and not line.startswith("del ")
        )
        continue
    assert not any(", u)" in i or ", u," in i for i in lines)

    no, nv = 3, 4
    n = no + nv
    rng = np.random.default_rng(7)
    arrays = dict(
        no=no,
        nv=nv,
        u=_dbbar_random(rng, (n, n, n, n)),
        f=rng.random((n, n)),
        t=_dbbar_random(rng, (nv, nv, no, no)),
    )
    env = _run_einsum(code, **arrays)
    ref_env = _run_einsum(ref_code, **arrays)
    assert np.allclose(env["r"], ref_env["r"])
    assert not any(i in env for i in formed)


def test_blocks_allow_shared_bounds_in_different_axes(spark_ctx):
    """Test block slicing for tensors mixing a boson and orbital ranges.

    The boson range and the occupied range both start from zero, which is
    harmless since they never share an axis of the same tensor.
    """

    from drudge import Drudge

    dr = Drudge(spark_ctx)
    no, nv, nb = Symbol("no"), Symbol("nv"), Symbol("nb")
    o = Range("O", 0, no)
    v = Range("V", no, no + nv)
    bo = Range("B", 0, nb)
    dr.set_dumms(o, [Symbol(n) for n in "ijkl"])
    dr.set_dumms(v, [Symbol(n) for n in "abcd"])
    dr.set_dumms(bo, [Symbol(n) for n in "xy"])
    dr.add_resolver_for_dumms()
    i, j = dr.names.O_dumms[:2]
    a, b = dr.names.V_dumms[:2]
    x = dr.names.B_dumms[0]

    g = IndexedBase("G")
    t1 = IndexedBase("t1")
    tdef = dr.define_einst(IndexedBase("w")[x, a, i], g[x, a, j] * t1[j, i])
    code = EinsumPrinter(blocks=[g], base_indent=0).doprint([tdef])
    assert "G_bvo = G[0:nb, no:no + nv, 0:no]" in code

    # A genuine conflict in one axis is still rejected.
    o0 = Range("V0", 0, nv)
    dr.set_dumms(o0, [Symbol("e"), Symbol("f")])
    dr.add_resolver_for_dumms()
    e = dr.names.V0_dumms[0]
    # Zero-based V0 in the axis where O appears in the first definition.
    bad = dr.define_einst(IndexedBase("w2")[x, i, e], g[x, j, e] * t1[j, i])
    with pytest.raises(ValueError):
        EinsumPrinter(blocks=[g], base_indent=0).doprint([tdef, bad])
