"""Test optimization of different special kinds of tensor computations."""

import pytest
from drudge import Drudge, Range, TensorDef
from sympy import symbols, Symbol, IndexedBase, conjugate

from gristmill import optimize, verify_eval_seq, get_flop_cost, ContrStrat


@pytest.fixture(scope="module")
def simple_drudge(spark_ctx):
    """Make simple drudge.

    This fixture gives a simple drudge with a simple range and a few dummies.
    """

    dr = Drudge(spark_ctx)

    n = symbols("n")
    r = Range("r", 0, n)
    dumms = symbols("a b c d e f g h")
    dr.set_dumms(r, dumms)
    dr.add_default_resolver(r)

    dr.n = n
    dr.r = r
    dr.ds = dumms

    return dr


def test_summation_shared_by_four_factors(simple_drudge):
    """Test a summation involved by more than two factors.

    In a classical tensor contraction, every summation is involved by exactly
    two factors.  Here a single summation is involved by four, which the
    parenthesization search used to mishandle: it could skip every candidate
    partition and leave the subproblem with no evaluation at all, crashing the
    interpreter outright under the ``OPT`` and ``GREEDY`` strategies.

    There is only one way to evaluate this, so all strategies have to agree.
    """

    dr = simple_drudge
    a = dr.ds[0]

    x, y, z, w = (IndexedBase(i) for i in ["x", "y", "z", "w"])
    targets = [dr.define_einst(Symbol("s"), x[a] * y[a] * z[a] * w[a])]

    costs = []
    for strat in ContrStrat:
        eval_seq = optimize(targets, contr_strat=strat)
        assert verify_eval_seq(eval_seq, targets)
        costs.append(get_flop_cost(eval_seq))

    assert all(i == costs[0] for i in costs)


def test_constriction_over_an_outer_product(simple_drudge):
    """Test factoring a common tensor out over an outer product.

    Both terms share A, so the good evaluation forms X + y z once and contracts
    it with A, costing 4 n^2.  Taking the locally optimal contraction of the
    second term instead gives 4 n^2 + 2 n + 1.

    The optimizer used to return either, depending on the order the vertices of
    the constriction graph happened to be numbered in, which followed the order
    a C++ unordered map was iterated.  See issue #43.
    """

    dr = simple_drudge
    n = dr.n
    a, b = dr.ds[:2]

    A, X = IndexedBase("A"), IndexedBase("X")
    y, z = IndexedBase("y"), IndexedBase("z")

    targets = [
        dr.define_einst(Symbol("e"), A[a, b] * X[a, b] + A[a, b] * y[a] * z[b])
    ]

    eval_seq = optimize(targets)

    assert verify_eval_seq(eval_seq, targets)
    assert get_flop_cost(eval_seq) == 4 * n**2


def test_constriction_over_three_terms(simple_drudge):
    """Test a constriction that has to take all three terms at once.

    Every term has both X and y, so the good evaluation forms the vector
    z - u + 2 y, contracts X with y once, and puts the two together, costing
    2 n^2 + 4 n.  Stopping after a constriction over only two of the terms
    leaves a second matrix-vector product and costs 4 n^2 + 4 n.

    Both were reachable before, decided by vertex numbering.  See issue #43.
    """

    dr = simple_drudge
    n = dr.n
    a, b = dr.ds[:2]

    X = IndexedBase("X")
    y, z, u = (IndexedBase(i) for i in ["y", "z", "u"])

    targets = [
        dr.define_einst(
            Symbol("e"),
            X[a, b] * z[a] * y[b]
            + 2 * X[a, b] * y[a] * y[b]
            - X[a, b] * u[a] * y[b],
        )
    ]

    eval_seq = optimize(targets)

    assert verify_eval_seq(eval_seq, targets)
    assert get_flop_cost(eval_seq) == 2 * n**2 + 4 * n


def test_constriction_search_is_not_pruned_unsoundly(simple_drudge):
    """Test a constriction the pivot pruning used to throw away.

    All three terms share P and all three involve u, so the good evaluation
    contracts P once and costs 4 n^2 + 4 n.  Pruning the search by the pivot
    rule loses that and leaves 5 n^2 + 3 n, which is worse by a whole power of
    n rather than by a constant.

    The pivot rule assumes taking a vertex into a biclique cannot lower its
    saving, which is true for ordinary maximal cliques and false here.  See
    issue #43.
    """

    dr = simple_drudge
    n = dr.n
    a, b = dr.ds[:2]

    P = IndexedBase("P")
    z, u, w = (IndexedBase(i) for i in ["z", "u", "w"])

    targets = [
        dr.define_einst(
            Symbol("res"),
            -P[a, b] * z[a] * u[b]
            + P[a, b] * u[a] * w[b]
            + 3 * P[a, b] * w[a] * u[b],
        )
    ]

    eval_seq = optimize(targets)

    assert verify_eval_seq(eval_seq, targets)
    assert get_flop_cost(eval_seq) == 4 * n**2 + 4 * n


def test_sum_with_an_index_free_product(simple_drudge):
    """A sum in which one term reduces to an index-free product.

    After the shared factor is pulled out, the second term leaves a scalar
    behind, and the optimizer sends a product with no index at all to the
    parenthesization extension.  That used to segfault the interpreter
    (issue #45), through a null owning handle in cpypp.
    """

    dr = simple_drudge
    a, b, c = dr.ds[:3]

    q_mat, v = IndexedBase("Q"), IndexedBase("V")
    q, s = IndexedBase("q"), IndexedBase("s")
    targets = [
        dr.define_einst(
            Symbol("res"),
            q_mat[a, b] * q[a] * s[b]
            + 2 * q_mat[a, b] * q_mat[a, b]
            - q_mat[a, b] * v[a, c] * v[c, b],
        )
    ]

    for strat in ContrStrat:
        eval_seq = optimize(targets, contr_strat=strat)
        assert verify_eval_seq(eval_seq, targets)


def test_simple_scalar_optimization(spark_ctx):
    """Test optimization of a simple scalar.

    There is not much optimization that can be done for simple scalars.  But we
    need to ensure that we get correct result here.
    """

    dr = Drudge(spark_ctx)

    a, b, r = symbols("a b r")
    targets = [dr.define(r, a * b)]
    eval_seq = optimize(targets)
    assert verify_eval_seq(eval_seq, targets)


def test_conjugation_optimization(simple_drudge):
    """Test optimization of expressions containing complex conjugate."""

    dr = simple_drudge

    a, b, c, d = dr.ds[:4]

    p = IndexedBase("p")
    x = IndexedBase("x")
    y = IndexedBase("y")
    z = IndexedBase("z")

    targets = [
        dr.define_einst(
            p[a, b],
            x[a, c] * conjugate(y[c, b]) + x[a, c] * conjugate(z[c, b]),
        )
    ]
    eval_seq = optimize(targets)
    assert verify_eval_seq(eval_seq, targets)


def test_optimization_handles_coeffcients(simple_drudge):
    """Test optimization of scalar intermediates scaled by coefficients.

    This test comes from PoST theory.  It tests the optimization of tensor
    evaluations with scalar intermediates scaled by a factor.
    """

    dr = simple_drudge

    a, b = dr.ds[:2]

    r = IndexedBase("r")
    eps = IndexedBase("epsilon")
    t = IndexedBase("t")

    targets = [
        dr.define(r[a, b], dr.sum(2 * eps[a] * t[a, b]) - 2 * eps[b] * t[a, b])
    ]
    eval_seq = optimize(targets)
    assert verify_eval_seq(eval_seq, targets)


def test_optimization_handles_scalar_intermediates(simple_drudge):
    """Test optimization of scalar intermediates scaling other tensors.

    This is set as a special test primarily since it would entail the same
    collectible giving residues with different ranges.
    """

    dr = simple_drudge

    r = dr.r
    a, b, c = dr.ds[:3]

    u = IndexedBase("u")
    eps = IndexedBase("epsilon")
    t = IndexedBase("t")
    s = IndexedBase("s")

    targets = [
        dr.define(
            u,
            (a, r),
            (b, r),
            dr.sum((c, r), 8 * s[a, b] * eps[c] * t[a])
            - 8 * s[a, b] * eps[a] * t[a],
        )
    ]
    eval_seq = optimize(targets)
    assert verify_eval_seq(eval_seq, targets)


def test_optimization_handles_nonlinear_factors(simple_drudge):
    """Test optimization of with nonlinear factors.

    Here a factor is the square of an indexed quantity.
    """

    dr = simple_drudge

    r = dr.r
    a, b, c, d = dr.ds[:4]

    u = symbols("u")
    s = IndexedBase("s")

    targets = [
        dr.define(
            u,
            dr.sum(
                (a, r),
                (b, r),
                (c, r),
                (d, r),
                32 * s[a, c] ** 2 * s[b, d] ** 2
                + 32 * s[a, c] * s[a, d] * s[b, c] * s[b, d],
            ),
        )
    ]
    eval_seq = optimize(targets)
    assert verify_eval_seq(eval_seq, targets)


def test_common_summation_intermediate_recognition(simple_drudge):
    """Test recognition of summation intermediate differing only in a scalar."""

    dr = simple_drudge

    a, b, c = dr.ds[:3]

    x = IndexedBase("x")
    y = IndexedBase("y")
    p = IndexedBase("p")
    q = IndexedBase("q")
    r = IndexedBase("r")
    s = IndexedBase("s")

    alpha = symbols("alpha")

    for c1, c2, c3, c4 in [
        (1, 1, 1, 1),
        (1, 1, 2, 2),
        (1, 1, -1, -1),
        (1, -2, -1, 2),
        (1, -1, -1, 1),
        (1, -alpha, 2, -2 * alpha),
    ]:
        targets = [
            dr.define_einst(
                r[a, b], c1 * p[a, c] * x[c, b] + c2 * p[a, c] * y[c, b]
            ),
            dr.define_einst(
                s[a, b], c3 * q[a, c] * x[c, b] + c4 * q[a, c] * y[c, b]
            ),
        ]

        eval_seq = optimize(targets)

        assert verify_eval_seq(eval_seq, targets)
        assert len(eval_seq) == 3


def test_removal_of_shallow_interms(simple_drudge):
    """Test if removal of shallow intermediates can be turned on/off."""

    dr = simple_drudge

    r = dr.r
    a, b, c, d = dr.ds[:4]

    x = IndexedBase("x")
    y = IndexedBase("y")
    z = IndexedBase("z")
    u = IndexedBase("u")

    targets = [
        dr.define(
            u,
            (a, r),
            (b, r),
            (c, r),
            dr.sum((d, r), x[a, d] * y[b, d] * z[c, d]),
        )
    ]

    for i in [True, False]:
        eval_seq = optimize(targets, remove_shallow=i)
        verify_eval_seq(eval_seq, targets)
        assert len(eval_seq) == (1 if i else 2)
        continue


def test_get_cost_on_zero_cost(simple_drudge):
    """Test correct behaviour of get_flop_cost at input with no FLOP cost."""

    dr = simple_drudge

    a, b = dr.ds[:2]

    x = IndexedBase("x")
    r = IndexedBase("y")

    targets = [dr.define_einst(x[a, b], r[a, b])]

    for i in [get_flop_cost(targets), get_flop_cost(targets, leading=True)]:
        assert i == 0
        continue


@pytest.mark.parametrize(
    "ext_names,sum_name",
    [
        ("a b", "c"),
        ("a b", "e"),
        ("b c", "d"),
        ("c d", "e"),
        ("f g", "h"),
        ("g h", "a"),
    ],
)
def test_verification_of_a_result_with_any_external_symbols(
    simple_drudge, ext_names, sum_name
):
    """Verification must not depend on which symbols the caller used.

    The optimization canonicalizes the external indices of a result onto the
    leading dummies of the range, so they need not be the symbols that were
    written.  The verification used to subtract the two definitions as they
    stood, read that renaming as a difference, and reject a correct answer
    for every choice of external symbols but the canonical one.

    The last case matters most: the caller's external symbols are the ones
    the optimization picked as a summation dummy inside the result, so
    lining the two up naively would capture that summation.
    """

    dr = simple_drudge
    i0, i1 = symbols(ext_names)
    k = symbols(sum_name)

    x = IndexedBase("X")
    y, z, u = (IndexedBase(i) for i in ["Y", "Z", "U"])
    res = IndexedBase("res")

    targets = [
        dr.define_einst(
            res[i0, i1], x[i0, k] * y[k, i1] + x[i0, k] * z[k, i1] + u[i0, i1]
        )
    ]

    eval_seq = optimize(targets)
    assert verify_eval_seq(eval_seq, targets)


def test_verification_still_rejects_a_wrong_result(simple_drudge):
    """The check has to keep failing on an answer that is actually wrong.

    Guards the fix above from being a check that passes everything, in both
    the canonical spelling of the external indices and another one.
    """

    dr = simple_drudge

    x = IndexedBase("X")
    y, z, u = (IndexedBase(i) for i in ["Y", "Z", "U"])
    res = IndexedBase("res")

    for ext_names, sum_name in [("a b", "c"), ("c d", "e")]:
        i0, i1 = symbols(ext_names)
        k = symbols(sum_name)
        targets = [
            dr.define_einst(
                res[i0, i1],
                x[i0, k] * y[k, i1] + x[i0, k] * z[k, i1] + u[i0, i1],
            )
        ]

        eval_seq = list(optimize(targets))
        assert verify_eval_seq(eval_seq, targets)

        # Scaling the result makes it wrong, whatever the indices are called.
        final = eval_seq[-1]
        eval_seq[-1] = TensorDef(final.base, final.exts, final.rhs * 2)
        with pytest.raises(ValueError):
            verify_eval_seq(eval_seq, targets)
