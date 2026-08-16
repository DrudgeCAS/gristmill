"""Tests that the optimization result does not depend on incidental order.

The constriction graph gets its vertices numbered in the order edges arrive,
which follows the order the libparenth memoir is walked.  Issue #43 was the
search reading those numbers, so that a different but equally valid walk gave
a different answer.  These tests permute the walk and require the same result
every time.  On the code before the fix, the ascending walk lost the effective
T of the CCSD energy equation and both matrix factorizations below.
"""

import random

import pytest
from drudge import Drudge, PartHoleDrudge, Range
from sympy import symbols, Symbol, IndexedBase, Rational

import gristmill.optimize
from gristmill import optimize, verify_eval_seq, get_flop_cost

_REAL_SORTED = sorted


def _looks_like_memoir(items):
    """If the items are the entries of a libparenth memoir.

    The memoir walk is the only sort in the optimizer that passes reverse=True
    over pairs keyed by tuples of integers.
    """
    if not items:
        return False
    for entry in items:
        if not (isinstance(entry, tuple) and len(entry) == 2):
            return False
        key = entry[0]
        if not (
                isinstance(key, tuple) and key
                and all(isinstance(i, int) for i in key)
        ):
            return False
    return True


def _walk_memoir(monkeypatch, order):
    """Make the optimizer walk the memoir in the given order.

    ``order`` maps the sorted list of memoir entries to a permutation of it.
    Every other call to ``sorted`` inside the optimizer is left alone.
    """

    def patched(iterable, **kwargs):
        items = list(iterable)
        if kwargs.get('reverse') and _looks_like_memoir(items):
            return order(_REAL_SORTED(items))
        return _REAL_SORTED(items, **kwargs)

    monkeypatch.setattr(gristmill.optimize, 'sorted', patched, raising=False)


def _shuffled(seed):
    def order(items):
        items = list(items)
        random.Random(seed).shuffle(items)
        return items

    return order


ORDERS = {
    'descending': lambda items: list(reversed(items)),
    'ascending': lambda items: list(items),
    'by_size': lambda items: _REAL_SORTED(
        items, key=lambda x: (len(x[0]), x[0])
    ),
    'shuffle_1': _shuffled(1),
    'shuffle_2': _shuffled(2),
    'shuffle_3': _shuffled(3),
}


@pytest.fixture(scope='module')
def problems(spark_ctx):
    """The problems whose answer used to depend on the walk order.

    Each entry gives the targets, the substitutions for the optimizer, and the
    number of definitions the right answer has.
    """

    res = {}

    # The CCSD energy equation, where the effective T has to be found.  Not
    # the CCSD energy equation exactly, the same as in the CC tests.
    dr = PartHoleDrudge(spark_ctx)
    p = dr.names
    a, b = p.V_dumms[:2]
    i, j = p.O_dumms[:2]
    u = dr.two_body
    t = IndexedBase('t')
    energy = dr.define_einst(
        Symbol('e'),
        u[i, j, a, b] * t[a, b, i, j] * Rational(1, 2)
        + u[i, j, a, b] * t[a, i] * t[b, j]
    )
    res['ccsd_energy'] = ([energy], {p.nv: p.no * 10}, 2)

    # Two matrix problems.
    dr = Drudge(spark_ctx)
    m = symbols('m')
    r = Range('M', 0, m)
    dumms = symbols('a b c d e f g h')
    dr.set_dumms(r, dumms)
    dr.add_default_resolver(r)
    a, b, c, e = dumms[0], dumms[1], dumms[2], dumms[4]
    u, x, y, z = (IndexedBase(i) for i in ['U', 'X', 'Y', 'Z'])

    target = dr.define_einst(
        IndexedBase('T')[a, b],
        u[a, b] * z[c, e] * x[e, c] + u[a, b] * z[c, e] * y[e, c]
    )
    res['disconnected_outer_product'] = ([target], {}, 3)

    target = dr.define_einst(
        Symbol('T'), x[b, a] * z[a, b] + y[a, b] * z[b, a]
    )
    res['needing_canonicalization'] = ([target], {}, 2)

    return res


@pytest.fixture(scope='module')
def reference(problems):
    """The result of the production walk, for comparison."""

    res = {}
    for name, (targets, substs, n_defs) in problems.items():
        eval_seq = optimize(targets, substs=substs)
        assert verify_eval_seq(eval_seq, targets)
        assert len(eval_seq) == n_defs
        res[name] = get_flop_cost(eval_seq)
    return res


@pytest.mark.parametrize('order', list(ORDERS))
def test_memoir_walk_order_does_not_change_the_result(
        monkeypatch, problems, reference, order
):
    """The result must be the same under any walk of the memoir."""

    _walk_memoir(monkeypatch, ORDERS[order])

    for name, (targets, substs, n_defs) in problems.items():
        eval_seq = optimize(targets, substs=substs)
        assert verify_eval_seq(eval_seq, targets)
        assert len(eval_seq) == n_defs, name
        assert get_flop_cost(eval_seq) == reference[name], name


def test_rand_constr_gives_random_constrictions(spark_ctx):
    """Random constriction has to be random.

    An earlier version of the search sorted the candidates after shuffling
    them, which silently turned ``rand_constr`` off.  Here a sum with many
    equally good constrictions is optimized under a few seeds, and more than
    one answer has to come out.  Every answer must still be correct.
    """

    dr = Drudge(spark_ctx)
    n = symbols('n')
    r = Range('r', 0, n)
    dr.set_dumms(r, symbols('a b c d'))
    dr.add_default_resolver(r)
    a, b = symbols('a b')

    xs = [IndexedBase('x%d' % i) for i in range(3)]
    ys = [IndexedBase('y%d' % i) for i in range(3)]
    amp = sum(
        (i + 2 * j + 1) * xs[i][a, b] * ys[j][a, b]
        for i in range(3) for j in range(3)
    )
    targets = [dr.define_einst(Symbol('r'), amp)]

    seen = set()
    for seed in range(8):
        random.seed(seed)
        eval_seq = optimize(targets, rand_constr=True)
        assert verify_eval_seq(eval_seq, targets)
        seen.add(tuple(str(i) for i in eval_seq))

    assert len(seen) > 1
