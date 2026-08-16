"""Tests of the parenthesization extension at its own interface.

The optimizer feeds every product node to ``parenth``.  Most of its behaviour
is covered through the optimizer, but the degenerate shapes are easier to pin
down here.
"""

import pytest

from gristmill._parenth import parenth


@pytest.mark.parametrize('mode', [0, 1, 2])
@pytest.mark.parametrize('if_incl', [False, True])
def test_factor_without_indices(mode, if_incl):
    """A single factor with no index at all gives the trivial memoir.

    No dimensions, no summations, one factor involving nothing.  This is what
    the optimizer hands in for an index-free product, and it used to bring the
    interpreter down: the empty dimension list gave an exhausted iterator whose
    null value was copied while the ``dims`` vector was built (issue #45).
    """

    res = parenth([], 0, [[]], mode, if_incl)

    assert list(res.keys()) == [(0,)]
    interm = res[(0,)]
    assert interm.sums == ()
    assert interm.exts == ()
    assert len(interm.evals) == 1
    ev = interm.evals[0]
    assert ev.ops == ((0,), ())
    assert ev.sums == ()
    assert ev.cost == 0


def test_factor_with_a_single_index():
    """The smallest non-degenerate problem, for comparison."""

    res = parenth([5], 1, [[0]], 1, True)

    assert list(res.keys()) == [(0,)]
    assert res[(0,)].sums == (0,)
    assert res[(0,)].exts == ()
