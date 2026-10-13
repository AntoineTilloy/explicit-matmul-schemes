"""Indexed access to the pinned current rational LITA construction.

Only O(n²) prototype labels are indexed, rather than materializing O(n³)
products. Coefficient forms are created lazily on the first term request.
"""
from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import comb
from pathlib import Path

from .sources import lita_current as source


def _combination_at(n, size, index):
    result = []
    low = 0
    for left in range(size - 1, -1, -1):
        for value in range(low, n):
            count = comb(n - value - 1, left)
            if index < count:
                result.append(value)
                low = value + 1
                break
            index -= count
        else:
            raise IndexError(index)
    return tuple(result)


def _triangle(co, views, vertices):
    # Same ordering, anchor merge and weights as the pinned generator.
    for i, j, k in permutations(vertices):
        for barred in (0, 1):
            if i < j < k or k < j < i:
                term = [source._cycle(v, i, j, k, barred) for v in views]
                term[2] *= -1 if barred else 1
                yield source._weighted(term, co, (i, j, k))
            if barred == 1 and i == 1 and k == 0 and 2 <= j < co.D:
                continue
            term = source._mixed_factors(views, i, j, k, barred)
            if barred == 0 and i == 0 and j == 1 and 2 <= k < co.D:
                other = source._mixed_factors(views, 1, k, 0, 1)
                if term[:2] != other[:2]:
                    raise ArithmeticError("anchor triangle alignment failed")
                term = term[0], term[1], term[2] + other[2]
            yield source._weighted(term, co, (i, j, k))


def _edge(co, views, ids):
    i, j = ids
    for left, right in ((i, j), (j, i)):
        recipes = tuple(source._edge_recipe(v, left, right) for v in views)
        for r in range(4):
            yield from source._cyclic(tuple(recipe[r] for recipe in recipes))


def _vertex(co, views, i):
    D, w = co.D, co.w
    recipes = []
    for v in views:
        Y, g = v.entry, v.gamma
        u, l = Y(1, i, i) + g, Y(2, i, i) + g
        k0 = v.half_sum(0, i, i) + w*(Y(1, i, D)+Y(2, D, i)) + (D-4)*g
        k3 = v.half_sum(3, i, i) + w*(Y(1, D, i)+Y(2, i, D)) + (D-4)*g
        recipes.append(((2*u, l, k0), (-2*l, u, k3)))
    for r in range(2):
        yield from source._cyclic(tuple(recipe[r] for recipe in recipes))


def _boundary(co, views, i):
    D = co.D
    for barred in (0, 1):
        b = 3*barred
        term = tuple(v.entry(b, i, D)+v.entry(b, D, i)
                     + ((D-2)*v.entry(b, D, D)-v.gamma)/(D-3) for v in views)
        yield (D-2)*(D-3)*term[0], term[1], (-1 if barred else 1)*term[2]
    for ids in ((i, i, D), (i, D, i), (i, D, D), (D, i, i), (D, i, D), (D, D, i)):
        for barred in (0, 1):
            yield source._weighted(source._mixed_factors(views, *ids, barred), co, ids)
    for left, right in ((i, D), (D, i)):
        for u, v, z in source._heptad(tuple(source._seeds(view, left, right) for view in views)):
            yield co.weight(left)*u, co.weight(right)*v, 2*z


class CurrentLITA:
    """Return exact U/V/W terms, with W in output row/column convention."""
    def __init__(self, n, source_sha256):
        source._check_dimension(n)
        actual = sha256(Path(source.__file__).read_bytes()).hexdigest()
        if source_sha256 != actual:
            raise ValueError("current LITA source SHA256 mismatch")
        self.n = n
        self.rank = source.lita_rank(n)
        self.source_sha256 = actual
        self._context = None

    def _initialize(self):
        if self._context is not None:
            return
        co = source._Coordinates(self.n)
        views = tuple(source._View(co, axis) for axis in range(3))
        self._context = co, views
        # Triangle products have no additional zero terms; the D-2 anchor
        # groups contain 15 rather than 16 products. All other families are
        # filtered on the full exact forms, just as lita() does.
        self._triangles = 16*comb(co.D+1, 3)-(co.D-2)
        self._groups = []
        self._ends = [self._triangles]
        families = ((_edge, combinations(range(co.D), 2)),
                    (_vertex, range(co.D)), (_boundary, range(co.D)),
                    (source._center_terms, (None,)))
        for function, keys in families:
            for key in keys:
                terms = function(co, views) if key is None else function(co, views, key)
                slots = tuple(i for i, forms in enumerate(terms) if all(forms))
                if slots:
                    self._groups.append((function, key, slots))
                    self._ends.append(self._ends[-1]+len(slots))
        if self._ends[-1] != self.rank:
            self._context = None
            raise ArithmeticError("indexed current LITA rank mismatch")

    def term(self, t):
        if type(t) is not int or not 0 <= t < self.rank:
            raise IndexError(t)
        self._initialize()
        co, views = self._context
        if t < self._triangles:
            anchor = 15*(co.D-2)
            if t < anchor:
                group, slot = divmod(t, 15)
            else:
                group, slot = divmod(t-anchor, 16)
                group += co.D-2
            terms = _triangle(co, views, _combination_at(co.D+1, 3, group))
        else:
            group = bisect_right(self._ends, t)-1
            function, key, slots = self._groups[group]
            slot = slots[t-self._ends[group]]
            terms = function(co, views) if key is None else function(co, views, key)
        forms = next(forms for index, forms in enumerate(terms) if index == slot)
        result = []
        for axis, form in enumerate(forms):
            row = {}
            for index, numerator in form.items():
                i, j = divmod(index, self.n)
                row[(j, i) if axis == 2 else (i, j)] = Fraction(numerator, co.grid)
            result.append(row)
        return tuple(result)
