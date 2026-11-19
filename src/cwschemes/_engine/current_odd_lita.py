"""Indexed access to the pinned rational odd-LITA construction for all odd N >= 9.

Only O(N^2) prototype groups are indexed rather than materializing O(N^3) products.
"""
from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import comb
from pathlib import Path

from .sources import lita_odd_current as source


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


def _make_coordinates(n):
    if type(n) is not int or n < 9 or n % 2 != 1:
        raise ValueError("N must be an odd integer >= 9")
    old_check = source._check_dimension
    source._check_dimension = lambda _: None
    try:
        return source._Coordinates(n)
    finally:
        source._check_dimension = old_check


def odd_lita_rank(n):
    if type(n) is not int or n < 9 or n % 2 != 1:
        raise ValueError("N must be an odd integer >= 9")
    m = (n - 3) // 2
    return 16 * comb(m, 3) + 68 * comb(m, 2) + 98 * m + 49


def _triangle(c, views, vertices):
    for i, j, k in permutations(vertices):
        for barred in (0, 1):
            if i < j < k or k < j < i:
                row = tuple(source._cycle(v, i, j, k, barred) for v in views)
                row = row[0], row[1], (-1 if barred else 1) * row[2]
                yield source._weighted(row, c, (i, j, k))
            yield source._weighted(source._mixed_factors(views, i, j, k, barred), c, (i, j, k))


def _edge(c, views, ids):
    i, j = ids
    for left, right in ((i, j), (j, i)):
        recipes = tuple(source._edge_recipe(v, left, right) for v in views)
        for r in range(4):
            yield from source._cyclic(tuple(recipe[r] for recipe in recipes))


def _vertex(c, views, i):
    D, w = c.D, c.w
    recipes = []
    for v in views:
        Y, g = v.entry, v.gamma
        u, ell = Y(1, i, i) + g, Y(2, i, i) + g
        k0 = v.half_sum(0, i, i) + w * (Y(1, i, D) + Y(2, D, i)) + (D - 4) * g - c.defects[0][i]
        k3 = v.half_sum(3, i, i) + w * (Y(1, D, i) + Y(2, i, D)) + (D - 4) * g - c.defects[3][i]
        recipes.append(((2 * u, ell, k0), (-2 * ell, u, k3)))
    for r in range(2):
        yield from source._cyclic(tuple(recipe[r] for recipe in recipes))


def _boundary(c, views, i):
    D = c.D
    for barred in (0, 1):
        b = 3 * barred
        term = tuple(v.entry(b, i, D) + v.entry(b, D, i) + ((D - 2) * v.entry(b, D, D) - v.gamma) / (D - 3) for v in views)
        yield (D - 2) * (D - 3) * term[0], term[1], (-1 if barred else 1) * term[2]
    for ids in ((i, i, D), (i, D, i), (i, D, D), (D, i, i), (D, i, D), (D, D, i)):
        for barred in (0, 1):
            yield source._weighted(source._mixed_factors(views, *ids, barred), c, ids)
    for left, right in ((i, D), (D, i)):
        for u, v, z in source._heptad(tuple(source._seeds(view, left, right) for view in views)):
            yield c.weight(left) * u, c.weight(right) * v, 2 * z


def _col_corr(c, views, i):
    for j in range(c.D + 1):
        if i == j:
            continue
        recipes = []
        for v in views:
            Y = v.entry
            recipes.append((
                (c.defects[0][i], Y(0, i, j) + Y(1, i, j), Y(2, i, j) - Y(2, i, i)),
                (-c.defects[3][i], Y(2, i, j) + Y(3, i, j), Y(1, i, j) - Y(1, i, i)),
            ))
        for r in range(2):
            for u, v, z in source._cyclic(tuple(recipe[r] for recipe in recipes)):
                yield c.weight(i) * u, c.weight(j) * v, z


class CurrentOddLITA:
    """Return exact U/V/W terms, with W in output row/column convention."""
    def __init__(self, n, source_sha256):
        actual = sha256(Path(source.__file__).read_bytes()).hexdigest()
        if source_sha256 != actual:
            raise ValueError("current odd LITA source SHA256 mismatch")
        self.n = n
        self.rank = odd_lita_rank(n)
        self.source_sha256 = actual
        self._context = None

    def _initialize(self):
        if self._context is not None:
            return
        c = _make_coordinates(self.n)
        view = source._View(c)
        views = (view, view, view)
        self._context = c, views
        self._triangles = 16 * comb(c.D + 1, 3)
        self._groups = []
        self._ends = [self._triangles]
        families = (
            (_edge, combinations(range(c.D), 2)),
            (_vertex, range(c.D)),
            (_boundary, range(c.D)),
            (source._center_terms, (None,)),
            (_col_corr, range(c.D + 1)),
        )
        for function, keys in families:
            for key in keys:
                terms = function(c, views) if key is None else function(c, views, key)
                slots = tuple(i for i, forms in enumerate(terms) if all(forms))
                if slots:
                    self._groups.append((function, key, slots))
                    self._ends.append(self._ends[-1] + len(slots))
        if self._ends[-1] != self.rank:
            self._context = None
            raise ArithmeticError("indexed current odd LITA rank mismatch")

    def term(self, t):
        if type(t) is not int or not 0 <= t < self.rank:
            raise IndexError(t)
        self._initialize()
        c, views = self._context
        if t < self._triangles:
            group, slot = divmod(t, 16)
            terms = _triangle(c, views, _combination_at(c.D + 1, 3, group))
        else:
            group = bisect_right(self._ends, t) - 1
            function, key, slots = self._groups[group]
            slot = slots[t - self._ends[group]]
            terms = function(c, views) if key is None else function(c, views, key)
        forms = next(forms for index, forms in enumerate(terms) if index == slot)
        result = []
        for axis, form in enumerate(forms):
            row = {}
            for index, numerator in form.items():
                i, j = divmod(index, self.n)
                row[(j, i) if axis == 2 else (i, j)] = Fraction(numerator, c.grid)
            result.append(row)
        return tuple(result)
