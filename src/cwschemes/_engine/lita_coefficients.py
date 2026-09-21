"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from functools import lru_cache
from itertools import permutations
from fractions import Fraction
from math import comb
from .sources import lita_round6 as L

class CoordinateCoefficient(L._Coordinates):

    def __init__(self, n, variable):
        super().__init__(n)
        self.row, self.col = divmod(variable, n)

    def entry(self, block, i, j):
        D = self.D
        br, bc = divmod(block, 2)
        r, c = (self.row, self.col)
        sign = -1 if c >= D else 1
        if i < D and j < D:
            value = self.grid * sign if r == br * D + i and c == bc * D + j else 0
        elif i < D:
            value = self.grid * sign // (2 * (D - 2)) if r == br * D + i else 0
        elif j < D:
            value = self.grid * sign // (2 * (D - 2)) if c == bc * D + j else 0
        else:
            value = self.grid * sign // (4 * (D - 2) ** 2)
        return L._Form({0: value}) if value else L._Form()

@lru_cache(6)
def context(n, variable):
    c = CoordinateCoefficient(n, variable)
    return (c, tuple((L._View(c, a) for a in range(3))))

def combination_at(n, k, index):
    result = []
    first = 0
    for left in range(k, 0, -1):
        for x in range(first, n - left + 1):
            count = comb(n - x - 1, left - 1)
            if index < count:
                result.append(x)
                first = x + 1
                break
            index -= count
        else:
            raise IndexError(index)
    return tuple(result)

def boundary(c, views, i):
    D, w = (c.D, c.w)
    f = tuple((L._cycle(v, i, i, D, 0) for v in views))
    g = tuple((L._cycle(v, D, D, i, 0) for v in views))
    yield (f[0], f[1], w * f[2] + w * w * g[2])
    f = tuple((L._cycle(v, i, i, D, 1) for v in views))
    g = tuple((L._cycle(v, D, D, i, 1) for v in views))
    yield (f[0] + w * g[0], f[1] + w * g[1], -w * f[2] / (1 + w))
    for ids in ((i, i, D), (i, D, i), (i, D, D), (D, i, i), (D, i, D), (D, D, i)):
        for barred in (0, 1):
            yield L._weighted(L._mixed_factors(views, *ids, barred), c, ids)
    for left, right in ((i, D), (D, i)):
        for u, v, z in L._heptad(tuple((L._seeds(view, left, right) for view in views))):
            yield (c.weight(left) * u, c.weight(right) * v, 2 * z)

def term(c, views, index):
    D = c.D
    triangles = 16 * comb(D + 1, 3)
    edges = 24 * comb(D, 2)
    if index < triangles:
        a, t = divmod(index, 16)
        vertices = combination_at(D + 1, 3, a)
        for i, j, k in permutations(vertices):
            for barred in (0, 1):
                if i < j < k or k < j < i:
                    if t == 0:
                        forms = [L._cycle(v, i, j, k, barred) for v in views]
                        forms[2] = (-1 if barred else 1) * forms[2]
                        return L._weighted(forms, c, (i, j, k))
                    t -= 1
                if t == 0:
                    return L._weighted(L._mixed_factors(views, i, j, k, barred), c, (i, j, k))
                t -= 1
        raise AssertionError(index)
    index -= triangles
    if index < edges:
        a, t = divmod(index, 24)
        i, j = combination_at(D, 2, a)
        r, shift = divmod(t, 3)
        recipes = tuple((L._edge_recipe(L._edge_forms(v, i, j))[r] for v in views))
        return tuple((recipes[a][(a + shift) % 3] for a in range(3)))
    index -= edges
    if index < 28 * D:
        i, t = divmod(index, 28)
        return next((forms for j, forms in enumerate(boundary(c, views, i)) if j == t))
    index -= 28 * D
    return next((forms for j, forms in enumerate(L._center_terms(c, views)) if j == index))

def entry(n, mode, index, row, col):
    rank = L.lita_rank(n)
    if mode not in (0, 1, 2) or not 0 <= index < rank or (not (0 <= row < n and 0 <= col < n)):
        raise IndexError((mode, index, row, col))
    variable = col * n + row if mode == 2 else row * n + col
    c, views = context(n, variable)
    return Fraction(term(c, views, index)[mode].get(0, 0), c.grid)
