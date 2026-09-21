"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
import itertools, json
from pathlib import Path
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from .fixed_sector_quotient import FixedSectorQuotient
from .general_square import Ordinary
from .square_scheme import STRASSEN, base_factors
ROOT = resource_root()
AXES = ((0, 1), (1, 2), (0, 2))

@lru_cache(None)
def small_rect(dims):
    if dims == (4, 4, 4):
        return (48, base_factors(4))
    if len(set(dims)) == 1 and dims[0] in (2, 3, 5, 6, 7):
        n = dims[0]
        s = Ordinary(dict(kind='lille', n=n, rank={2: 7, 3: 23, 5: 93, 6: 153, 7: 250}[n]))
        return (s.rank, [[[s.factor_entry(h, t, i, j) for t in range(s.rank)] for i in range(n) for j in range(n)] for h in range(3)])
    for shape in [(3, 3, 4), (3, 4, 4)]:
        if sorted(dims) != list(shape):
            continue
        d = json.loads((ROOT / 'seeds' / ('rect' + 'x'.join(map(str, shape)) + '.json')).read_text())
        perm = next((p for p in itertools.permutations(range(3)) if all((shape[i] == dims[p[i]] for i in range(3)))))
        factors = []
        for mode, (a, b) in enumerate(AXES):
            old = next((h for h, ab in enumerate(AXES) if {perm[x] for x in ab} == {a, b}))
            u, v = AXES[old]
            rows = []
            for i in range(dims[a]):
                for j in range(dims[b]):
                    coords = {a: i, b: j}
                    rows.append(d['factors'][old][coords[perm[u]] * shape[v] + coords[perm[v]]])
            factors.append(rows)
        return (d['rank'], factors)
    I, K, J = dims
    rank = I * K * J
    factors = [[[0] * rank for _ in range(dims[a] * dims[b])] for a, b in AXES]
    for i in range(I):
        for k in range(K):
            for j in range(J):
                t = (i * K + k) * J + j
                for mode, var in enumerate((i * K + k, k * J + j, i * J + j)):
                    factors[mode][var][t] = 1
    return (rank, factors)

class RectOrdinary:

    def __init__(self, dims, leaf, recipe=None, allow_split=True):
        self.square = None
        self.dims = tuple(dims)
        self.leaf = leaf
        self.h = leaf.n
        assert all((x % self.h == 0 for x in dims))
        self.small = tuple((x // self.h for x in dims))
        self.outer_cube = None
        self.factors = None
        if recipe:
            self.kind = recipe['kind']
            self.rank = recipe['rank']
            if self.kind == 'split':
                self.pieces = [(tuple(p['offset']), RectOrdinary(p['dims'], leaf, p['recipe'])) for p in recipe['pieces']]
                self.piece_ends = []
                total = 0
                for _, child in self.pieces:
                    total += child.rank
                    self.piece_ends.append(total)
                assert total == self.rank
                return
            if self.kind == 'square':
                self.square = Ordinary(recipe['square_recipe'])
                assert self.rank == self.square.rank
                return
            if recipe.get('outer_cube_recipe'):
                self.outer_cube = Ordinary(recipe['outer_cube_recipe'])
                self.outer_rank = self.outer_cube.rank
            elif sorted(self.small) in ([3, 3, 4], [3, 4, 4]):
                self.outer_rank, self.factors = small_rect(self.small)
            else:
                self.outer_rank = self.small[0] * self.small[1] * self.small[2]
            assert self.rank == self.outer_rank * leaf.rank
            return
        raise ValueError('An explicit rectangular recipe is required')

    def certificate(self):
        out = dict(dims=self.dims, kind=self.kind, square_recipe=self.square.recipe if self.square else None, outer_cube_recipe=self.outer_cube.recipe if self.outer_cube else None, leaf_recipe=self.leaf.recipe, rank=self.rank)
        if self.kind == 'split':
            out['pieces'] = [dict(offset=offset, dims=c.dims, recipe=c.certificate()) for offset, c in self.pieces]
        return out

    def factor_entry(self, mode, t, row, col):
        a, b = AXES[mode]
        if not (0 <= row < self.dims[a] and 0 <= col < self.dims[b]):
            return F(0)
        if self.kind == 'split':
            which = bisect_right(self.piece_ends, t)
            offset, child = self.pieces[which]
            local = t - (self.piece_ends[which - 1] if which else 0)
            return child.factor_entry(mode, local, row - offset[a], col - offset[b])
        if self.square:
            return self.square.factor_entry(mode, t, row, col)
        ot, it = divmod(t, self.leaf.rank)
        br, ir = divmod(row, self.h)
        bc, ic = divmod(col, self.h)
        if self.outer_cube:
            value = self.outer_cube.factor_entry(mode, ot, br, bc)
        elif self.factors:
            value = F(self.factors[mode][br * self.small[b] + bc][ot])
        else:
            I, K, J = self.small
            i, rest = divmod(ot, K * J)
            k, j = divmod(rest, J)
            value = F((br, bc) == [(i, k), (k, j), (i, j)][mode])
        return value * self.leaf.factor_entry(mode, it, ir, ic)

def unequal_maps(a, b):
    M = [[], [], []]

    def add(mode, call, *rects):
        while len(M[mode]) <= call:
            M[mode].append([])
        M[mode][call] = list(rects)
    for mode in (0, 1):
        add(mode, 0, (0, 0, a, a, 1), (a, a, b, b, 1))
    add(2, 0, (0, 0, a, a, 1), (a, a, b, b, 1))
    add(0, 1, (a, 0, b, a, 1), (a, a, b, b, 1))
    add(1, 1, (0, 0, a, a, 1))
    add(2, 1, (a, 0, b, a, 1), (a, a, b, b, -1))
    add(0, 2, (0, 0, a, a, 1))
    add(1, 2, (0, a, a, b, 1), (a, a, b, b, -1))
    add(2, 2, (0, a, a, b, 1), (a, a, b, b, 1))
    add(0, 3, (a, a, b, b, 1))
    add(1, 3, (a, 0, b, a, 1), (0, 0, b, a, -1))
    add(2, 3, (0, 0, b, a, 1), (a, 0, b, a, 1))
    add(0, 4, (0, a, a, b, 1), (0, 0, a, b, 1))
    add(1, 4, (a, a, b, b, 1))
    add(2, 4, (0, 0, a, b, -1), (0, a, a, b, 1))
    add(0, 5, (a, 0, b, a, 1), (0, 0, b, a, -1))
    add(1, 5, (0, 0, a, b, 1), (0, a, a, b, 1))
    add(2, 5, (a, a, b, b, 1))
    add(0, 6, (0, a, a, b, 1), (a, a, b, b, -1))
    add(1, 6, (a, 0, b, a, 1), (a, a, b, b, 1))
    add(2, 6, (0, 0, a, a, 1))
    return (M, [(a, a, a), (b, a, a), (a, a, b), (b, b, a), (a, b, b), (b, a, b), (a, b, a)])

class FixedQuotientSquare:
    format = 'cw-fixed-quotient-square-v1'

    def __init__(self, packet, a, b, leaf_recipe, target=None, ordinary_recipes=None):
        self.packet = packet
        self.a = a
        self.b = b
        self.inner = a + b
        assert packet.n == self.inner and len(packet.free) == 6
        self.full_n = 2 * self.inner
        self.n = target or self.full_n
        assert self.n <= self.full_n
        self.leaf_recipe = leaf_recipe
        self.leaf = Ordinary(leaf_recipe)
        self.maps, self.call_dims = unequal_maps(a, b)
        c = a - b
        assert c >= 0
        whole = lambda call: (call, (0, 0, 0), self.call_dims[call])
        if packet.characters == (-1, -1, 1):
            self.assigned = [whole(6), whole(5), whole(3), (2, (0, 0, 0), (b, a, b))]
            self.remaining = [whole(0), whole(1), whole(4), (2, (b, 0, 0), (c, a, b))]
        elif packet.characters == (1, 1, 1):
            self.assigned = [whole(0), (3, (0, 0, 0), (b, b, b)), whole(6), (4, (0, 0, 0), (b, b, b))]
            self.remaining = [whole(1), whole(2), whole(5), (3, (0, 0, b), (b, b, c)), (4, (b, 0, 0), (c, b, b))]
        else:
            raise ValueError('Unsupported square assignment')
        assert len(packet.products) == 10
        for p, part in zip(packet.products[6:], self.assigned):
            assert all((x <= y for x, y in zip(part[2], p['dims'])))
        self.remaining = [part for part in self.remaining if all(part[2])]
        self.ordinary = [RectOrdinary(part[2], self.leaf, ordinary_recipes[i] if ordinary_recipes else None) for i, part in enumerate(self.remaining)]
        self.ends = [packet.rank]
        for s in self.ordinary:
            self.ends.append(self.ends[-1] + s.rank)
        self.rank = self.ends[-1]

    def pullback(self, mode, part, row, col):
        call, offset, dims = part
        a, b = AXES[mode]
        out = []
        for r, c, h, w, sign in self.maps[mode][call]:
            if r <= row < r + h and c <= col < c + w:
                i, j = (row - r - offset[a], col - c - offset[b])
                if 0 <= i < dims[a] and 0 <= j < dims[b]:
                    out.append((i, j, sign))
        return out

    def factor_entry(self, mode, t, row, col):
        if mode not in (0, 1, 2) or not (0 <= t < self.rank and 0 <= row < self.n and (0 <= col < self.n)):
            raise IndexError((mode, t, row, col))
        br, ir = divmod(row, self.inner)
        bc, ic = divmod(col, self.inner)
        outer = STRASSEN[mode][2 * br + bc]
        value = F(0)
        if t < self.packet.rank:
            for h in range(6):
                if outer[h]:
                    value += outer[h] * self.packet.factor_entry(mode, t, h, ir, ic)
            if outer[6]:
                for h, part in enumerate(self.assigned, 6):
                    for i, j, sign in self.pullback(mode, part, ir, ic):
                        value += outer[6] * sign * self.packet.factor_entry(mode, t, h, i, j)
        elif outer[6]:
            which = bisect_right(self.ends, t) - 1
            local = t - self.ends[which]
            for i, j, sign in self.pullback(mode, self.remaining[which], ir, ic):
                value += outer[6] * sign * self.ordinary[which].factor_entry(mode, local, i, j)
        return value

    def certificate(self):
        p = self.packet
        return dict(format=self.format, coefficient_field='Q', n=self.n, full_n=self.full_n, rank=self.rank, q=p.q, m=p.m, rows=p.rows, permutation=p.perm, characters=p.characters, label_flips=p.label_flips, label_fixed=p.label_fixed, a=self.a, b=self.b, leaf_recipe=self.leaf_recipe, packet_rank=p.rank, fixed_term_count=p.fixed_terms, remaining_ranks=[s.rank for s in self.ordinary], ordinary_recipes=[s.certificate() for s in self.ordinary], coefficient_generator='fixed_quotient_square.load_certificate')

def load_certificate(value):
    d = json.loads(Path(value).read_text()) if isinstance(value, (str, Path)) else value
    assert d['format'] == FixedQuotientSquare.format
    p = FixedSectorQuotient(d['q'], d['rows'], d['m'], d['permutation'], tuple(d['characters']), d.get('label_flips', ()), d.get('label_fixed'))
    s = FixedQuotientSquare(p, d['a'], d['b'], d['leaf_recipe'], d['n'], d.get('ordinary_recipes'))
    assert (s.rank, p.rank, p.fixed_terms) == (d['rank'], d['packet_rank'], d['fixed_term_count'])
    assert [x.rank for x in s.ordinary] == d['remaining_ranks']
    return s
