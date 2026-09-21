"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from .tailored_completion import TailoredPacket

class FixedSectorQuotient:

    def __init__(self, q, rows, m, permutation, characters=(1, 1, 1), label_flips=(), label_fixed=None):
        self.q = q
        self.rows = list(map(tuple, rows))
        self.m = m
        self.k = 3 * m
        self.n = q ** m
        self.perm = tuple(permutation)
        self.characters = tuple(characters)
        self.label_flips = tuple(label_flips)
        self.label_fixed = q % 2 if label_fixed is None else label_fixed
        assert 0 <= self.label_fixed <= q and (q - self.label_fixed) % 2 == 0
        assert all((self.perm[i] == i for i in self.label_flips))
        assert len(self.perm) == self.k and all((self.perm[self.perm[i]] == i for i in range(self.k)))
        assert characters[0] * characters[1] * characters[2] == 1
        self.base = TailoredPacket(q, self.rows, m)
        self.spec = self.base.spec
        assert self.spec['common'] == 0
        assert all((self.spec['missing'][i] == self.spec['missing'][self.perm[i]] for i in range(self.k)))
        rowids = {r: i for i, r in enumerate(self.rows)}
        self.products = []
        self.free = []
        self.fixed = []
        for h, row in enumerate(self.rows):
            g = rowids[tuple((self.mask(z) for z in row))]
            if h < g:
                self.free.append(h)
                self.products.append(dict(kind='free', row=h, dims=(self.n,) * 3))
            elif h == g:
                self.fixed.append(h)
        for h in self.fixed:
            masks = (self.rows[h][1], self.rows[h][2], self.rows[h][0])
            for si in (1, -1):
                signs = (si, characters[0] * si, characters[0] * characters[1] * si)
                dims = tuple((len(self.eigenbasis(mask)[0 if sign == 1 else 1]) for mask, sign in zip(masks, signs)))
                if all(dims):
                    self.products.append(dict(kind='fixed', row=h, dims=dims, signs=signs, masks=masks))
        digits = self.spec['digits']
        code = {tuple(ds): h for h, ds in digits.items()}
        self.blocks = []
        self.ends = []
        rank = fixed_terms = 0
        for H in self.spec['kept']:
            ds = digits[H]
            g = code[self.permute(ds)]
            O = [i for i, t in enumerate(ds) if t == 1]
            if H > g:
                continue
            if H < g:
                self.blocks.append((H, None, O))
                rank += q ** len(O)
                self.ends.append(rank)
                continue
            components = []
            for i in O:
                if i < self.perm[i]:
                    components.append(('pair', i, self.perm[i], q, q * (q - 1) // 2, q * q))
                elif i == self.perm[i]:
                    components.append(('flip' if i in self.label_flips else 'single', i, i, self.label_fixed if i in self.label_flips else q, (q - self.label_fixed) // 2 if i in self.label_flips else 0, q))
            fixed_count = 1
            for x in components:
                fixed_count *= x[3]
            fixed_terms += fixed_count
            for pivot, x in enumerate(components):
                radices = [y[3] if j < pivot else y[4] if j == pivot else y[5] for j, y in enumerate(components)]
                size = 1
                for r in radices:
                    size *= r
                if size:
                    self.blocks.append((H, (components, pivot), radices))
                    rank += size
                    self.ends.append(rank)
            if characters == (1, 1, 1) and fixed_count:
                self.blocks.append((H, (components, None), [x[3] for x in components]))
                rank += fixed_count
                self.ends.append(rank)
        self.fixed_terms = fixed_terms
        self.rank = rank
        assert rank == (self.base.rank + (fixed_terms if characters == (1, 1, 1) else -fixed_terms)) // 2

    def label_image(self, x):
        return x if x == 0 or x > self.q - self.label_fixed else x + 1 if x % 2 else x - 1

    def variable_image(self, v):
        out = list(self.permute(v))
        for i in self.label_flips:
            out[i] = self.label_image(out[i])
        return tuple(out)

    def mask(self, z):
        return sum((1 << self.perm[i] for i in range(self.k) if z >> i & 1))

    def permute(self, v):
        out = [None] * self.k
        for i, x in enumerate(v):
            out[self.perm[i]] = x
        return tuple(out)

    @lru_cache(None)
    def eigenbasis(self, mask):
        positions = [i for i in range(self.k) if mask >> i & 1]
        assert self.mask(mask) == mask

        def image(t):
            v = [0] * self.k
            for i in reversed(positions):
                t, x = divmod(t, self.q)
                v[i] = x + 1
            v = self.variable_image(v)
            r = 0
            for i in positions:
                r = r * self.q + v[i] - 1
            return r
        fixed = []
        pairs = []
        for i in range(self.n):
            j = image(i)
            if i == j:
                fixed.append(((i, F(1)),))
            elif i < j:
                pairs.append((i, j))
        plus = fixed + [tuple(((i, F(1)) for i in pair)) for pair in pairs]
        minus = [((i, F(1)), (j, F(-1))) for i, j in pairs]
        return (plus, minus)

    def eigenvector(self, mask, sign, index, dual=False):
        v = self.eigenbasis(mask)[0 if sign == 1 else 1][index]
        return tuple(((i, c / 2 if dual and len(v) == 2 else c) for i, c in v))

    @lru_cache(512)
    def term(self, t):
        if not 0 <= t < self.rank:
            raise IndexError(t)
        b = bisect_right(self.ends, t)
        local = t - (self.ends[b - 1] if b else 0)
        H, data, reps = self.blocks[b]
        labels = [0] * self.k
        multiplicity = 2
        if data is None:
            for i in reversed(reps):
                local, x = divmod(local, self.q)
                labels[i] = x + 1
        else:
            components, pivot = data
            values = [0] * len(reps)
            for j in reversed(range(len(reps))):
                local, values[j] = divmod(local, reps[j])
            for j, (kind, i, k, fc, nc, ac) in enumerate(components):
                x = values[j]
                state = 'fixed' if pivot is None or j < pivot else 'pivot' if j == pivot else 'all'
                if kind == 'single':
                    labels[i] = x + 1
                elif kind == 'flip':
                    labels[i] = self.q - self.label_fixed + x + 1 if state == 'fixed' else 2 * x + 1 if state == 'pivot' else x + 1
                elif state == 'fixed':
                    labels[i] = labels[k] = x + 1
                elif state == 'all':
                    u, v = divmod(x, self.q)
                    labels[i] = u + 1
                    labels[k] = v + 1
                else:
                    u = 0
                    while x >= self.q - u - 1:
                        x -= self.q - u - 1
                        u += 1
                    labels[i] = u + 1
                    labels[k] = u + 2 + x
            if pivot is None:
                multiplicity = 1
        return (H, tuple(labels), multiplicity)

    def original(self, mode, H, labels, variable):
        Z = sum(((x == 0) << i for i, x in enumerate(variable)))
        if Z not in self.spec['families'][mode]:
            return F(0)
        value = 1
        for i, (t, x) in enumerate(zip(self.spec['digits'][H], variable)):
            if t == 1:
                if x and x != labels[i]:
                    return F(0)
            else:
                if bool(x) != (mode == self.spec['missing'][i]):
                    return F(0)
                if mode == 2:
                    value = -value
        return F(value)

    def factor_entry(self, mode, t, product, row, col):
        p = self.products[product]
        I, K, J = p['dims']
        shape = [(I, K), (K, J), (I, J)][mode]
        if mode not in (0, 1, 2) or not (0 <= row < shape[0] and 0 <= col < shape[1]):
            raise IndexError((mode, product, row, col))
        H, labels, multiplicity = self.term(t)
        h = p['row']
        if p['kind'] == 'free':
            v = self.base.matrix_variable(h, mode, row, col)
            value = self.original(mode, H, labels, v) + self.characters[mode] * self.original(mode, H, labels, self.variable_image(v))
            if mode == 2:
                value /= 2
        else:
            a, b = [(0, 1), (1, 2), (0, 2)][mode]
            masks = p['masks']
            signs = p['signs']
            left = self.eigenvector(masks[a], signs[a], row, dual=mode == 2)
            right = self.eigenvector(masks[b], signs[b], col, dual=mode != 2)
            value = F(0)
            for i, x in left:
                for j, y in right:
                    value += x * y * self.original(mode, H, labels, self.base.matrix_variable(h, mode, i, j))
        if mode == 2:
            value *= multiplicity
        return value
