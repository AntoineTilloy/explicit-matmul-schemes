"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
from .restrictions import audit
from .scheme import CompactPacket

@lru_cache(128)
def layout(rows, m, force_symmetric=(), allow_cross=False):
    if not allow_cross:
        audit(rows, m)
    else:
        full = (1 << 3 * m) - 1
        assert rows and all((all((Z.bit_count() == m and 0 <= Z <= full for Z in r)) and r[0] | r[1] | r[2] == full and (not (r[0] & r[1] or r[0] & r[2] or r[1] & r[2])) for r in rows))
        assert all((len({r[h] for r in rows}) == len(rows) for h in range(3)))
    k = 3 * m
    full = (1 << k) - 1
    families = [sorted({r[h] for r in rows}) for h in range(3)]
    unions = []
    for family in families:
        union = 0
        for Z in family:
            union |= Z
        unions.append(union)
    common = unions[0] & unions[1] & unions[2] | sum((1 << i for i in force_symmetric))
    missing = [next((h for h in range(3) if not unions[h] >> i & 1), None) for i in range(k)]
    digits = {}
    ordinary = {}
    requirements = [{} for _ in range(3)]
    kept = []
    for ds in product(*[(0, 1, 2) if common >> i & 1 else (1, 2) for i in range(k)]):
        H = sum((t * 3 ** (k - 1 - i) for i, t in enumerate(ds)))
        P = sum(((t == 0) << i for i, t in enumerate(ds)))
        O = sum(((t == 1) << i for i, t in enumerate(ds)))
        S = full ^ (P | O)
        req = [P | sum(((t == 2 and (not common >> i & 1) and (missing[i] is not None) and (missing[i] != h)) << i for i, t in enumerate(ds))) for h in range(3)]
        if all((any((Z & r == r for Z in family)) for r, family in zip(req, families))):
            kept.append(H)
            digits[H] = ds
            ordinary[H] = O
            for h in range(3):
                requirements[h][H] = req[h]
    return dict(k=k, full=full, families=families, common=common, missing=missing, digits=digits, ordinary=ordinary, requirements=requirements, kept=kept)

def coefficient(q, spec, mode, Z, I, H, scalar=F, p=None):
    if Z & I or spec['requirements'][mode][H] & Z != spec['requirements'][mode][H]:
        return scalar(0)
    O = spec['ordinary'][H]
    if O & I != I:
        return scalar(0)
    ds = spec['digits'][H]
    C = spec['common']
    full = spec['full']
    P = sum(((t == 0) << i for i, t in enumerate(ds)))
    S = full ^ (P | O)

    def frac(a, b):
        return a * pow(b, -1, p) % p if p else scalar(a, b)

    def power(x, n):
        return pow(x, n, p) if p else x ** n
    value = power(frac(2, q), (O & C & ~Z & ~I).bit_count()) * power(frac(1, q), (O & ~C & ~Z & ~I).bit_count())
    if mode == 2:
        value *= power(scalar(-2), (S & C & ~Z).bit_count()) * power(frac(-4 * q, 3), (S & C & Z).bit_count()) * power(frac(q, 3), P.bit_count()) * power(scalar(-1), (S & ~C).bit_count())
    else:
        value *= power(frac(3, 2 * q), (S & C & ~Z).bit_count())
    return value % p if p else value

def make_state(q, rows, m, exact=False, force_symmetric=(), allow_cross=False):
    if exact:
        from .rational_state import ExactState
        from .rational_arrays import Q
        s = ExactState.__new__(ExactState)
        scalar = Q
        p = None
        dtype = object
    else:
        raise ValueError('Only exact certificate replay is supported')
    spec = layout(tuple(map(tuple, rows)), m, tuple(force_symmetric), allow_cross)
    s.q = q
    s.rows = rows
    s.m = m
    s.k = 3 * m
    s.families = spec['families']
    s.digits = spec['digits']
    s.ordinary = spec['ordinary']
    s.kept = spec['kept'].copy()
    s.reindex()
    s.coefs = [[] for _ in range(3)]
    for mode in range(3):
        for I, ps in enumerate(s.patterns):
            s.coefs[mode].append(np.array([[coefficient(q, spec, mode, Z, I, H, scalar, p) for H in ps] for Z in s.families[mode]], dtype=dtype))
    return s

class TailoredPacket:
    compression = 'tailored-completion-v1'

    def __init__(self, q, rows, m):
        self.q = q
        self.rows = rows
        self.m = m
        self.k = 3 * m
        self.patched = True
        self.spec = layout(tuple(map(tuple, rows)), m)
        self.base = CompactPacket(q, rows, m, patched=True)
        self.kept = self.spec['kept']
        self.ends = []
        total = 0
        for H in self.kept:
            total += q ** self.spec['ordinary'][H].bit_count()
            self.ends.append(total)
        self.rank = total
        self.plan = dict(self.base.plan, rank=total)

    def matrix_variable(self, *args):
        return self.base.matrix_variable(*args)

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2) or not 0 <= index < self.rank:
            raise IndexError((mode, index))
        if len(variable) != self.k or any((not 0 <= x <= self.q for x in variable)):
            raise ValueError('invalid variable')
        Z = sum(((x == 0) << i for i, x in enumerate(variable)))
        if Z not in self.spec['families'][mode]:
            return F(0)
        j = bisect_right(self.ends, index)
        H = self.kept[j]
        local = index - (self.ends[j - 1] if j else 0)
        ds = self.spec['digits'][H]
        labels = {}
        for i in reversed(range(self.k)):
            if ds[i] == 1:
                local, labels[i] = divmod(local, self.q)
        value = F(1)
        for i, (t, x) in enumerate(zip(ds, variable)):
            if self.spec['common'] >> i & 1:
                if t == 0:
                    value *= F(int(x == 0)) * (F(self.q, 3) if mode == 2 else 1)
                elif t == 1:
                    value *= F(1) if x == 0 else F(int(labels[i] == x - 1)) + F(1, self.q)
                else:
                    value *= (F(1) if x == 0 else F(3, 2 * self.q)) * (F(-4 * self.q, 3) if mode == 2 else 1)
            elif t == 1:
                value *= F(int(x == 0 or labels[i] == x - 1))
            else:
                value *= F(int(x != 0 if mode == self.spec['missing'][i] else x == 0)) * (-1 if mode == 2 else 1)
            if not value:
                return value
        return value
