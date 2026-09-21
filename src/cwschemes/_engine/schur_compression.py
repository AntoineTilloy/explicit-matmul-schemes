"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from fractions import Fraction as F
from itertools import product
from math import prod
from functools import lru_cache
from bisect import bisect_right
from .restrictions import completion, audit
from .patches import plan_group, group_cost
TABLES = {}

@lru_cache(None)
def _h_patterns(fa, fb, k):
    pairs = {(a ^ b, a & b) for a in fa for b in fb}
    return [h for h in range(1 << k) if any((not sym & ~h and (not common & h) for sym, common in pairs))]

def h_patterns(fa, fb, k):
    return _h_patterns(tuple(sorted(fa)), tuple(sorted(fb)), k)

def compress_mask(mask, positions):
    return sum(((mask >> i & 1) << j for j, i in enumerate(positions)))

def group_candidate(rows, m, s, pair, q):
    positions = [i for i in range(3 * m) if not s >> i & 1]
    families = [sorted({compress_mask(r[h], positions) for r in rows if r[h] & s == s}) for h in range(3)]
    hs = h_patterns(families[pair[0]], families[pair[1]], len(positions))
    return dict(pure=s, positions=positions, pair=pair, families=families, allowed_h=hs, rank=sum((q ** h.bit_count() for h in hs)), raw_rank=(q + 1) ** len(positions))

def compression_plan(rows, m, q):
    assert q >= 1
    w = audit(rows, m)
    groups = []
    for s in w['surviving_pure_subsets']:
        patched = plan_group(rows, m, s)
        pc = group_cost(q, patched)
        candidates = [group_candidate(rows, m, s, p, q) for p in [(0, 1), (0, 2), (1, 2)]]
        best = min(candidates, key=lambda c: c['rank'])
        groups.append(dict(pure=s, patched_rank=pc, compressed_rank=best['rank'], compression=best, use_compression=best['rank'] < pc))
    return dict(rank=sum((min(g['patched_rank'], g['compressed_rank']) for g in groups)), patched_rank=sum((g['patched_rank'] for g in groups)), groups=groups)

class CompressedGroup:

    def __init__(self, q, rows, m, pure=0, pair=(0, 1), cache=None):
        assert q >= 1
        self.q = q
        self.m = m
        self.pair = tuple(pair)
        self.third = 3 - sum(pair)
        self.spec = group_candidate(rows, m, pure, pair, q)
        self.positions = self.spec['positions']
        self.k = len(self.positions)
        self.allowed = self.spec['allowed_h']
        self.allowed_set = set(self.allowed)
        self.rank = self.spec['rank']
        self.families = [set(f) for f in self.spec['families']]
        self.ends = []
        offset = 0
        for h in self.allowed:
            offset += q ** h.bit_count()
            self.ends.append(offset)
        self.forms, self.weights = completion(q)
        self.tables = {} if cache is None else cache

    def term(self, index):
        assert 0 <= index < self.rank
        j = bisect_right(self.ends, index)
        h = self.allowed[j]
        if j:
            index -= self.ends[j - 1]
        result = [self.q] * self.k
        for i in reversed(range(self.k)):
            if h >> i & 1:
                index, result[i] = divmod(index, self.q)
        return result

    def orbit_setup(self, zero_mask):
        q = self.q
        sizes = [2 if zero_mask >> i & 1 or q == 1 else 3 for i in range(self.k)]
        states = list(product(*(range(s) for s in sizes)))
        removed = [i for i, v in enumerate(states) if sum(((x != 0) << j for j, x in enumerate(v))) not in self.allowed_set]
        mats = [[[4 * q, -3 * q], [-3, 3]] if zero_mask >> i & 1 else [[4, -3], [-3, 3]] if q == 1 else [[4 * q, -3, -3 * (q - 1)], [-3, 3, 0], [-3, 0, 3]] for i in range(self.k)]
        gs = [[-4 * q * q, 3 * q] if zero_mask >> i & 1 else [-6, 6] if q == 1 else [-6 * q, 3 * q + 3, 3] for i in range(self.k)]
        return (sizes, states, removed, mats, gs)

    def build_table(self, zero_mask):
        import sys
        import json, hashlib
        from pathlib import Path
        key = json.dumps([1, self.q, self.k, self.allowed, zero_mask], separators=(',', ':'))
        digest = hashlib.sha256(key.encode()).hexdigest()
        folder = resource_root() / 'coefficient_cache/schur'
        path = folder / (digest + '.json')
        if digest in TABLES:
            self.tables[zero_mask] = TABLES[digest]
            return self.tables[zero_mask]
        if path.exists():
            data = json.loads(path.read_text())
            assert data['key'] == key
            answer = [F(a, b) for a, b in data['values']]
            TABLES[digest] = answer
            self.tables[zero_mask] = answer
            return answer
        from sympy.polys.matrices import DomainMatrix
        from sympy import ZZ
        sizes, states, removed, mats, gs = self.orbit_setup(zero_mask)
        n = len(removed)
        a = [[prod((mats[j][states[r][j]][states[c][j]] for j in range(self.k))) for c in removed] for r in removed]
        b = [[prod((gs[j][states[r][j]] for j in range(self.k)))] for r in removed]
        if n:
            am = DomainMatrix.from_list(a, ZZ)
            bm = DomainMatrix.from_list(b, ZZ)
            numerator, denominator = am.solve_den(bm)
            values = numerator.to_Matrix()
            den = int(denominator) * self.q ** self.k
            x = [F(int(values[i, 0]), den) for i in range(n)]
        else:
            x = []
        vector = [F(0)] * len(states)
        for r, value in zip(removed, x):
            vector[r] = value
        for axis, mat in enumerate(mats):
            stride = prod(sizes[axis + 1:])
            block = stride * sizes[axis]
            out = [F(0)] * len(states)
            for start in range(0, len(states), block):
                for suffix in range(stride):
                    for i in range(sizes[axis]):
                        out[start + i * stride + suffix] = sum((F(mat[i][j], 3) * vector[start + j * stride + suffix] for j in range(sizes[axis])))
            vector = out
        scale = (3 * self.q) ** self.k
        answer = [F(prod((gs[j][v[j]] for j in range(self.k))), scale) - vector[i] for i, v in enumerate(states)]
        assert all((answer[i] == 0 for i in removed))
        self.tables[zero_mask] = answer
        TABLES[digest] = answer
        folder.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(dict(key=key, values=[(x.numerator, x.denominator) for x in answer]), separators=(',', ':')) + '\n')
        temporary.replace(path)
        return answer

    def factor_entry(self, mode, index, variable):
        assert len(variable) == 3 * self.m
        s = self.spec['pure']
        if any((variable[i] != 0 for i in range(3 * self.m) if s >> i & 1)):
            return F(0)
        v = [variable[i] for i in self.positions]
        mask = sum(((x == 0) << i for i, x in enumerate(v)))
        if mask not in self.families[mode]:
            return F(0)
        term = self.term(index)
        if mode != self.third:
            return prod((self.forms[t][x] for t, x in zip(term, v)))
        if mask not in self.tables:
            self.build_table(mask)
        address = 0
        for t, x in zip(term, v):
            width = 2 if x == 0 or self.q == 1 else 3
            state = 0 if t == self.q else 1 if x == 0 or t == x - 1 else 2
            address = address * width + state
        return self.tables[mask][address] * self.weights[-1] ** s.bit_count()

class SchurPacket:
    """Use exact Schur compression whenever it beats the existing group patch."""
    compression = 'pair-schur-v1'

    def __init__(self, q, rows, m):
        from .scheme import CompactPacket
        self.q = q
        self.rows = rows
        self.m = m
        self.patched = True
        self.base = CompactPacket(q, rows, m, patched=True)
        self.compression_plan = compression_plan(rows, m, q)
        self.plan = dict(self.base.plan, rank=self.compression_plan['rank'])
        self.rank = self.plan['rank']
        self.groups = []
        self._ends = []
        total = 0
        for j, g in enumerate(self.compression_plan['groups']):
            if g['use_compression']:
                spec = g['compression']
                obj = CompressedGroup(q, rows, m, g['pure'], spec['pair'])
                self.groups.append(obj)
                total += obj.rank
            else:
                self.groups.append(None)
                total += g['patched_rank']
            self._ends.append(total)
        assert total == self.rank

    def matrix_variable(self, *args):
        return self.base.matrix_variable(*args)

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2) or not 0 <= index < self.rank:
            raise IndexError((mode, index))
        if len(variable) != 3 * self.m or any((not 0 <= x <= self.q for x in variable)):
            raise ValueError('invalid variable')
        gi = bisect_right(self._ends, index)
        local = index - (self._ends[gi - 1] if gi else 0)
        g = self.groups[gi]
        if g is not None:
            return g.factor_entry(mode, local, variable)
        start = self.base._ends[gi - 1] if gi else 0
        return self.base.factor_entry(mode, start + local, variable)
