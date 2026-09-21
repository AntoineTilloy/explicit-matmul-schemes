"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from fractions import Fraction as F
from bisect import bisect_right
from math import prod
from pathlib import Path
from functools import lru_cache
import hashlib, json
from .schur_compression import CompressedGroup
from .scheme import CompactPacket
ROOT = resource_root()
CACHE = {}

@lru_cache(32)
def proof_data(path, mtime_ns, size):
    """Share immutable proof data among repeated copies of the same packet."""
    raw = Path(path).read_bytes()
    d = json.loads(raw)
    steps = []
    for step in d['steps']:
        steps.append(dict(free=3 - sum(step['pair']), kept=frozenset(step['kept']), blocks={int(I): dict(basis=b['basis'], deleted=b['deleted'], coefficients=[[F(c) for c in row] for row in b['coefficients']]) for I, b in step['blocks'].items()}))
    return (raw, d, steps)

class IteratedOrbitPacket:
    compression = 'orbit-iterated-v1'

    def __init__(self, q, rows, m, relations):
        self.q = q
        self.rows = rows
        self.m = m
        self.patched = True
        self.relations_path = relations
        path = (ROOT / relations).resolve()
        stat = path.stat()
        self.proof_bytes, self.proof, self.steps = proof_data(str(path), stat.st_mtime_ns, stat.st_size)
        d = self.proof
        assert d['format'] == 'cw-iterated-orbit-relations-v1' and d['complete'] and (d['q'] == q) and (d['m'] == m) and (q >= 2)
        assert sorted(map(tuple, d['rows'])) == sorted(map(tuple, rows))
        self.pure = d.get('pure', 0)
        self.group = CompressedGroup(q, rows, m, self.pure)
        self.k = self.group.k
        self.base = CompactPacket(q, rows, m, patched=True)
        assert self.pure == 0 and self.base.plan['witness']['surviving_pure_subsets'] == [0]
        self.kept = d['steps'][-1]['kept']
        self._ends = []
        total = 0
        for H in self.kept:
            total += q ** H.bit_count()
            self._ends.append(total)
        self.rank = total
        assert total == d['rank']
        self.plan = dict(self.base.plan, rank=total)
        self.tables = {}

    def matrix_variable(self, *args):
        return self.base.matrix_variable(*args)

    def term(self, index):
        j = bisect_right(self._ends, index)
        H = self.kept[j]
        local = index - (self._ends[j - 1] if j else 0)
        digits = [self.q] * self.k
        for i in reversed(range(self.k)):
            if H >> i & 1:
                local, digits[i] = divmod(local, self.q)
        return digits

    def build_table(self, mode, zero):
        key = (mode, zero)
        digest = hashlib.sha256(self.proof_bytes + json.dumps(key).encode()).hexdigest()
        if digest in CACHE:
            self.tables[key] = CACHE[digest]
            return CACHE[digest]
        folder = ROOT / 'results/round6/iterated_cache'
        path = folder / (digest + '.json')
        if path.exists():
            data = json.loads(path.read_text())
            assert data['digest'] == digest
            table = [F(a, b) for a, b in data['values']]
            CACHE[digest] = table
            self.tables[key] = table
            return table
        q = self.q
        k = self.k
        values = {}
        for I in range(1 << k):
            if I & zero:
                continue
            v = {}
            for H in range(1 << k):
                if H & I != I:
                    continue
                factors = []
                for i in range(k):
                    if I >> i & 1:
                        continue
                    if zero >> i & 1:
                        factors.append(F(-4 * q, 3) if mode == 2 and (not H >> i & 1) else F(1))
                    elif H >> i & 1:
                        factors.append(F(2, q))
                    else:
                        factors.append(F(-2) if mode == 2 else F(3, 2 * q))
                v[H] = prod(factors)
            values[I] = v
        for step in self.steps:
            updated = {I: {H: v for H, v in values[I].items() if H in step['kept']} for I in values}
            if step['free'] == mode:
                for I, b in step['blocks'].items():
                    if I not in values:
                        continue
                    inputs = [values[I].get(H, F(0)) for H in b['deleted']]
                    for H, row in zip(b['basis'], b['coefficients']):
                        updated[I][H] += sum((c * v for c, v in zip(row, inputs) if c and v))
            values = updated
        sizes = [2 if zero >> i & 1 else 3 for i in range(k)]
        table = [F(0)] * prod(sizes)
        for I, vals in values.items():
            for H, value in vals.items():
                address = 0
                for i, s in enumerate(sizes):
                    address = address * s + (2 if I >> i & 1 else 1 if H >> i & 1 else 0)
                table[address] = value
        for axis, s in enumerate(sizes):
            if s == 2:
                continue
            stride = prod(sizes[axis + 1:])
            block = 3 * stride
            out = table.copy()
            for start in range(0, len(table), block):
                for tail in range(stride):
                    average = table[start + stride + tail]
                    standard = table[start + 2 * stride + tail]
                    out[start + stride + tail] = average + F(q - 1, q) * standard
                    out[start + 2 * stride + tail] = average - F(1, q) * standard
            table = out
        folder.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(dict(digest=digest, mode=mode, zero=zero, values=[(v.numerator, v.denominator) for v in table]), separators=(',', ':')) + '\n')
        tmp.replace(path)
        CACHE[digest] = table
        self.tables[key] = table
        return table

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2) or not 0 <= index < self.rank:
            raise IndexError((mode, index))
        if len(variable) != 3 * self.m or any((not 0 <= x <= self.q for x in variable)):
            raise ValueError('invalid variable')
        if any((variable[i] != 0 for i in range(3 * self.m) if self.pure >> i & 1)):
            return F(0)
        variable = [variable[i] for i in self.group.positions]
        zero = sum(((x == 0) << i for i, x in enumerate(variable)))
        if zero not in self.group.families[mode]:
            return F(0)
        key = (mode, zero)
        if key not in self.tables:
            self.build_table(mode, zero)
        address = 0
        for t, x in zip(self.term(index), variable):
            width = 2 if x == 0 else 3
            state = 0 if t == self.q else 1 if x == 0 or t == x - 1 else 2
            address = address * width + state
        return self.tables[key][address] * (F(self.q, 3) ** self.pure.bit_count() if mode == 2 else 1)
