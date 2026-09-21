"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from fractions import Fraction as F
from bisect import bisect_right
from math import prod
from pathlib import Path
import hashlib, json
from .schur_compression import CompressedGroup
from .scheme import CompactPacket
ROOT = resource_root()
CACHE = {}

class SecondSchurGroup:
    compression = 'pair-schur-second-v1'

    def __init__(self, q, rows, m, relations):
        self.q = q
        self.rows = rows
        self.m = m
        self.patched = True
        self.relations_path = relations
        path = ROOT / relations
        self.relations_bytes = path.read_bytes()
        self.proof = json.loads(self.relations_bytes)
        d = self.proof
        assert d['format'] == 'cw-second-orbit-relations-v1' and d['q'] == q and (d['m'] == m)
        assert sorted(map(tuple, d['rows'])) == sorted(map(tuple, rows))
        self.pure = d.get('pure', 0)
        self.first = CompressedGroup(q, rows, m, self.pure)
        assert d['allowed_first'] == self.first.allowed
        self.base = CompactPacket(q, rows, m, patched=True)
        self.kept = d['kept']
        self.kept_set = set(self.kept)
        self.ends = []
        offset = 0
        for H in self.kept:
            offset += q ** H.bit_count()
            self.ends.append(offset)
        self.rank = offset
        assert self.rank == d['rank']
        self._ends = self.ends
        self.plan = dict(self.base.plan, rank=self.rank)
        self.first_starts = {H: self.first.ends[j - 1] if j else 0 for j, H in enumerate(self.first.allowed)}
        self.tables = {}

    def matrix_variable(self, *args):
        return self.base.matrix_variable(*args)

    def term(self, index):
        j = bisect_right(self.ends, index)
        H = self.kept[j]
        local = index - (self.ends[j - 1] if j else 0)
        return (H, local, self.first_starts[H] + local)

    def build_middle_table(self, zero):
        q = self.q
        k = self.first.k
        d = self.proof
        digest = hashlib.sha256(self.relations_bytes + str(zero).encode()).hexdigest()
        if digest in CACHE:
            self.tables[zero] = CACHE[digest]
            return CACHE[digest]
        folder = ROOT / 'results/round6/second_cache'
        path = folder / (digest + '.json')
        if path.exists():
            data = json.loads(path.read_text())
            assert data['digest'] == digest
            table = [F(a, b) for a, b in data['values']]
            CACHE[digest] = table
            self.tables[zero] = table
            return table
        sizes = [2 if zero >> i & 1 else 3 for i in range(k)]
        table = [F(0)] * prod(sizes)

        def raw(I, H):
            return prod((F(1) if zero >> i & 1 else F(2, q) if H >> i & 1 else F(3, 2 * q) for i in range(k) if not I >> i & 1))
        for I in range(1 << k):
            if I & zero:
                continue
            values = {H: raw(I, H) for H in self.kept if I & H == I}
            block = d['blocks'].get(str(I))
            if block:
                inputs = [raw(I, H) for H in block['deleted']]
                for H, row in zip(block['basis'], block['coefficients']):
                    values[H] += sum((F(c) * v for c, v in zip(row, inputs)))
            for H, value in values.items():
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
        tmp.write_text(json.dumps(dict(digest=digest, values=[(v.numerator, v.denominator) for v in table]), separators=(',', ':')) + '\n')
        tmp.replace(path)
        CACHE[digest] = table
        self.tables[zero] = table
        return table

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2) or not 0 <= index < self.rank:
            raise IndexError((mode, index))
        if len(variable) != 3 * self.m or any((not 0 <= x <= self.q for x in variable)):
            raise ValueError('invalid variable')
        H, local, first_index = self.term(index)
        if mode != 1:
            return self.first.factor_entry(mode, first_index, variable)
        if any((variable[i] != 0 for i in range(3 * self.m) if self.pure >> i & 1)):
            return F(0)
        variable = [variable[i] for i in self.first.positions]
        zero = sum(((x == 0) << i for i, x in enumerate(variable)))
        if zero not in self.first.families[1]:
            return F(0)
        if zero not in self.tables:
            self.build_middle_table(zero)
        term = self.first.term(first_index)
        address = 0
        for t, x in zip(term, variable):
            s = 2 if x == 0 else 3
            state = 0 if t == self.q else 1 if x == 0 or t == x - 1 else 2
            address = address * s + state
        return self.tables[zero][address]

class SecondSchurPacket(SecondSchurGroup):

    def __init__(self, q, rows, m, relations):
        super().__init__(q, rows, m, relations)
        assert self.pure == 0 and self.base.plan['witness']['surviving_pure_subsets'] == [0]

class MultiSchurPacket:
    """Choose the cheapest certified group compression or existing patch."""
    compression = 'pair-schur-multi-v1'

    def __init__(self, q, rows, m, relations):
        from .schur_compression import SchurPacket
        self.q = q
        self.rows = rows
        self.m = m
        self.patched = True
        self.relations_path = relations
        self.base = SchurPacket(q, rows, m)
        self.groups = []
        self._ends = []
        total = 0
        self.relation_map = {json.loads((ROOT / path).read_text()).get('pure', 0): path for path in relations}
        for i, s in enumerate(self.base.plan['witness']['surviving_pure_subsets']):
            cost = self.base._ends[i] - (self.base._ends[i - 1] if i else 0)
            group = SecondSchurGroup(q, rows, m, self.relation_map[s]) if s in self.relation_map else None
            if group is not None and group.rank < cost:
                self.groups.append(group)
                total += group.rank
            else:
                self.groups.append(None)
                total += cost
            self._ends.append(total)
        self.rank = total
        self.plan = dict(self.base.plan, rank=total)

    def matrix_variable(self, *args):
        return self.base.matrix_variable(*args)

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2) or not 0 <= index < self.rank:
            raise IndexError((mode, index))
        i = bisect_right(self._ends, index)
        t = index - (self._ends[i - 1] if i else 0)
        if self.groups[i] is not None:
            return self.groups[i].factor_entry(mode, t, variable)
        return self.base.factor_entry(mode, t + (self.base._ends[i - 1] if i else 0), variable)
