"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib, json
from .tailored_completion import layout, TailoredPacket
from .compact_orbit_packet import proof_data, CompactOrbitPacket
ROOT = resource_root()

@lru_cache(64)
def metadata(rows, m, coordinates):
    spec = layout(rows, m)
    digits = spec['digits'].copy()
    ordinary = spec['ordinary'].copy()
    kept = spec['kept'].copy()
    for c in coordinates:
        bit = 1 << c
        nextid = max(digits) + 1
        new = []
        for H in kept:
            new.append(H)
            if ordinary[H] & bit:
                J = nextid
                nextid += 1
                ds = list(digits[H])
                ds[c] = 3
                digits[J] = tuple(ds)
                ordinary[J] = ordinary[H] ^ bit
                new.append(J)
        kept = sorted(new)
    return (spec, digits, ordinary, kept)

def initial_coefficient(q, spec, digits, coordinates, mode, Z, E, I, H):
    if (Z | E) & I:
        return F(0)
    value = F(1)
    split = set(coordinates)
    for i, t in enumerate(digits[H]):
        zero = bool(Z >> i & 1)
        exceptional = bool(E >> i & 1)
        common = bool(spec['common'] >> i & 1)
        if I >> i & 1:
            if t != 1:
                return F(0)
            continue
        if common:
            if t == 0:
                factor = F(int(zero)) * (F(q, 3) if mode == 2 else 1)
            elif t in (1, 3):
                factor = F(1) if zero else F(1, q) + (F(int(exceptional)) if t == 3 else F(0) if exceptional else F(1, q - 1) if i in split else F(1, q))
            else:
                factor = (F(1) if zero else F(3, 2 * q)) * (F(-4 * q, 3) if mode == 2 else 1)
        elif t in (1, 3):
            factor = F(1) if zero else F(int(exceptional)) if t == 3 else F(0) if exceptional else F(1, q - 1) if i in split else F(1, q)
        else:
            factor = F(int(not zero if mode == spec['missing'][i] else zero)) * (-1 if mode == 2 else 1)
        value *= factor
        if not value:
            break
    return value
TABLES = {}

class RefinedOrbitPacket:
    compression = 'refined-orbit-v1'

    def __init__(self, q, rows, m, relations):
        self.q = q
        self.rows = rows
        self.m = m
        self.k = 3 * m
        self.patched = True
        self.relations_path = relations
        path = (ROOT / relations).resolve()
        stat = path.stat()
        self.proof_bytes, self.proof, self.steps = proof_data(str(path), stat.st_mtime_ns, stat.st_size)
        d = self.proof
        assert d['format'] == 'cw-refined-orbit-relations-v1' and d['complete'] and (d['q'] == q) and (d['m'] == m) and (sorted(map(tuple, rows)) == sorted(map(tuple, d['rows'])))
        self.coordinates = tuple(d['coordinates'])
        self.qs = [q - (i in self.coordinates) for i in range(self.k)]
        self.joint_base = 'base_proof' in d
        if self.joint_base:
            self.base = CompactOrbitPacket(q, rows, m, d['base_proof'])
            self.spec = dict(families=[sorted(f) for f in self.base.families])
            self.digits = dict(enumerate(self.base.digits))
            self.ordinary = dict(enumerate(self.base.ordinary))
            self.initial_kept = self.base.kept.copy()
            self.origins = {H: H for H in self.initial_kept}
            for c in self.coordinates:
                new = []
                nextid = max(self.digits) + 1
                for H in self.initial_kept:
                    new.append(H)
                    if self.ordinary[H] >> c & 1:
                        N = nextid
                        nextid += 1
                        ds = list(self.digits[H])
                        ds[c] = 3
                        self.digits[N] = tuple(ds)
                        self.ordinary[N] = self.ordinary[H] ^ 1 << c
                        self.origins[N] = self.origins[H]
                        new.append(N)
                self.initial_kept = sorted(new)
        else:
            self.spec, self.digits, self.ordinary, self.initial_kept = metadata(tuple(map(tuple, rows)), m, self.coordinates)
            self.base = TailoredPacket(q, rows, m)
        self.kept = d['steps'][-1]['kept']
        self.ends = []
        total = 0
        for H in self.kept:
            size = 1
            for i in range(self.k):
                if self.ordinary[H] >> i & 1:
                    size *= self.qs[i]
            total += size
            self.ends.append(total)
        self.rank = total
        assert total == d['rank']
        self.plan = dict(self.base.plan, rank=total)
        self.tables = {}

    def matrix_variable(self, *args):
        return self.base.matrix_variable(*args)

    def term(self, index):
        if not 0 <= index < self.rank:
            raise IndexError(index)
        j = bisect_right(self.ends, index)
        H = self.kept[j]
        local = index - (self.ends[j - 1] if j else 0)
        labels = [1 if t == 3 else 0 for t in self.digits[H]]
        for i in reversed(range(self.k)):
            if self.ordinary[H] >> i & 1:
                local, x = divmod(local, self.qs[i])
                labels[i] = x + 1 + (i in self.coordinates)
        return (H, labels)

    def build_table(self, mode, Z, E):
        key = (mode, Z, E)
        digest = hashlib.sha256(self.proof_bytes + json.dumps(key).encode()).hexdigest()
        if digest in TABLES:
            self.tables[key] = TABLES[digest]
            return TABLES[digest]
        values = {I: {} for I in range(1 << self.k) if not I & (Z | E)}
        if self.joint_base:
            base = self.base.build_table(mode, Z)
            split = sum((1 << i for i in self.coordinates))
            for H in self.initial_kept:
                for J, v in base[self.origins[H]].items():
                    refined = J & split
                    fixed = J & ~split
                    if fixed & (Z | E):
                        continue
                    terms = {fixed: v}
                    for c in self.coordinates:
                        if not refined >> c & 1:
                            continue
                        exceptional_term = self.digits[H][c] == 3
                        scalar = (1 - F(1, self.q) if E >> c & 1 else -F(1, self.q)) if exceptional_term else -F(1, self.q) if E >> c & 1 else F(1, self.q - 1) - F(1, self.q)
                        new = {I: value * scalar for I, value in terms.items()}
                        if not exceptional_term and (not (Z | E) >> c & 1):
                            new.update({I | 1 << c: value for I, value in terms.items()})
                        terms = new
                    for I, value in terms.items():
                        if value:
                            values[I][H] = values[I].get(H, F(0)) + value
        else:
            for I, v in values.items():
                for H in self.initial_kept:
                    if self.ordinary[H] & I == I:
                        value = initial_coefficient(self.q, self.spec, self.digits, self.coordinates, mode, Z, E, I, H)
                        if value:
                            v[H] = value
        for step in self.steps:
            updated = {I: {H: value for H, value in vs.items() if H in step['kept']} for I, vs in values.items()}
            if step['free'] == mode:
                for I, b in step['blocks'].items():
                    if I not in values:
                        continue
                    inputs = [(j, values[I][H]) for j, H in enumerate(b['deleted']) if H in values[I]]
                    if not inputs:
                        continue
                    for H, row in zip(b['basis'], b['coefficients']):
                        value = updated[I].get(H, F(0)) + sum((row[j] * v for j, v in inputs if row[j]))
                        if value:
                            updated[I][H] = value
                        else:
                            updated[I].pop(H, None)
            values = updated
        table = {H: {} for H in self.kept}
        for I, vs in values.items():
            for H, value in vs.items():
                table[H][I] = value
        TABLES[digest] = table
        self.tables[key] = table
        return table

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2):
            raise IndexError(mode)
        if len(variable) != self.k or any((not 0 <= x <= self.q for x in variable)):
            raise ValueError('invalid variable')
        H, term = self.term(index)
        Z = sum(((x == 0) << i for i, x in enumerate(variable)))
        E = sum(((variable[i] == 1) << i for i in self.coordinates))
        if Z not in self.spec['families'][mode]:
            return F(0)
        key = (mode, Z, E)
        if key not in self.tables:
            self.build_table(*key)
        value = F(0)
        for I, c in self.tables[key][H].items():
            for i in range(self.k):
                if I >> i & 1:
                    c *= F(int(term[i] == variable[i])) - F(1, self.qs[i])
            value += c
        return value
