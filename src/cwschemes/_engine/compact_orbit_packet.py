"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib, json
from .scheme import CompactPacket
from .triangular_full_compression import allowed
from .rational_state import initial_free
ROOT = resource_root()

@lru_cache(16)
def proof_data(path, mtime_ns, size):
    raw = Path(path).read_bytes()
    d = json.loads(raw)
    steps = [dict(free=3 - sum(step['pair']), kept=frozenset(step['kept']), blocks={int(I): dict(basis=b['basis'], deleted=b['deleted'], coefficients=[[F(x) for x in row] for row in b['coefficients']]) for I, b in step['blocks'].items()}) for step in d['steps']]
    return (raw, d, steps)
TABLES = {}

class CompactOrbitPacket:
    compression = 'joint-orbit-v1'

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
        assert d['format'] == 'cw-compact-orbit-relations-v1' and d['complete'] and (d['q'] == q) and (d['m'] == m)
        assert sorted(map(tuple, rows)) == sorted(map(tuple, d['rows']))
        self.base = CompactPacket(q, rows, m, patched=True)
        self.families = [{r[h] for r in rows} for h in range(3)]
        self.initial_kept = allowed(rows, self.k, (0, 1))
        self.kept = d['steps'][-1]['kept']
        self.digits = list(product(range(3), repeat=self.k))
        self.ordinary = [sum(((v == 1) << i for i, v in enumerate(ds))) for ds in self.digits]
        self.ends = []
        total = 0
        for H in self.kept:
            total += q ** self.ordinary[H].bit_count()
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
        digits = [self.q + 1 if t == 0 else self.q for t in self.digits[H]]
        for i in reversed(range(self.k)):
            if self.ordinary[H] >> i & 1:
                local, digits[i] = divmod(local, self.q)
        return (H, digits)

    def build_table(self, mode, zero):
        key = (mode, zero)
        digest = hashlib.sha256(self.proof_bytes + json.dumps(key).encode()).hexdigest()
        if digest in TABLES:
            self.tables[key] = TABLES[digest]
            return TABLES[digest]
        path = ROOT / 'coefficient_cache/joint' / f'{digest}.json'
        if path.exists():
            d = json.loads(path.read_text())
            assert d['digest'] == digest
            table = {int(H): {int(I): F(a, b) for I, a, b in values} for H, values in d['values']}
            TABLES[digest] = table
            self.tables[key] = table
            return table
        q = self.q
        k = self.k
        values = {I: {} for I in range(1 << k) if not I & zero}
        if mode == 2:
            raw, codes, standards = initial_free(q, k, zero, self.initial_kept)
            for value, H, I in zip(raw, codes, standards):
                if int(H) in self.initial_kept:
                    values[int(I)][int(H)] = F(str(value))
        else:
            for I in values:
                for H in self.initial_kept:
                    if self.ordinary[H] & I != I:
                        continue
                    value = F(1)
                    for i, t in enumerate(self.digits[H]):
                        if I >> i & 1 or zero >> i & 1:
                            continue
                        value *= F(0) if t == 0 else F(2, q) if t == 1 else F(3, 2 * q)
                    values[I][H] = value
        for step in self.steps:
            updated = {I: {H: value for H, value in v.items() if H in step['kept']} for I, v in values.items()}
            if step['free'] == mode:
                for I, b in step['blocks'].items():
                    if I not in values:
                        continue
                    inputs = [values[I].get(H, F(0)) for H in b['deleted']]
                    for H, row in zip(b['basis'], b['coefficients']):
                        updated[I][H] = updated[I].get(H, F(0)) + sum((c * v for c, v in zip(row, inputs) if c and v))
            values = updated
        table = {H: {} for H in self.kept}
        for I, vs in values.items():
            for H, value in vs.items():
                if value:
                    table[H][I] = value
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix('.tmp')
        tmp.write_text(json.dumps(dict(digest=digest, values=[[H, [[I, v.numerator, v.denominator] for I, v in vs.items()]] for H, vs in table.items()]), separators=(',', ':')) + '\n')
        tmp.replace(path)
        TABLES[digest] = table
        self.tables[key] = table
        return table

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2):
            raise IndexError(mode)
        if len(variable) != self.k or any((not 0 <= x <= self.q for x in variable)):
            raise ValueError('invalid variable')
        H, term = self.term(index)
        zero = sum(((x == 0) << i for i, x in enumerate(variable)))
        if zero not in self.families[mode]:
            return F(0)
        key = (mode, zero)
        if key not in self.tables:
            self.build_table(mode, zero)
        result = F(0)
        for I, value in self.tables[key][H].items():
            for i in range(self.k):
                if I >> i & 1:
                    value *= F(int(term[i] == variable[i] - 1)) - F(1, self.q)
            result += value
        return result
