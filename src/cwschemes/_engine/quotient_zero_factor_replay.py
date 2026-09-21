"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np
from flint import fmpq_mat
from .compressed_fixed_quotient import CompressedFixedQuotient
from .quotient_invariant_gram import coordinate_kernel, representation
from .tailored_completion import layout

class SingleGram:

    def __init__(self, q, rows, m, groups, bases, mode=1):
        spec = layout(tuple(map(tuple, rows)), m)
        self.q = q
        self.groups = groups
        self.mode = mode
        self.character = (-1, -1, 1)[mode]
        self.sectors = sorted(spec['families'][mode])
        self.tables = []
        for group, bs in zip(groups, bases):
            labels = sorted({tuple(map(int, x)) for b in bs for x in b['labels']})
            lookup = {x: i for i, x in enumerate(labels)}
            V = np.zeros((len(bs), len(labels)), dtype=np.int64)
            for i, b in enumerate(bs):
                for x, v in zip(b['labels'], b['vector']):
                    V[i, lookup[tuple(x)]] = v
            labels = np.array(labels, dtype=np.int64)
            cache = {}
            tables = []
            for sector in self.sectors:
                zero = sum(((sector >> pos & 1) << k for k, pos in enumerate(group)))
                for twist in range(2):
                    key = (zero, twist)
                    if key not in cache:
                        K = coordinate_kernel(q, labels, labels, spec['missing'][group[0]], mode, zero, twist, len(group))
                        assert int(np.abs(K).max()) * len(labels) ** 2 < 2 ** 63
                        cache[key] = V @ K @ V.T
                    tables.append(cache[key])
            self.tables.append(tables)

    def entry(self, a, b):
        value = 0
        for sector in range(2 * len(self.sectors)):
            value += self.character ** (sector % 2) * math.prod((int(table[sector][i, j]) for table, i, j in zip(self.tables, a['indices'], b['indices'])))
        return 2 * a['factor'] * b['factor'] * value

def form_norms(adapter, removed):
    q = adapter.q
    metadata = adapter.metadata
    engine = SingleGram(q, metadata['rows'], metadata['m'], metadata['groups'], adapter.bases)
    norms = []
    covered = 0
    for key, block in enumerate(metadata['blocks']):
        rs = block['records']
        K = [i for i, r in enumerate(rs) if r['orbit'] in removed]
        if not K:
            continue
        relation = adapter.proof['relations'].get(str(key))
        J = relation['deleted'] if relation else []
        B = relation['basis'] if relation else []
        if J:
            metric = fmpq_mat([[adapter.metric(rs[i], rs[j]) for j in J] for i in J])
            C = fmpq_mat(relation['coefficients'])
        for k in K:
            if J:
                v = fmpq_mat([[adapter.metric(rs[b], rs[k])] for b in B])
                T = metric.solve(C.transpose() * v)
                weights = [F(1)] + [F(str(T[i, 0])) for i in range(len(J))]
            else:
                weights = [F(1)]
            indices = [k] + J
            nonzero = [(i, w) for i, w in zip(indices, weights) if w]
            indices = [i for i, _ in nonzero]
            weights = [w for _, w in nonzero]
            H = [[engine.entry(rs[i], rs[j]) for j in indices] for i in indices]
            norm = sum((weights[i] * H[i][j] * weights[j] for i in range(len(indices)) for j in range(len(indices))), F(0))
            assert norm == 0, (q, key, k, norm)
            norms.append(dict(block=key, record=k, indices=indices, weights=list(map(str, weights)), integer_gram=H, exact_squared_norm='0'))
            covered += block['dimension']
    return (norms, covered)

def verify_frozen(path):
    from .compressed_fixed_quotient import input_path
    proof = json.loads(input_path(path).read_text())
    assert proof['complete'] and proof['certificate_method'] == 'exact-single-factor-zero-norm'
    assert proof['pair'] == [1, 2] and proof['zero_mode'] == 1 and (not proof['relations'])
    history = proof['history']
    assert len(history) == 2
    previous = [json.loads(input_path(p).read_text()) for p in history]
    assert [p['pair'] for p in previous] == [[0, 1], [0, 2]] and all((p['complete'] for p in previous))
    assert previous[1]['history'] == [history[0]]
    assert all((p['q'] == proof['q'] and p['original_rank'] == proof['original_rank'] for p in previous))
    removed = set(proof['deleted_orbits'])
    assert not removed.intersection(previous[0]['deleted_orbits'] + previous[1]['deleted_orbits'])
    adapter = CompressedFixedQuotient(history[1], allow_stage=True)
    metadata = json.loads(input_path(Path(proof['source']).with_suffix('.json')).read_text())
    for key in ('q', 'm', 'rows', 'groups', 'orbits', 'blocks'):
        assert metadata[key] == adapter.metadata[key]
    norms, covered = form_norms(adapter, removed)
    assert norms == proof['zero_norm_certificates']
    assert covered == proof['covered_irreducible_dimension'] == proof['savings'] == proof['q'] ** 3
    assert proof['rank'] == proof['original_rank'] - covered
    assert proof['cumulative_rank'] == proof['original_rank'] - sum((p['savings'] for p in previous)) - covered
    return dict(q=proof['q'], exact_zero_norms=len(norms), covered_columns=covered)
