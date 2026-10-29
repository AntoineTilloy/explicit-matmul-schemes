#!/usr/bin/env python3
"""Register the exact Strassen-tensor (S_q) packet schemes of October 2026.

Sizes 8,388,608 (replaced), 13,947,137,604 (new) and 13,468,840,704 (replaced, padded).
Reads the hash-bound study certificates (support rows, Ruzsa-Szemeredi sets, parameters)
and replays every support and count check before any repository write. No search is run.

Usage, from the repository root, with the package sources on the path:
    PYTHONPATH=src python checks/sq/export_sq.py --certificates DIR [--write]
"""
import argparse
import gzip
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import sys
import time

REPO = Path(__file__).resolve().parents[2]
INPUTS = {
    'n8388608/S_N18.txt': '505739507ccab290d9955c84e81a0ddb15fc1bed501b9637b48626b07a10f386',
    'n8388608/rs_sets.json': '4a4de73f647a53f1fc632e27c8dd6f6fb0a41d2d5f18d26a47899418165d6c70',
    'n8388608/params.json': '62acc78b218fbf7d11a73c747eecab0ab185eb48678e8e3fe3573da390b899f8',
    'n13468840704/S_N27.txt': '818cbb8b589f797c6d90063f79c2c8eb565dd8f35dfacd64dd6237591535f50d',
    'n13468840704/rs_sets.json': 'c04869275bb1caaa87e64d3f080629ad6510564e7c7f6ebf128d95d1e981022b',
    'n13468840704/params.json': 'b05da0d1b9139a2b7db143f696d74e67e5a78ade7e07d91c5a7a88293bbd9a3a',
}
# Root object of the catalogue size-32 entry; its certificate is reused as the outer scheme.
EXPLICIT32_ROOT = 'a7fb19d69d8e078fb0d6d614ade10417225a186f4fbf94a8d4292f3935995e05'
LITA_COMMIT = '472a79fbc8248423a97bebb3aacbc89e769b9b3e'


def sha(raw): return hashlib.sha256(raw).hexdigest()
def canonical(data): return (json.dumps(data, sort_keys=True, separators=(',', ':')) + '\n').encode()
def pretty(data): return json.dumps(data, indent=2) + '\n'


def read_bound(base, logical):
    raw = (base / logical).read_bytes()
    if sha(raw) != INPUTS[logical]:
        raise ValueError(f'Input hash mismatch: {logical}')
    return raw


def rs_set(R, m):
    """Tripartite Ruzsa-Szemeredi: injective projections and no nontrivial triangle."""
    assert all(0 <= x < m for r in R for x in r) and len(set(R)) == len(R)
    for i, j in ((0, 1), (1, 2), (2, 0)):
        assert len({(r[i], r[j]) for r in R}) == len(R)
    S = set(R)
    ab = {(r[0], r[1]) for r in R}
    ca = {(r[2], r[0]) for r in R}
    by_b = {}
    for r in R:
        by_b.setdefault(r[1], []).append(r[2])
    for a, b in ab:
        for c in by_b[b]:
            assert (c, a) not in ca or (a, b, c) in S, ('triangle', a, b, c)


def rows_from_rs(sets, m):
    """Row (t,a,b,c): A = S0(a) + K2 - S2(c), B = K0 - S0(a) + S1(b), C = K1 - S1(b) + S2(c)."""
    K = [((1 << m) - 1) << (m * j) for j in range(3)]
    out = []
    for t in range(m + 1):
        labels = [sum(1 << i for i in s) for s in combinations(range(m), t)]
        R = [tuple(r) for r in sets[str(len(labels))]]
        rs_set(R, len(labels))
        for a, b, c in R:
            s0, s1, s2 = labels[a], labels[b] << m, labels[c] << 2 * m
            out.append((s0 | (K[2] & ~s2), (K[0] & ~s0) | s1, (K[1] & ~s1) | s2))
    return out


def cross_triples(rows, N):
    """Independent of the package audit: every (u,v) with disjoint A_u, B_v, then C lookup."""
    import numpy as np
    full = (1 << N) - 1
    A = np.array([r[0] for r in rows], dtype=np.int64)
    B = np.array([r[1] for r in rows], dtype=np.int64)
    C = np.array([r[2] for r in rows], dtype=np.int64)
    order = np.argsort(C)
    Cs = C[order]
    diagonal = off = 0
    for u in range(len(rows)):
        ok = (A[u] & B) == 0
        want = full ^ (A[u] | B)
        pos = np.minimum(np.searchsorted(Cs, want), len(rows) - 1)
        for v in np.nonzero(ok & (Cs[pos] == want))[0]:
            w = int(order[pos[v]])
            if u == v == w:
                diagonal += 1
            else:
                off += 1
    return dict(diagonal=diagonal, off_diagonal=off)


def lita36_outer():
    from cwschemes._engine.sources import lita_current
    source = sha(Path(lita_current.__file__).read_bytes())
    assert source == sha((REPO / 'sources/lita_current.py').read_bytes())
    rank = lita_current.lita_rank(36)
    node = dict(children=[], commit=LITA_COMMIT, operation='lita_current', rank=rank,
                source_sha256=source, tensor=[36, 36, 36])
    key = sha(json.dumps(node, sort_keys=True, separators=(',', ':')).encode())
    certificate = dict(format='exact-rectangular-recipe-v1', coefficient_field='Q', n=36, rank=rank,
                       tensor=[36, 36, 36], root=key, nodes={key: node},
                       construction=f'Current LITA row-column aggregation scheme at size 36 (github.com/khoruzhii/lita, revision {LITA_COMMIT[:7]}, rank {rank}), indexed from the pinned generator sources/lita_current.py.',
                       primitive_identity_basis='Pinned current LITA construction theorem with hash-bound indexed coefficients; the size-36 tensor was also checked exactly over Q by an independent CRT replay. See checks/sq/independent.json.')
    return dict(kind='cw', n=36, certificate=certificate)


def explicit32_outer():
    raw = gzip.decompress((REPO / f'src/cwschemes/data/objects/{EXPLICIT32_ROOT}.gz').read_bytes())
    assert sha(raw) == EXPLICIT32_ROOT
    certificate = json.loads(raw)
    assert (certificate['format'], certificate['n'], certificate['rank']) == ('explicit-qcsr-square-v1', 32, 14196)
    return dict(kind='cw', n=32, certificate=certificate)


def comparison_13947137604():
    """A specific published construction: current LITA 36 (pinned) times three size-729 closure recipes."""
    sys.path.insert(0, str(REPO / 'comparisons'))
    import build
    seeds = json.loads((REPO / 'comparisons/primitives.json').read_text())
    ranks, recipes = build.close(729, seeds)
    leaf = recipes[729]
    inner = dict(kind='product', n=729 ** 3, left=leaf, right=dict(kind='product', n=729 ** 2, left=leaf, right=leaf))
    recipe = dict(kind='product', n=36 * 729 ** 3, left=dict(kind='lita-even', n=36), right=inner)
    assert recipe['n'] == 36 * 9 ** 9 and build.evaluate(recipe, seeds) == build.lita(36) * ranks[729] ** 3
    return recipe, build.evaluate(recipe, seeds)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificates', type=Path, required=True)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    start = time.monotonic()
    sys.path.insert(0, str(REPO / 'src'))
    from cwschemes import _resources
    for logical in _resources.index()['paths']:
        if logical.startswith('seeds/'):
            _resources.materialize(logical)
    from cwschemes._engine import strassen_packet
    from cwschemes._engine.general_square import MixedSquare, load_certificate

    supports = {}
    for name in ('n8388608', 'n13468840704'):
        params = json.loads(read_bound(args.certificates, f'{name}/params.json'))
        N, q = params['N'], params['q']
        m = N // 3
        lines = read_bound(args.certificates, f'{name}/S_N{N}.txt').decode().split('\n')
        header = tuple(map(int, lines[0].split()))
        rows = [tuple(map(int, line.split())) for line in lines[1:] if line.strip()]
        assert header == (N, len(rows))
        sets = json.loads(read_bound(args.certificates, f'{name}/rs_sets.json'))
        rebuilt = rows_from_rs(sets, m)
        assert sorted(rebuilt) == sorted(rows), 'stored rows differ from the RS witness'
        witness = strassen_packet.audit(rows, m)
        independent = cross_triples(rows, N)
        assert independent == dict(diagonal=len(rows), off_diagonal=0)
        supports[N] = dict(q=q, m=m, rows=rows, sets=sets, record=dict(
            power=N, products=len(rows), stored_row_order='certificate order (defines the call assignment)',
            rows_sha256=sha(canonical([list(r) for r in rows])),
            rs_sizes={k: len(v) for k, v in sets.items()}, rows_rebuilt_from_rs_sets=True,
            package_audit=witness, independent_cross_triple_enumeration=independent,
            local_identity_coefficients=strassen_packet.verify_local(q)))

    def spec(N):
        s = supports[N]
        return dict(local_tensor=strassen_packet.LOCAL_TENSOR, q=s['q'], m=s['m'],
                    rows=[list(r) for r in s['rows']], rank=(s['q'] + 1) ** N)

    builds = [
        (8388608, explicit32_outer(), 18, 40),
        (13947137604, lita36_outer(), 27, 2),
        (13468840704, lita36_outer(), 27, 2),
    ]
    catalogue = json.loads((REPO / 'catalogue.json').read_text())
    old = {r['n']: r for r in catalogue}
    objects, changes = {}, []
    for n, outer, N, copies in builds:
        q, m = supports[N]['q'], supports[N]['m']
        leaf = dict(kind='classical', n=q ** m)
        square = MixedSquare(outer, leaf, [spec(N)] * copies, n, leftover_products=0)
        rank = copies * (q + 1) ** N
        assert square.rank == rank and square.leaves == 0 and n <= square.full_n == outer['n'] * q ** m
        assert sum(square.used) == square.outer.rank and (copies - 1) * len(supports[N]['rows']) < square.outer.rank
        data = square.certificate()
        restored = load_certificate(json.loads(canonical(data)))
        assert (restored.n, restored.rank, restored.used) == (n, rank, square.used)
        raw = canonical(data)
        digest = sha(raw)
        objects[digest] = raw
        previous = old.get(str(n))
        if previous is not None:
            assert rank < int(previous['rank'])
        changes.append(dict(n=str(n), old_rank=previous['rank'] if previous else None, rank=str(rank),
                            full_n=str(square.full_n), outer=f"size {outer['n']}, {square.outer.rank} products",
                            packets=copies, packet_power=N, q=q, rows_per_packet=len(supports[N]['rows']),
                            packet_products_used=square.used, unused_rows=copies * len(supports[N]['rows']) - square.outer.rank,
                            root_object=digest, sub_2_7=rank ** 10 < n ** 27,
                            effective_exponent=math.log(rank) / math.log(n)))
    comparison, comparison_rank = comparison_13947137604()

    def record(change):
        n = int(change['n'])
        if n == 8388608:
            construction = 'Strassen-tensor S_q packets (October 2026): explicit size-32 outer scheme (14,196 products) with 40 power-18 S_8 packets of 358 independent products'
        elif n == 13947137604:
            construction = 'Strassen-tensor S_q packets (October 2026): current LITA size-36 outer scheme (19,666 products) with 2 power-27 S_9 packets of 9,842 independent products'
        else:
            construction = 'Padding of the size-13,947,137,604 Strassen-tensor S_q packet scheme (October 2026)'
        return dict(format='cwschemes-catalogue-entry-v1', n=change['n'], rank=change['rank'], field='Q',
                    root_object=change['root_object'], effective_exponent=change['effective_exponent'],
                    provenance=dict(construction=construction,
                                    packet_support='Tripartite Ruzsa-Szemeredi sets per weight class; witnesses in checks/sq/rs_sets.json'),
                    verification='Exact finite S_q packet supports and rational coefficient access included; cwschemes verify replays the support audit, exact local identity, packet counts and sampled exact packet entries. Independent checks recorded in checks/sq/. See docs/verification.md for check scope.')

    new_records = {c['n']: record(c) for c in changes}
    updated = [new_records.get(r['n'], r) for r in catalogue]
    position = next(i for i, r in enumerate(updated) if r['n'] == '13468840704') + 1
    updated.insert(position, new_records['13947137604'])
    assert [int(r['n']) for r in updated] == sorted(int(r['n']) for r in updated)
    large_path = REPO / 'comparisons/large_recipes.json'
    large = json.loads(large_path.read_text())
    large_items = []
    for key, value in large.items():
        large_items.append((key, value))
        if key == '13468840704':
            large_items.append(('13947137604', comparison))
    large = dict(large_items)

    rs_witness = {str(N): dict(power=N, m=s['m'], products=len(s['rows']),
                               convention='Labels of weight class t are the t-subsets of range(m) in lexicographic order; row (t,a,b,c) has A = S0(a) + K2 - S2(c), B = K0 - S0(a) + S1(b), C = K1 - S1(b) + S2(c), with K_j the coordinates of kind j and S_j the label subset placed on K_j.',
                               rs_sets=s['sets'])
                  for N, s in sorted(supports.items())}
    export = dict(status='PASS', script_sha256=sha(Path(__file__).read_bytes()),
                  input_hashes=INPUTS, previous_catalogue_sha256=sha((REPO / 'catalogue.json').read_bytes()),
                  previous_index_sha256=sha((REPO / 'src/cwschemes/data/index.json').read_bytes()),
                  supports={str(N): s['record'] for N, s in sorted(supports.items())},
                  changes=changes,
                  comparison_13947137604=dict(rank=str(comparison_rank),
                                              exponent=math.log(comparison_rank) / math.log(13947137604),
                                              scope='specific published construction: pinned LITA 36 times three closure recipes at 729'),
                  catalogue_entries_before=len(catalogue), catalogue_entries_after=len(updated),
                  seconds=time.monotonic() - start)
    if args.write:
        index_path = REPO / 'src/cwschemes/data/index.json'
        index = json.loads(index_path.read_text())
        for digest, raw in objects.items():
            path = REPO / f'src/cwschemes/data/objects/{digest}.gz'
            packed = gzip.compress(raw, compresslevel=9, mtime=0)
            path.write_bytes(packed)
            index['objects'][digest] = dict(bytes=len(raw), compressed_bytes=len(packed))
            index['paths'][f'proofs/{digest}.json'] = digest
        index_path.write_text(json.dumps(index, indent=2) + '\n')
        for path in (REPO / 'catalogue.json', REPO / 'src/cwschemes/data/catalogue.json'):
            path.write_text(pretty(updated))
        for r in new_records.values():
            (REPO / f"schemes/{r['n']}.json").write_text(pretty(r))
        large_path.write_text(pretty(large))
        (REPO / 'checks/sq/rs_sets.json').write_text(pretty(rs_witness))
        (REPO / 'checks/sq/export.json').write_text(pretty(export))
    print(json.dumps(dict(status=export['status'], changes=[(c['n'], c['rank'], c['root_object']) for c in changes],
                          seconds=export['seconds'], written=args.write)))


if __name__ == '__main__':
    main()
