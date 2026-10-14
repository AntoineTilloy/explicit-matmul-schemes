#!/usr/bin/env python3
"""SOURCE ONLY until review. Independent finite full-balanced capacity replay.

Ordinary 2-D profile DP followed by pair enumeration, no Kronecker encoding
or research-module imports. Execute on Atlas only after authorization.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import resource
import sys
import time

GRADING_SHA = 'c2199f0c9b59129ac18eeb71b8da891d433a1a5852a9b3db2f02a7f8dd72fd0c'
CATALOGUE_SHA = '36922b9293dba32191be65b09f2ba5dd23038f11c778aa33ca9c2c8d0a4df628'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def profile_counts(histogram, t, guard):
    current = {(0, 0): 1}
    for _ in range(t):
        guard()
        nxt = Counter()
        for (a, b), count in current.items():
            for (x, y), multiplicity in histogram.items():
                nxt[a+x, b+y] += count*multiplicity
        current = dict(nxt)
    assert sum(current.values()) == sum(histogram.values())**t
    return current


def pair_capacity(counts, total, guard):
    """Orientations (a,b,c), (b,c,a), (c,a,b), as coefficient reader.

    a1+b2+c3=W and b1+c2+a3=W imply
    (a3,b3)=(a2+b2-b1, a1+b1-a2). The third equation follows.
    """
    answer = 0
    items = list(counts.items())
    for (a, b), first in items:
        guard()
        for (d, e), second in items:
            third_a, third_b = d+e-b, a+b-d
            if min(third_a, third_b) < 0 or third_a+third_b > total:
                continue
            third = counts.get((third_a, third_b), 0)
            if third:
                answer += first*second*third
    assert answer >= sum(c**3 for c in counts.values())
    assert answer <= sum(counts.values())**3
    return answer


def controls(guard):
    records = []
    for rows, w in (([(2, 0), (0, 1), (0, 1)], 2),
                    ([(1, 1), (2, 0), (0, 2)], 3)):
        for t in (1, 2):
            sequences = list(product(rows, repeat=t))
            totals = [(sum(p[0] for p in seq), sum(p[1] for p in seq)) for seq in sequences]
            counts = profile_counts(Counter(rows), t, guard)
            assert counts == Counter(totals)
            old_orientation = reverse_orientation = 0
            W = w*t
            for (a, b), (d, e), (g, h) in product(totals, repeat=3):
                c, f, i = W-a-b, W-d-e, W-g-h
                old_orientation += (a+e+i == b+f+g == c+d+h == W)
                reverse_orientation += (a+f+h == b+d+i == c+e+g == W)
            got = pair_capacity(counts, W, guard)
            assert got == old_orientation == reverse_orientation
            records.append(dict(width=w, t=t, direct_sequence_triples=len(totals)**3,
                                balanced=got, both_cyclic_orientations='PASS'))
    # Independent local rational coefficient expansion, not a helper import.
    for q in (2, 3):
        forms = [[F(1)] + [F(i == j)+F(1, q) for j in range(q)] for i in range(q)]
        forms += [[F(1)] + [F(3, 2*q)]*q]
        for x, y, z in product(range(q+1), repeat=3):
            actual = sum(f[x]*f[y]*f[z] for f in forms[:q]) - F(4*q, 3)*forms[q][x]*forms[q][y]*forms[q][z]
            zeros = sum(v == 0 for v in (x, y, z))
            if zeros == 3:
                expected = -F(q, 3)
            elif zeros == 2:
                expected = F(0)
            elif zeros == 1:
                nonzero = [v for v in (x, y, z) if v]
                expected = F(nonzero[0] == nonzero[1])
            else:
                expected = F(x == y == z) + F((x == y)+(x == z)+(y == z), q) - F(1, 2*q*q)
            assert actual == expected
    return dict(exhaustive_cyclic_cases=records, local_Q_identity='PASS')


def replay_grading(target, guard):
    s, w = target['puzzle']
    rows = target['rows']
    assert len(rows) == len(set(rows)) == s
    assert all(len(row) == w and set(row) <= set('123') for row in rows)
    masks = [[frozenset(j for j, ch in enumerate(row) if ch == str(mode+1))
              for row in rows] for mode in range(3)]
    assert all(len(set(mode)) == s for mode in masks)
    potential = target['certificate']['potentials']
    assert len(potential) == 3 and all(len(mode) == s and all(type(v) is int for v in mode) for mode in potential)
    degree = unwanted = wanted = 0
    for i in range(s):
        guard()
        for j, k in product(range(s), repeat=2):
            owners = (i, j, k)
            if any(sum(pos in masks[mode][owners[mode]] for mode in range(3)) == 2 for pos in range(w)):
                continue
            grade = sum(potential[mode][owners[mode]] for mode in range(3))
            if i == j == k:
                assert grade == 0
                wanted += 1
            else:
                assert grade > 0
                unwanted += 1
                degree = max(degree, grade)
    assert wanted == s and degree == target['certificate']['degree']
    return dict(degree=degree, wanted=wanted, unwanted=unwanted), Counter(
        (row.count('1'), row.count('2')) for row in rows)


def main():
    p = argparse.ArgumentParser()
    for name in ('grading', 'result', 'catalogue', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--expected-result-sha', required=True)
    p.add_argument('--max-pairs', type=int, default=8000000)
    a = p.parse_args()
    assert 1 <= a.max_pairs <= 8000000
    resource.setrlimit(resource.RLIMIT_CPU, (9, 10))
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    sys.set_int_max_str_digits(0)
    started = time.process_time()
    def guard():
        if time.process_time()-started > 8:
            raise TimeoutError('Internal 8 CPU-second guard')
    assert sha(a.grading) == GRADING_SHA and sha(a.catalogue) == CATALOGUE_SHA
    assert sha(a.result) == a.expected_result_sha
    grading, pricing, catalogue = [json.loads(path.read_text()) for path in (a.grading, a.result, a.catalogue)]
    assert pricing['grading_sha256'] == GRADING_SHA and pricing['catalogue_sha256'] == CATALOGUE_SHA
    lita = [x for x in catalogue if int(x['n']) == 32 and int(x['rank']) == 14196 and x['field'] == 'Q']
    assert lita
    report = dict(status='RUNNING', source_sha256=sha(Path(__file__)),
                  input_hashes={name: sha(getattr(a, name)) for name in ('grading', 'result', 'catalogue')},
                  scope='Independent capacity, grading and finite-price replay; no coefficient-generator audit',
                  limits=dict(CPU_seconds=10, MiB=512, pairs=a.max_pairs),
                  pairs_used=0, cases=[], items=[], inherited_aggregate_items=[])
    def save():
        report['CPU_seconds'] = time.process_time()-started
        temporary = a.out.with_suffix('.tmp')
        temporary.write_text(json.dumps(report, indent=2)+'\n')
        temporary.replace(a.out)
    save()
    try:
        report['controls'] = controls(guard)
        save()
        replayed, capacity_cache = {}, {}
        for label in ('first_new_hit', 'best_natural', 'best_padded'):
            item = pricing.get(label)
            if item is None:
                continue
            if not item['construction'].startswith(('All balanced row tuples', 'Full balanced cyclic packet')):
                report['inherited_aggregate_items'].append(label)
                continue
            s, w = item['puzzle']
            t, q = item['t'], item['q']
            assert type(t) is int and t >= 1 and type(q) is int and q >= 2
            if (s, w) not in replayed:
                target = next(x for x in grading['targets'] if x['puzzle'] == [s, w])
                replayed[s, w] = replay_grading(target, guard)
            replay, histogram = replayed[s, w]
            key = (s, w, t)
            if key not in capacity_cache:
                counts = profile_counts(histogram, t, guard)
                assert all(min(a0, b0) >= 0 and a0+b0 <= w*t for a0, b0 in counts)
                needed = len(counts)**2
                if report['pairs_used']+needed > a.max_pairs:
                    report['status'] = 'UNKNOWN_PAIR_LIMIT'
                    report['pending'] = dict(puzzle=[s, w], t=t, states=len(counts), pairs=needed)
                    save()
                    return
                capacity = pair_capacity(counts, w*t, guard)
                report['pairs_used'] += needed
                capacity_cache[key] = capacity
                table = [[a0, b0, c] for (a0, b0), c in sorted(counts.items())]
                report['cases'].append(dict(puzzle=[s, w], t=t, states=len(counts), pairs=needed,
                    capacity=capacity, grading=replay,
                    profile_table_sha256=hashlib.sha256(json.dumps(table, separators=(',', ':')).encode()).hexdigest()))
                save()
            capacity = capacity_cache[key]
            degree = replay['degree']
            assert item['square_capacity'] == capacity and item['degree'] == degree
            assert item['base_copies'] == 3*t and item['power'] == 3*w*t
            factor = 3*t*degree+1
            batches = item.get('batches', 1)
            assert type(batches) is int and batches >= 1
            packet_rank = factor*(q+1)**(3*w*t)
            if 'packet_rank' in item:
                assert item['packet_rank'] == packet_rank
            rank = batches*packet_rank
            assert item['interpolation_factor'] == factor and item['rank'] == rank
            assert item['inner_side'] == q**(w*t)
            outer = item['outer']
            public = next(x for x in catalogue if x['root_object'] == outer['public_root'] and x['field'] == 'Q')
            if 'public_base_n' in outer:
                assert int(outer['public_base_n']) == int(public['n'])
            if 'public_base_rank' in outer:
                assert int(outer['public_base_rank']) == int(public['rank'])
            h, i, j = outer['public_base_copies'], outer['lita32_copies'], outer['strassen2_copies']
            assert all(type(v) is int and v >= 0 for v in (h, i, j))
            side = int(public['n'])**h * 32**i * 2**j
            calls = int(public['rank'])**h * 14196**i * 7**j
            assert outer['side'] == side and outer['rank'] == calls
            assert batches == (calls+capacity-1)//capacity
            assert (batches-1)*capacity < calls <= batches*capacity
            natural, cropped = side*q**(w*t), item['padded_n']
            assert natural == item['natural_n'] and cropped <= natural
            assert rank**20 < natural**53 and (cropped-1)**53 <= rank**20 < cropped**53
            report['items'].append(dict(item=label, puzzle=[s, w], t=t, q=q, public_root=public['root_object'],
                                        batches=batches, packet_rank=packet_rank,
                                        exact_capacity='PASS', exact_batch_count='PASS',
                                        exact_price='PASS', exact_crop='PASS'))
            save()
        report['status'] = 'PASS' if report['items'] else 'NO_NEW_BALANCED_CANDIDATE'
    except TimeoutError as error:
        report['status'] = 'UNKNOWN_CPU_LIMIT'
        report['error'] = str(error)
    save()


if __name__ == '__main__':
    main()
