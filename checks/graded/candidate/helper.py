#!/usr/bin/env python3
"""Bounded constructive pricing using equal aggregate profiles in cyclic blocks.

SOURCE ONLY until reviewed. Scientific execution must be on Atlas.
No numerical tolerances certify positives or global minimum claims.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import resource
import sys
import time
from collections import Counter
from pathlib import Path

START = time.process_time()
sys.set_int_max_str_digits(0)
LIBRARY_SHA = '89a3f3ceb8af898c122b354ecbe3f89464d28fee17d1e386049c51792378238d'
GRADING_SHA = 'c2199f0c9b59129ac18eeb71b8da891d433a1a5852a9b3db2f02a7f8dd72fd0c'
CATALOGUE_SHA = '36922b9293dba32191be65b09f2ba5dd23038f11c778aa33ca9c2c8d0a4df628'
OUTER_SIDE = 12310160190
OUTER_RANK = 1463616394485527559895695360
OUTER_ROOT = 'f21e6922b24860603e77cccc387c07b1764b96c3b12b44d267675169b52acfaf'


def save(path, obj):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(obj, indent=2) + '\n')
    temp.replace(path)


def integer_root(value, degree):
    if value < 2:
        return value
    x = 1 << ((value.bit_length() + degree - 1) // degree)
    while True:
        y = ((degree - 1) * x + value // x**(degree - 1)) // degree
        if y >= x:
            assert x**degree <= value < (x+1)**degree
            return x
        x = y


def controls():
    # Repeated profiles exercise the multiplicity factors in the recurrence.
    profiles = [(2, 0), (0, 1), (0, 1)]
    histogram = Counter(profiles)
    current = {(0, 0): 1}
    for t in range(1, 4):
        nxt = Counter()
        for (a, b), count in current.items():
            for (x, y), multiplicity in histogram.items():
                nxt[a+x, b+y] += count*multiplicity
        brute = Counter((sum(p[0] for p in seq), sum(p[1] for p in seq))
                        for seq in itertools.product(profiles, repeat=t))
        assert nxt == brute
        capacity = sum(count**3 for count in nxt.values())
        tuples = list(itertools.product(profiles, repeat=t))
        totals = [(sum(p[0] for p in seq), sum(p[1] for p in seq)) for seq in tuples]
        assert capacity == sum(a == b == c for a, b, c in itertools.product(totals, repeat=3))
        current = nxt
    for degree in (1, 2, 3, 53):
        for value in range(100):
            got = integer_root(value, degree)
            assert got**degree <= value < (got+1)**degree
    return 'PASS: weighted profile DP, triple-capacity enumeration, exact integer roots'


def replay(target):
    assert target['status'] == 'CERTIFIED_POSITIVE_GRADING'
    rows = target['rows']
    s, w = target['puzzle']
    assert len(rows) == s and all(len(row) == w and set(row) <= set('123') for row in rows)
    potential = target['certificate']['potentials']
    assert len(potential) == 3 and all(len(p) == s and all(type(x) is int for x in p) for p in potential)
    masks = [[sum(1 << j for j, ch in enumerate(row) if ch == str(h+1))
              for row in rows] for h in range(3)]
    assert all(len(set(mode)) == s for mode in masks)
    degree, edges = 0, 0
    for i, j, k in itertools.product(range(s), repeat=3):
        A, B, C = masks[0][i], masks[1][j], masks[2][k]
        if A&B == A&C == B&C:
            edges += 1
            grade = potential[0][i] + potential[1][j] + potential[2][k]
            if i == j == k:
                assert grade == 0
            else:
                assert grade > 0
                degree = max(degree, grade)
    assert degree == target['certificate']['degree']
    return s, w, degree, edges, Counter((row.count('1'), row.count('2')) for row in rows)


def outer_options(capacity):
    # All tensor products of the pinned best public base, LITA32 and at most
    # four Strassen levels. Five Strassen levels are dominated by LITA32.
    assert 14196 < 7**5
    options = []
    h, rank_h, side_h = 0, 1, 1
    while rank_h <= capacity:
        i, rank_i, side_i = 0, rank_h, side_h
        while rank_i <= capacity:
            rank, side = rank_i, side_i
            for j in range(5):
                if rank <= capacity:
                    options.append(dict(side=side, rank=rank, public_base_copies=h,
                                        public_root=OUTER_ROOT, lita32_copies=i,
                                        strassen2_copies=j))
                rank *= 7
                side *= 2
            i += 1
            rank_i *= 14196
            side_i *= 32
        h += 1
        rank_h *= OUTER_RANK
        side_h *= OUTER_SIDE
    return sorted(options, key=lambda x: (x['side'], x['rank']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grading', type=Path, required=True)
    ap.add_argument('--library', type=Path, required=True)
    ap.add_argument('--catalogue', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--max-t', type=int, default=64)
    ap.add_argument('--max-q', type=int, default=64)
    ap.add_argument('--max-transitions', type=int, default=4000000)
    args = ap.parse_args()
    assert 1 <= args.max_t <= 128 and 2 <= args.max_q <= 128
    assert 1 <= args.max_transitions <= 8000000
    resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    raw = args.grading.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == GRADING_SHA
    assert hashlib.sha256(args.library.read_bytes()).hexdigest() == LIBRARY_SHA
    catalogue_raw = args.catalogue.read_bytes()
    assert hashlib.sha256(catalogue_raw).hexdigest() == CATALOGUE_SHA
    catalogue = json.loads(catalogue_raw)
    outer = next(x for x in catalogue if int(x['n']) == OUTER_SIDE)
    assert int(outer['rank']) == OUTER_RANK and outer['field'] == 'Q' and outer['root_object'] == OUTER_ROOT
    lita = next(x for x in catalogue if int(x['n']) == 32)
    assert int(lita['rank']) == 14196 and lita['field'] == 'Q'
    spec = importlib.util.spec_from_file_location('reviewed_pricer', args.library)
    library = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(library)
    args.out.mkdir(parents=True, exist_ok=True)
    data = json.loads(raw)
    result = dict(status='RUNNING', source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  grading_sha256=GRADING_SHA, library_sha256=LIBRARY_SHA,
                  catalogue_sha256=CATALOGUE_SHA, public_outer=outer,
                  controls=controls(), local_G_controls=library.local_controls(),
                  scope='Equal aggregate profiles in three globally cyclic blocks; exact capacities and positive checks; bounded heuristic pricing, no global minimum claim',
                  ranges=dict(max_t=args.max_t, max_q=args.max_q, max_transitions=args.max_transitions),
                  transitions=0, first_hit=None, best_natural=None, best_padded=None, puzzles=[])

    def checkpoint():
        result['CPU_seconds'] = time.process_time() - START
        save(args.out/'result.json', result)

    stop = False
    for wanted in ((35, 8), (23, 7), (14, 6), (8, 5), (5, 4)):
        target = next(x for x in data['targets'] if tuple(x['puzzle']) == wanted)
        s, w, d, edges, histogram = replay(target)
        record = dict(puzzle=[s, w], degree=d, structural_edges=edges, replay='PASS',
                      profile_histogram=[[a, b, c] for (a, b), c in sorted(histogram.items())],
                      completed_t=0, status='RUNNING', iterations=[])
        result['puzzles'].append(record)
        current = {(0, 0): 1}
        for t in range(1, args.max_t+1):
            if result['best_natural'] and result['best_padded']:
                upper = max(result['best_natural']['natural_n'], result['best_padded']['padded_n'])
                # Every future packet costs at least 3^(3wt), since q>=2.
                # Thus even an arbitrarily good outer cannot improve either
                # saved size once this exact strict-exponent bound is met.
                if 3**(60*w*t) >= upper**53:
                    record['status'] = 'COMPLETE_SIZE_PRUNED'
                    break
            needed = len(current)*len(histogram)
            if (time.process_time()-START > 8.0 or
                    result['transitions']+needed > args.max_transitions or len(current) > 100000):
                record['status'] = 'UNKNOWN_REMAINDER_BUDGET'
                result['status'] = 'PARTIAL_BUDGET'
                stop = True
                break
            checkpoint()
            nxt = Counter()
            for (a, b), count in current.items():
                for (x, y), multiplicity in histogram.items():
                    nxt[a+x, b+y] += count*multiplicity
            result['transitions'] += needed
            current = nxt
            assert sum(current.values()) == s**t
            capacity = sum(count**3 for count in current.values())
            assert capacity >= max(current.values())**3
            options = outer_options(capacity)
            largest = options[-1]
            factor = 3*t*d+1
            iteration = dict(t=t, states=len(current), square_capacity=capacity,
                             largest_outer=largest, interpolation_factor=factor)
            record['iterations'].append(iteration)
            for q in range(2, args.max_q+1):
                if time.process_time()-START > 8.5:
                    record['status'] = 'UNKNOWN_REMAINDER_BUDGET'
                    result['status'] = 'PARTIAL_BUDGET'
                    stop = True
                    break
                logrank = math.log(factor) + 3*w*t*math.log(q+1)
                lognatural = math.log(largest['side']) + w*t*math.log(q)
                if logrank >= 2.65*lognatural+1e-10:
                    continue
                if result['best_natural'] and result['best_padded']:
                    bound = max(math.log(result['best_natural']['natural_n']),
                                math.log(result['best_padded']['padded_n']))
                    if logrank/2.65 > bound+1e-10:
                        continue
                rank = factor*(q+1)**(3*w*t)
                inner = q**(w*t)
                rank20 = rank**20
                if rank20 >= (largest['side']*inner)**53:
                    continue
                padded = integer_root(rank20, 53)+1
                chosen = next(o for o in options if o['side']*inner >= padded)
                natural = chosen['side']*inner
                assert chosen['rank'] <= capacity and rank20 < natural**53
                assert (padded-1)**53 <= rank20 < padded**53
                item = dict(puzzle=[s, w], degree=d, t=t, q=q, base_copies=3*t,
                            power=3*w*t, square_capacity=capacity, inner_side=inner,
                            outer=chosen, interpolation_factor=factor, rank=rank,
                            natural_n=natural, padded_n=padded,
                            natural_exponent=math.log(rank)/math.log(natural), exact_threshold='PASS',
                            construction='Direct sum over equal aggregate-profile triples; globally cyclic blocks; interpolation; pinned public-base/LITA32/Strassen outer family')
                if result['first_hit'] is None:
                    result['first_hit'] = item
                if result['best_natural'] is None or (natural, rank) < (result['best_natural']['natural_n'], result['best_natural']['rank']):
                    result['best_natural'] = item
                if result['best_padded'] is None or (padded, rank) < (result['best_padded']['padded_n'], result['best_padded']['rank']):
                    result['best_padded'] = item
                checkpoint()
            if stop:
                break
            record['completed_t'] = t
            checkpoint()
        else:
            record['status'] = 'COMPLETE_BOUNDED_SCREEN'
        if stop:
            break
    else:
        result['status'] = 'COMPLETE_BOUNDED_SCREEN'
    checkpoint()


if __name__ == '__main__':
    main()
