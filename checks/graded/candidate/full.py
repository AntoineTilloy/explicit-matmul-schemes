#!/usr/bin/env python3
"""SOURCE ONLY. Exact full balanced capacity via carry-free Kronecker encoding.

Requires reviewed execution on Atlas. All scientific positives use integer
comparisons. Bounds and heuristic q screens do not claim a global minimum.
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
HELPER_SHA = 'f5b30a56c675c6698bf85ec0ff2abcbe7fba1f3e0bff94ea9d6bc9248066df92'
AGGREGATE_SHA = 'c0fc8ad67d6929d3cba17400d22594ed47e291cb451fbd6d0c648dc1e0c7e7d8'


def read_bound(path, sha):
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == sha
    return json.loads(raw)


def save(path, obj):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(obj, indent=2)+'\n')
    temporary.replace(path)


def encoded_capacity(histogram, w, t, integer):
    s = sum(histogram.values())
    digit_bits = (s**(3*t)).bit_length()
    stride = 3*w*t+1
    digit_mask = (integer(1) << digit_bits)-1
    factors = []
    for orientation in range(3):
        encoded = integer(0)
        for (a, b), multiplicity in histogram.items():
            c = w-a-b
            assert min(a, b, c) >= 0
            x, y = ((a, b), (c, a), (b, c))[orientation]
            encoded += integer(multiplicity) << (digit_bits*(x+stride*y))
        factors.append(encoded)
    polynomial = (factors[0]*factors[1]*factors[2])**t
    offset = digit_bits*w*t*(stride+1)
    count = int((polynomial >> offset) & digit_mask)
    assert 0 <= count <= s**(3*t)
    return count, dict(digit_bits=digit_bits, stride=stride,
                       extraction_offset_bits=offset, encoded_bit_bound=digit_bits*stride**2)


def tiny_controls(integer):
    # Includes repeated profiles and profiles with two zero coordinates.
    for profiles, w in (([(2, 0), (0, 1), (0, 1)], 2),
                        ([(1, 1), (2, 0), (0, 2)], 3)):
        histogram = Counter(profiles)
        for t in (1, 2):
            sequences = list(itertools.product(profiles, repeat=t))
            totals = [(sum(p[0] for p in seq), sum(p[1] for p in seq)) for seq in sequences]
            expected = 0
            for (a1, b1), (a2, b2), (a3, b3) in itertools.product(totals, repeat=3):
                c1, c2, c3 = w*t-a1-b1, w*t-a2-b2, w*t-a3-b3
                expected += (a1+c2+b3 == b1+a2+c3 == c1+b2+a3 == w*t)
            got, _ = encoded_capacity(histogram, w, t, integer)
            assert got == expected
            counts = Counter(totals)
            assert got >= sum(value**3 for value in counts.values())
    return 'PASS: all balanced cyclic triples exhaustively checked for four tiny cases'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--helper', type=Path, required=True)
    ap.add_argument('--aggregate', type=Path, required=True)
    ap.add_argument('--grading', type=Path, required=True)
    ap.add_argument('--catalogue', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--max-t', type=int, default=64)
    ap.add_argument('--max-q', type=int, default=64)
    ap.add_argument('--max-case-bits', type=int, default=32000000)
    ap.add_argument('--max-total-bits', type=int, default=256000000)
    args = ap.parse_args()
    assert 1 <= args.max_t <= 128 and 2 <= args.max_q <= 128
    assert 1 <= args.max_case_bits <= 32000000 and 1 <= args.max_total_bits <= 512000000
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    args.out.mkdir(parents=True, exist_ok=True)
    assert hashlib.sha256(args.helper.read_bytes()).hexdigest() == HELPER_SHA
    spec = importlib.util.spec_from_file_location('reviewed_aggregate_helper', args.helper)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    aggregate = read_bound(args.aggregate, AGGREGATE_SHA)
    grading = read_bound(args.grading, helper.GRADING_SHA)
    catalogue = read_bound(args.catalogue, helper.CATALOGUE_SHA)
    outer = next(x for x in catalogue if int(x['n']) == helper.OUTER_SIDE)
    assert int(outer['rank']) == helper.OUTER_RANK and outer['root_object'] == helper.OUTER_ROOT and outer['field'] == 'Q'
    lita = next(x for x in catalogue if int(x['n']) == 32)
    assert int(lita['rank']) == 14196 and lita['field'] == 'Q'
    result = dict(status='RUNNING', source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  helper_sha256=HELPER_SHA, aggregate_sha256=AGGREGATE_SHA,
                  grading_sha256=helper.GRADING_SHA, catalogue_sha256=helper.CATALOGUE_SHA,
                  scope='All balanced tuples in three globally cyclic blocks; bounded pricing, no global minimum claim',
                  limits=dict(max_t=args.max_t, max_q=args.max_q, max_case_bits=args.max_case_bits,
                              max_total_bits=args.max_total_bits),
                  cumulative_encoded_bit_bound=0, first_new_hit=None,
                  best_natural=aggregate['best_natural'], best_padded=aggregate['best_padded'], puzzles=[])

    def checkpoint():
        result['CPU_seconds'] = time.process_time()-START
        save(args.out/'result.json', result)

    checkpoint()
    try:
        import flint
    except ImportError:
        result['status'] = 'UNKNOWN_FLINT_UNAVAILABLE'
        checkpoint()
        return
    result['flint_version'] = flint.__version__
    result['controls'] = tiny_controls(flint.fmpz)
    result['helper_controls'] = helper.controls()
    stop = False
    for wanted in ((23, 7), (35, 8), (14, 6), (8, 5), (5, 4)):
        target = next(x for x in grading['targets'] if tuple(x['puzzle']) == wanted)
        s, w, d, edges, histogram = helper.replay(target)
        prior = next(x for x in aggregate['puzzles'] if tuple(x['puzzle']) == wanted)
        lower_counts = {x['t']: x['square_capacity'] for x in prior['iterations']}
        record = dict(puzzle=[s, w], degree=d, grading_replay='PASS', structural_edges=edges,
                      status='RUNNING', completed_t=0, cases=[])
        result['puzzles'].append(record)
        for t in range(1, args.max_t+1):
            upper = max(result['best_natural']['natural_n'], result['best_padded']['padded_n'])
            if 3**(60*w*t) >= upper**53:
                record['status'] = 'COMPLETE_SIZE_PRUNED'
                break
            bits = (s**(3*t)).bit_length()*(3*w*t+1)**2
            if bits > args.max_case_bits:
                record['status'] = 'UNKNOWN_REMAINDER_CASE_BIT_LIMIT'
                break
            if (time.process_time()-START > 17.5 or
                    result['cumulative_encoded_bit_bound']+bits > args.max_total_bits):
                record['status'] = 'UNKNOWN_REMAINDER_GLOBAL_BUDGET'
                result['status'] = 'PARTIAL_GLOBAL_BUDGET'
                stop = True
                break
            record['pending_t'] = t
            checkpoint()
            capacity, encoding = encoded_capacity(histogram, w, t, flint.fmpz)
            assert encoding['encoded_bit_bound'] == bits
            result['cumulative_encoded_bit_bound'] += bits
            if t in lower_counts:
                assert capacity >= lower_counts[t]
            options = helper.outer_options(capacity)
            largest = options[-1]
            factor = 3*t*d+1
            record['cases'].append(dict(t=t, capacity=capacity, aggregate_lower=lower_counts.get(t),
                                        encoding=encoding, largest_outer=largest))
            for q in range(2, args.max_q+1):
                if time.process_time()-START > 18.0:
                    record['status'] = 'UNKNOWN_REMAINDER_GLOBAL_BUDGET'
                    result['status'] = 'PARTIAL_GLOBAL_BUDGET'
                    stop = True
                    break
                logrank = math.log(factor)+3*w*t*math.log(q+1)
                lognatural = math.log(largest['side'])+w*t*math.log(q)
                if logrank >= 2.65*lognatural+1e-10:
                    continue
                bound = max(math.log(result['best_natural']['natural_n']), math.log(result['best_padded']['padded_n']))
                if logrank/2.65 > bound+1e-10:
                    continue
                rank = factor*(q+1)**(3*w*t)
                inner = q**(w*t)
                rank20 = rank**20
                if rank20 >= (largest['side']*inner)**53:
                    continue
                padded = helper.integer_root(rank20, 53)+1
                chosen = next(o for o in options if o['side']*inner >= padded)
                natural = chosen['side']*inner
                assert chosen['rank'] <= capacity and rank20 < natural**53
                assert (padded-1)**53 <= rank20 < padded**53
                item = dict(puzzle=[s, w], degree=d, t=t, q=q, base_copies=3*t,
                            power=3*w*t, square_capacity=capacity, inner_side=inner, outer=chosen,
                            interpolation_factor=factor, rank=rank, natural_n=natural, padded_n=padded,
                            natural_exponent=math.log(rank)/math.log(natural), exact_threshold='PASS',
                            construction='All balanced row tuples in three globally cyclic blocks; interpolation; pinned public-base/LITA32/Strassen outer family')
                if result['first_new_hit'] is None:
                    result['first_new_hit'] = item
                if (natural, rank) < (result['best_natural']['natural_n'], result['best_natural']['rank']):
                    result['best_natural'] = item
                if (padded, rank) < (result['best_padded']['padded_n'], result['best_padded']['rank']):
                    result['best_padded'] = item
                checkpoint()
            if stop:
                break
            record['completed_t'] = t
            record.pop('pending_t', None)
            checkpoint()
        else:
            record['status'] = 'COMPLETE_BOUNDED_SCREEN'
        if stop:
            break
    else:
        result['status'] = ('PARTIAL_CASE_BIT_LIMIT' if any('UNKNOWN' in p['status'] for p in result['puzzles'])
                            else 'COMPLETE_BOUNDED_SCREEN')
    checkpoint()


if __name__ == '__main__':
    main()
