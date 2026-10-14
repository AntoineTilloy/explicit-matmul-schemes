#!/usr/bin/env python3
"""SOURCE ONLY: bounded catalogue-power outer packing and full balanced pricing.

One catalogue base at a time, with LITA32/Strassen fillers; up to 64 packets.
Exact positives, no global minimality or numerical-screen exclusion claims.
"""
import argparse
import bisect
import hashlib
import importlib.util
import json
import math
import resource
import sys
import time
from pathlib import Path

START = time.process_time()
sys.set_int_max_str_digits(0)
FULL_SHA = '32c7eea410bcd5f3d1e9f039066a8081f7404427a6889c846cc0f50964189a9e'
HELPER_SHA = 'f5b30a56c675c6698bf85ec0ff2abcbe7fba1f3e0bff94ea9d6bc9248066df92'
CACHED_SHA = '0548ffa84ee20c7270f5635abfac7fa3e5f46b1368ac065cec11ffc60f7578be'


def load_module(path, sha, name):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == sha
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bound_json(path, sha):
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == sha
    return json.loads(raw)


def save(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(path)


def outer_groups(capacity, bases, max_batches, budget):
    """Keep lowest rank at each side, then partition by exact packet count."""
    maxrank = capacity*max_batches
    choices = {}
    partial = False

    def fillers(base, h, rank_h, side_h):
        nonlocal partial
        i, rank_i, side_i = 0, rank_h, side_h
        while rank_i <= maxrank:
            rank, side = rank_i, side_i
            for j in range(5):
                if rank <= maxrank:
                    budget['candidates'] += 1
                    if budget['candidates'] > budget['limit'] or (budget['candidates'] % 4096 == 0 and time.process_time()-START > 17.0):
                        partial = True
                        return
                    item = (rank, base['root_object'], h, i, j, int(base['n']), int(base['rank']))
                    if side not in choices or rank < choices[side][0]:
                        choices[side] = item
                rank *= 7
                side *= 2
            i += 1
            rank_i *= 14196
            side_i *= 32

    fillers(dict(root_object=None, n=1, rank=1), 0, 1, 1)
    for base in bases:
        if partial:
            break
        n, rank_base = int(base['n']), int(base['rank'])
        h, rank_h, side_h = 1, rank_base, n
        while rank_h <= maxrank:
            fillers(base, h, rank_h, side_h)
            if partial:
                break
            h += 1
            rank_h *= rank_base
            side_h *= n
    groups = {}
    for side, rest in choices.items():
        batches = (rest[0]+capacity-1)//capacity
        assert 1 <= batches <= max_batches
        groups.setdefault(batches, []).append((side,)+rest)
    for entries in groups.values():
        entries.sort()
    return groups, partial, len(choices)


def outer_record(entry):
    side, rank, root, h, i, j, base_n, base_rank = entry
    assert side == base_n**h*32**i*2**j
    assert rank == base_rank**h*14196**i*7**j
    return dict(side=side, rank=rank, public_base_copies=h, public_root=root,
                public_base_n=base_n, public_base_rank=base_rank,
                lita32_copies=i, strassen2_copies=j)


def packing_controls():
    # Exhaustively compare the grouped family with explicit small exponents.
    bases = [dict(n=3, rank=23, root_object='test3'), dict(n=4, rank=49, root_object='test4')]
    budget = dict(candidates=0, limit=100000)
    groups, partial, _ = outer_groups(100, bases, 3, budget)
    assert not partial
    expected = {}
    for base in [dict(n=1, rank=1, root_object=None)]+bases:
        for h in range(0 if base['n'] == 1 else 1, 1 if base['n'] == 1 else 5):
            for i in range(3):
                for j in range(5):
                    rank = base['rank']**h*14196**i*7**j
                    side = base['n']**h*32**i*2**j
                    if rank <= 300:
                        expected[side] = min(rank, expected.get(side, rank))
    actual = {entry[0]: entry[1] for entries in groups.values() for entry in entries}
    assert actual == expected
    for batches, entries in groups.items():
        assert all((e[1]+99)//100 == batches for e in entries)
        assert all(outer_record(e)['rank'] == e[1] for e in entries)
    return 'PASS: exact bounded-family enumeration, side deduplication, ceil batch grouping'


def main():
    ap = argparse.ArgumentParser()
    for name in ('full-source', 'helper', 'cached', 'grading', 'catalogue', 'out'):
        ap.add_argument('--'+name, type=Path, required=True)
    ap.add_argument('--max-q', type=int, default=64)
    ap.add_argument('--max-t', type=int, default=64)
    ap.add_argument('--max-batches', type=int, default=64)
    ap.add_argument('--max-outer-candidates', type=int, default=10000000)
    args = ap.parse_args()
    assert 2 <= args.max_q <= 128 and 1 <= args.max_t <= 128 and 1 <= args.max_batches <= 64
    assert 1 <= args.max_outer_candidates <= 10000000
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
    full = load_module(args.full_source, FULL_SHA, 'reviewed_full')
    helper = load_module(args.helper, HELPER_SHA, 'reviewed_helper')
    cached = bound_json(args.cached, CACHED_SHA)
    grading = bound_json(args.grading, helper.GRADING_SHA)
    catalogue = bound_json(args.catalogue, helper.CATALOGUE_SHA)
    bases = [x for x in catalogue if x['field'] == 'Q' and int(x['n']) > 1]
    assert all(int(x['rank']) > 1 and isinstance(x['root_object'], str) and len(x['root_object']) == 64 for x in bases)
    # Most efficient bases first makes a CPU-truncated family useful quickly.
    bases.sort(key=lambda x: math.log(int(x['rank']))/math.log(int(x['n'])))
    lita = next(x for x in bases if int(x['n']) == 32)
    assert int(lita['rank']) == 14196 and 14196 < 7**5
    args.out.mkdir(parents=True, exist_ok=True)
    result = dict(status='RUNNING', source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  full_source_sha256=FULL_SHA, helper_sha256=HELPER_SHA, cached_sha256=CACHED_SHA,
                  grading_sha256=helper.GRADING_SHA, catalogue_sha256=helper.CATALOGUE_SHA,
                  catalogue_bases=len(bases), max_t=args.max_t, max_q=args.max_q, max_batches=args.max_batches,
                  scope='Each single catalogue base power with LITA32/Strassen fillers; up to64 shared packets; exact positives; bounded heuristic pricing, no global minimum',
                  limits=dict(case_bits=128000000, total_bits=2000000000, outer_candidates=args.max_outer_candidates),
                  cumulative_new_encoded_bits=0, outer_candidates=0, cases=[],
                  best_natural=cached['best_natural'], best_padded=cached['best_padded'], first_new_hit=None)

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
    result['controls'] = [packing_controls(), full.tiny_controls(flint.fmpz), helper.controls()]
    cache = {(tuple(p['puzzle']), c['t']): c['capacity'] for p in cached['puzzles'] for c in p['cases']}
    puzzles = {}
    queue = []
    for target in grading['targets']:
        s, w, d, edges, histogram = helper.replay(target)
        puzzles[s, w] = (d, edges, histogram)
        for t in range(1, args.max_t+1):
            # q2, one packet: rigorous lower price and useful processing priority.
            lower_rank = (3*t*d+1)*3**(3*w*t)
            upper = max(result['best_natural']['natural_n'], result['best_padded']['padded_n'])
            if lower_rank**20 >= upper**53:
                break
            queue.append((lower_rank, s, w, t))
    queue.sort()
    result['queued_cases'] = len(queue)
    budget = dict(candidates=0, limit=args.max_outer_candidates)
    for _, s, w, t in queue:
        if time.process_time()-START > 17.0:
            result['status'] = 'PARTIAL_CPU_BUDGET'
            break
        d, edges, histogram = puzzles[s, w]
        lower_rank = (3*t*d+1)*3**(3*w*t)
        upper = max(result['best_natural']['natural_n'], result['best_padded']['padded_n'])
        if lower_rank**20 >= upper**53:
            continue
        record = dict(puzzle=[s, w], t=t, degree=d, grading_replay='PASS', structural_edges=edges, status='PENDING')
        result['cases'].append(record)
        checkpoint()
        if ((s, w), t) in cache:
            capacity = cache[((s, w), t)]
            record['capacity_source'] = 'hash-bound full-balanced prior result'
        else:
            bits = (s**(3*t)).bit_length()*(3*w*t+1)**2
            if bits > 128000000 or result['cumulative_new_encoded_bits']+bits > 2000000000:
                record['status'] = 'UNKNOWN_BIT_BUDGET'
                continue
            capacity, encoding = full.encoded_capacity(histogram, w, t, flint.fmpz)
            assert encoding['encoded_bit_bound'] == bits
            result['cumulative_new_encoded_bits'] += bits
            record['encoding'] = encoding
            record['capacity_source'] = 'exact carry-free Kronecker extraction'
        record['capacity'] = capacity
        groups, partial, count = outer_groups(capacity, bases, args.max_batches, budget)
        result['outer_candidates'] = budget['candidates']
        record['outer_choices'] = count
        record['outer_family_complete'] = not partial
        factor = 3*t*d+1
        stop = False
        for q in range(2, args.max_q+1):
            base_logrank = math.log(factor)+3*w*t*math.log(q+1)
            inner_log = w*t*math.log(q)
            base_rank = inner = None
            for batches, entries in sorted(groups.items()):
                if time.process_time()-START > 18.0:
                    stop = True
                    break
                logrank = base_logrank+math.log(batches)
                if logrank >= 2.65*(inner_log+math.log(entries[-1][0]))+1e-10:
                    continue
                bound = max(math.log(result['best_natural']['natural_n']), math.log(result['best_padded']['padded_n']))
                if logrank/2.65 > bound+1e-10:
                    continue
                if base_rank is None:
                    base_rank = factor*(q+1)**(3*w*t)
                    inner = q**(w*t)
                rank = batches*base_rank
                rank20 = rank**20
                if rank20 >= (entries[-1][0]*inner)**53:
                    continue
                padded = helper.integer_root(rank20, 53)+1
                needed_side = (padded+inner-1)//inner
                index = bisect.bisect_left(entries, (needed_side,))
                chosen = outer_record(entries[index])
                natural = chosen['side']*inner
                assert batches == (chosen['rank']+capacity-1)//capacity
                assert rank20 < natural**53 and (padded-1)**53 <= rank20 < padded**53
                item = dict(puzzle=[s, w], degree=d, t=t, q=q, base_copies=3*t, power=3*w*t,
                            square_capacity=capacity, batches=batches, inner_side=inner, outer=chosen,
                            interpolation_factor=factor, rank=rank, natural_n=natural, padded_n=padded,
                            natural_exponent=math.log(rank)/math.log(natural), exact_threshold='PASS',
                            construction='Full balanced cyclic packet; ceil outer rank/capacity packets; arbitrary pinned catalogue-base power with LITA32/Strassen fillers')
                if result['first_new_hit'] is None:
                    result['first_new_hit'] = item
                if (natural, rank) < (result['best_natural']['natural_n'], result['best_natural']['rank']):
                    result['best_natural'] = item
                if (padded, rank) < (result['best_padded']['padded_n'], result['best_padded']['rank']):
                    result['best_padded'] = item
                checkpoint()
            if stop:
                break
        record['status'] = 'PARTIAL' if partial or stop else 'COMPLETE_BOUNDED_SCREEN'
        checkpoint()
        if partial or stop:
            result['status'] = 'PARTIAL_OUTER_OR_CPU_BUDGET'
            break
    else:
        result['status'] = ('PARTIAL_BIT_BUDGET' if any(x['status'] == 'UNKNOWN_BIT_BUDGET' for x in result['cases'])
                            else 'COMPLETE_BOUNDED_SCREEN')
    checkpoint()


if __name__ == '__main__':
    main()
