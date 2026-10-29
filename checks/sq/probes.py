#!/usr/bin/env python3
"""Nonzero coefficient probes of the S_q packet entries through the published loader.

Each probe picks an outer call, a nonzero outer coefficient and a packet term whose local
factors are nonzero on that call's row, then compares Scheme.coefficient with a separate
evaluation: row wiring and local S_q factors re-derived here from docs/construction.md,
outer coefficients read from the explicit size-32 arrays or the current-LITA term forms.
Usage, from the repository root:
    PYTHONPATH=src python checks/sq/probes.py [--out checks/sq/probes.json]
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import random
import time

import cwschemes


def wiring(row, mode, r, c, q, m):
    """Variable: 0 on the leg's zero mask; row digits on B (U, W) or C (V); column digits on C, A, A."""
    A, B, C = row
    rmask, cmask = ((B, C), (C, A), (B, A))[mode]
    x = [0] * (3 * m)
    for mask, index in ((rmask, r), (cmask, c)):
        positions = [i for i in range(3 * m) if mask >> i & 1]
        digits = []
        for _ in positions:
            index, d = divmod(index, q)
            digits.append(d + 1)
        for i, d in zip(reversed(positions), digits):
            x[i] = d
    return x


def local(q, kind, mode, choice, x):
    """S(a,b,c) = sum_t (a0+at)(b0+bt)ct - a0 b0 sum_t ct; kind k reads a, b, c from legs k, k+1, k+2."""
    argument = 'abc'[(mode - kind) % 3]
    if choice < q:
        t = choice + 1
        return int(x == t) if argument == 'c' else int(x in (0, t))
    sign = -1 if mode == 2 else 1
    return sign * (int(x >= 1) if argument == 'c' else int(x == 0))


def outer_value(engine, mode, call, r, c):
    child = getattr(engine.outer, 'child', None)
    if hasattr(child, 'arrays'):
        a, axis = child.arrays, 'uvw'[mode]
        start, end = int(a[axis + '_indptr'][call]), int(a[axis + '_indptr'][call + 1])
        for k in range(start, end):
            if int(a[axis + '_indices'][k]) == child.n * r + c:
                return Fraction(int(a[axis + '_numerators'][k]), int(a[axis + '_denominators'][k]))
        return Fraction(0)
    return child.root.lita.term(call)[mode].get((r, c), Fraction(0))


def nonzero_positions(engine, mode, call):
    child = engine.outer.child
    if hasattr(child, 'arrays'):
        a, axis = child.arrays, 'uvw'[mode]
        start, end = int(a[axis + '_indptr'][call]), int(a[axis + '_indptr'][call + 1])
        return [divmod(int(a[axis + '_indices'][k]), child.n) for k in range(start, end)]
    return sorted(child.root.lita.term(call)[mode])


def probe(scheme, rng):
    engine = scheme._engine
    b = rng.randrange(len(engine.packets))
    p = engine.packets[b]
    q, m, s = p.q, p.m, p.q ** p.m
    h = rng.randrange(engine.used[b])
    mode = rng.randrange(3)
    call = engine.starts[b] + h
    while True:
        br, bc = rng.choice(nonzero_positions(engine, mode, call))
        ir, ic = rng.randrange(s), rng.randrange(s)
        if br * s + ir < scheme.n and bc * s + ic < scheme.n:
            break
    x = wiring(p.rows[h], mode, ir, ic, q, m)
    # A label t needs choice t - 1; a scalar index (only on legs a, b of S) survives every choice.
    choices = [v - 1 if v else rng.randrange(q + 1) for v in x]
    assert all(local(q, i // m, mode, c, v) for i, (c, v) in enumerate(zip(choices, x)))
    t = 0
    for c in choices:
        t = t * (q + 1) + c
    term = (engine.ends[b - 1] if b else 0) + t
    expected = Fraction(0)
    for g in range(engine.used[b]):
        y = wiring(p.rows[g], mode, ir, ic, q, m)
        value = 1
        for i in range(3 * m):
            value *= local(q, i // m, mode, choices[i], y[i])
            if not value:
                break
        if value:
            expected += outer_value(engine, mode, engine.starts[b] + g, br, bc) * value
    got = scheme.coefficient('UVW'[mode], term, br * s + ir, bc * s + ic)
    return dict(factor='UVW'[mode], term=str(term), row=br * s + ir, col=bc * s + ic, batch=b, packet_row=h,
                value=str(got), independent=str(expected), ok=got == expected and got != 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--probes', type=int, default=8)
    args = parser.parse_args()
    out = []
    for n in (8388608, 13947137604, 13468840704):
        start = time.monotonic()
        scheme = cwschemes.load(n)
        rng = random.Random(n)
        records = [probe(scheme, rng) for _ in range(args.probes)]
        assert all(r['ok'] for r in records), records
        out.append(dict(n=str(n), rank=str(scheme.rank), probes=records, passed=True, seconds=time.monotonic() - start))
    record = dict(status='PASS', script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), entries=out)
    if args.out:
        args.out.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(dict(status='PASS', probes=sum(len(e['probes']) for e in out))))


if __name__ == '__main__':
    main()
