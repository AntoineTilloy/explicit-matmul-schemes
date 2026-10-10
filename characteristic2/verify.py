#!/usr/bin/env python3
"""Independent exact primitive, Kronecker and masked-tensor checks over F2."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from crop_accessor import FrozenBasis
from class_recompress_accessor import ClassRecompressed


def basis(path):
    return FrozenBasis(path, hashlib.sha256(path.read_bytes()).hexdigest())


def support(b, mode, t):
    entries = list(b.sparse(mode, t))
    assert all(int(num) == int(den) == 1 for _, num, den in entries)
    return np.asarray([int(k) for k, _, _ in entries], dtype=np.int64)


def tensor(forms):
    """All input-coordinate pairs, with every output coefficient packed."""
    result = {}
    for u, v, w in forms:
        aa = []
        while u:
            low = u & -u
            aa.append(low.bit_length() - 1)
            u -= low
        while v:
            low = v & -v
            b = low.bit_length() - 1
            v -= low
            for a in aa:
                key = (a, b)
                result[key] = result.get(key, 0) ^ w
    return {k: v for k, v in result.items() if v}


def full_primitive(b):
    forms = [tuple(sum(1 << int(k) for k in support(b, mode, t))
                   for mode in range(3)) for t in range(b.rank)]
    expected = {(i*b.n+k, k*b.n+j): 1 << (i*b.n+j)
                for i in range(b.n) for j in range(b.n) for k in range(b.n)}
    assert tensor(forms) == expected
    return b.n ** 4


def kronecker(child, left, right):
    assert child.n == left.n * right.n and child.rank == left.rank * right.rank
    checked = 0
    for mode in range(3):
        ls = [support(left, mode, t) for t in range(left.rank)]
        rs = [support(right, mode, t) for t in range(right.rank)]
        for t in range(child.rank):
            a, b = divmod(t, right.rank)
            ar, ac = np.divmod(ls[a], left.n)
            br, bc = np.divmod(rs[b], right.n)
            coords = ((ar[:, None]*right.n+br)*child.n
                      + ac[:, None]*right.n+bc).ravel()
            got = support(child, mode, t)
            assert np.array_equal(np.sort(coords), np.sort(got))
            checked += len(got)
    return checked


def before_forms(c, key):
    """Rebuild clipped outer factors directly from frozen sparse coordinates."""
    blocks = c.n // c.base.inner.n
    orders = c.base.spec.get('vertex_orders')
    if orders is None:
        orders = [c.base.spec['outer_order']] * 3
    positions = [{original: new for new, original in enumerate(order)}
                 for order in orders]
    edges = ((0, 1), (1, 2), (0, 2))
    forms = []
    for t in range(c.base.outer.rank):
        factors = []
        for mode in range(3):
            flags = (key >> (3*mode)) & 7
            mask = 0
            for index in support(c.base.outer, mode, t):
                r, col = divmod(int(index), c.base.outer.n)
                a, b = edges[mode]
                r, col = positions[a][r], positions[b][col]
                alive = ((r < blocks and col < blocks)
                         or (r == blocks and col < blocks and flags & 1)
                         or (r < blocks and col == blocks and flags & 2)
                         or (r == blocks and col == blocks and flags & 4))
                if alive:
                    mask |= 1 << (r*c.base.outer.n+col)
            factors.append(mask)
        forms.append(tuple(factors))
    return forms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    start = time.monotonic()
    b4, b16, b64 = [basis(root/'bases'/f'{n}at.npz') for n in (4, 16, 64)]
    primitive_pairs = full_primitive(b4)
    coefficient_checks = dict(
        size16=kronecker(b16, b4, b4), size64=kronecker(b64, b16, b4))
    results = []
    for n in (61, 243, 969):
        c = ClassRecompressed(root/f'F2_{n}_class_recompressed.json')
        checked = 0
        for j, key in enumerate(c.keys):
            old = tensor(before_forms(c, key))
            new = tensor(c._outer(j, int(t)) for t in c.outer_ids[j])
            assert old == new
            checked += 1
        row = dict(n=n, field='F2', rank=c.rank, masked_outer_tensors=checked,
                   complete_allocations=True, exact=True)
        results.append(row)
        print(json.dumps(row), flush=True)
    result = dict(exact=True, primitive_input_pairs=primitive_pairs,
                  kronecker_coefficient_checks=coefficient_checks,
                  points=results, seconds=time.monotonic()-start)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
