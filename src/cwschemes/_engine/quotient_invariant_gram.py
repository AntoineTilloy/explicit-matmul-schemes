"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from functools import lru_cache
from itertools import product
import json
import math
import numpy as np
from .tailored_completion import layout

def local_bases(q, width):
    """Refine OO into diagonal and off-diagonal transitive label orbits."""
    assert q >= 4
    result = []
    v = np.zeros(q, dtype=np.int64)
    v[0] = 1
    v[1] = -1
    states = range(5) if width == 2 else range(2 + q % 2)
    for state in states:
        if width == 1:
            paired = q - q % 2
            labels = np.array([[i + 1] for i in range(paired)] if state == 1 else [[q]] if state == 2 else [[0]])
            choices = [(0, 0, np.ones(len(labels), dtype=np.int64))]
            if state == 1:
                even = np.zeros(paired, dtype=np.int64)
                even[:4] = [1, 1, -1, -1]
                odd = v[:paired].copy()
                choices += [(1, 0, even), (2, 0, odd)]
        else:
            ls = {0: [(0, 0)], 1: [(i + 1, 0) for i in range(q)], 2: [(0, i + 1) for i in range(q)], 3: [(i + 1, i + 1) for i in range(q)], 4: [(i + 1, j + 1) for i in range(q) for j in range(q) if i != j]}
            labels = np.array(ls[state], dtype=np.int64)
            choices = [(0, 0, np.ones(len(labels), dtype=np.int64))]
            if state in (1, 2, 3):
                choices += [(1, 0, v.copy())]
            if state == 4:
                sym = np.zeros((q, q), dtype=np.int64)
                for i, j, value in [(0, 1, 1), (2, 3, 1), (0, 2, -1), (1, 3, -1)]:
                    sym[i, j] = sym[j, i] = value
                alt = np.zeros((q, q), dtype=np.int64)
                for i, j in [(0, 1), (1, 2), (2, 0)]:
                    alt[i, j] = 1
                    alt[j, i] = -1
                choices += [(1, 0, v[labels[:, 0] - 1]), (1, 1, v[labels[:, 1] - 1]), (2, 0, sym[labels[:, 0] - 1, labels[:, 1] - 1]), (3, 0, alt[labels[:, 0] - 1, labels[:, 1] - 1])]
        for irrep, copy, vector in choices:
            result.append(dict(state=state, irrep=irrep, copy=copy, labels=labels, vector=vector))
    index = {(r['state'], r['irrep'], r['copy']): i for i, r in enumerate(result)}
    for r in result:
        state, irrep, copy = (r['state'], r['irrep'], r['copy'])
        gs = {1: 2, 2: 1}.get(state, state) if width == 2 else state
        gc = 1 - copy if width == 2 and state == 4 and (irrep == 1) else copy
        r['image'] = index[gs, irrep, gc]
        r['sign'] = -1 if width == 2 and irrep == 3 or (width == 1 and irrep == 2) else 1
    return result

def state_options(ds, q=None):
    if len(ds) == 1:
        return (1, 2) if ds[0] == 1 and q is not None and q % 2 else (int(ds[0] == 1),)
    return {(2, 2): (0,), (1, 2): (1,), (2, 1): (2,), (1, 1): (3, 4)}[tuple(ds)]

def representation(q, rows, m, groups, domain_character=-1):
    spec = layout(tuple(map(tuple, rows)), m)
    assert not spec['common']
    assert sorted((i for group in groups for i in group)) == list(range(3 * m))
    bases = [local_bases(q, len(group)) for group in groups]
    for group in groups:
        assert len({spec['missing'][i] for i in group}) == 1
    states = sorted({st for ds in spec['digits'].values() for st in product(*(state_options([ds[i] for i in g], q) for g in groups))})
    orbit_list = []
    orbit_index = {}
    for st in states:
        gst = tuple(({1: 2, 2: 1}.get(s, s) if len(g) == 2 else s for g, s in zip(groups, st)))
        assert gst in states
        if st > gst:
            continue
        raw_size = math.prod(((1, q, q, q, q * (q - 1))[s] if len(g) == 2 else (1, q - q % 2, 1)[s] for g, s in zip(groups, st)))
        fixed = st == gst and all((s in (0, 3) if len(g) == 2 else s != 1 for g, s in zip(groups, st)))
        size = 0 if fixed else raw_size // 2 if st == gst else raw_size
        if not size:
            continue
        degree = sum(((0, 1, 1, 2, 2)[s] if len(g) == 2 else int(s != 0) for g, s in zip(groups, st)))
        oid = len(orbit_list)
        orbit_list.append(dict(states=st, image=gst, size=size, degree=degree))
        orbit_index[st] = orbit_index[gst] = oid
    block_list = []
    entries = []
    for irreps in product(*(range(4) if len(g) == 2 else range(3) for g in groups)):
        records = []
        for st in states:
            if st not in orbit_index:
                continue
            choices = [[i for i, b in enumerate(bs) if b['state'] == s and b['irrep'] == ir] for bs, s, ir in zip(bases, st, irreps)]
            for indices in product(*choices):
                image = tuple((bs[i]['image'] for bs, i in zip(bases, indices)))
                sign = math.prod((bs[i]['sign'] for bs, i in zip(bases, indices)))
                if indices > image or (indices == image and sign != domain_character):
                    continue
                records.append(dict(indices=indices, image=image, sign=sign, factor=2 if indices != image else 1, orbit=orbit_index[st]))
        if records:
            dimension = math.prod(((1, q - 1, q * (q - 3) // 2, (q - 1) * (q - 2) // 2)[ir] if len(g) == 2 else (1, q // 2 - 1, q // 2)[ir] for g, ir in zip(groups, irreps)))
            block_list.append(dict(irreps=irreps, dimension=dimension, records=records))
            entries.extend((r['indices'] for r in records))
    assert sum((b['dimension'] * len(b['records']) for b in block_list)) == sum((o['size'] for o in orbit_list))
    return (spec, bases, orbit_list, block_list)

def coordinate_kernel(q, left, right, missing, mode, zero, twist, width):
    """Original factor inner product in one coordinate group, before projection."""
    right = right.copy()
    if twist:
        if width == 2:
            right = right[:, ::-1]
        else:
            x = right[:, 0]
            right[:, 0] = np.where((x == 0) | (q % 2 == 1) & (x == q), x, np.where(x % 2, x + 1, x - 1))
    K = np.ones((len(left), len(right)), dtype=np.int64)
    for i in range(width):
        x = left[:, i, None]
        y = right[None, :, i]
        special = (x == 0) | (y == 0)
        if zero >> i & 1:
            if mode == missing:
                K *= ~special
        else:
            if mode != missing:
                K *= ~special
            K *= np.where((x != 0) & (y != 0), x == y, np.where((x == 0) & (y == 0), q, 1))
    return K

def local_tables(q, spec, bases, groups, pair, prime):
    """Each entry is an exact integer local contraction, then reduced mod p."""
    family_a, family_b = [spec['families'][mode] for mode in pair]
    sectors = [(za, zb, ta, tb) for za in family_a for zb in family_b for ta in range(2) for tb in range(2)]
    max_basis = max(map(len, bases))
    out = np.zeros((len(groups), max_basis, max_basis, len(sectors)), dtype=np.int64)
    for gi, (group, bs) in enumerate(zip(groups, bases)):
        width = len(group)
        missing = spec['missing'][group[0]]

        @lru_cache(None)
        def kernel(s, t, mode, z, twist):
            left = next((r['labels'] for r in bs if r['state'] == s))
            right = next((r['labels'] for r in bs if r['state'] == t))
            return coordinate_kernel(q, left, right, missing, mode, z, twist, width)

        @lru_cache(None)
        def contract(i, j, za, zb, ta, tb):
            l, r = (bs[i], bs[j])
            a = kernel(l['state'], r['state'], pair[0], za, ta)
            b = kernel(l['state'], r['state'], pair[1], zb, tb)
            value = int(l['vector'] @ (a * b) @ r['vector'])
            return value % prime if prime else value
        for i, left in enumerate(bs):
            for j, right in enumerate(bs):
                if left['irrep'] != right['irrep']:
                    continue
                for k, (za, zb, ta, tb) in enumerate(sectors):
                    a = sum(((za >> pos & 1) << bit for bit, pos in enumerate(group)))
                    b = sum(((zb >> pos & 1) << bit for bit, pos in enumerate(group)))
                    out[gi, i, j, k] = contract(i, j, a, b, ta, tb)
        print('local group', group, 'nonzeros', int(np.count_nonzero(out[gi])), flush=True)
    return (out, sectors)

def assemble(tables, blocks, pair, prime):
    if not blocks:
        return []
    from numba import njit, prange
    offsets = np.array([0] + list(np.cumsum([len(b['records']) for b in blocks])), dtype=np.int64)
    indices = np.array([r['indices'] for b in blocks for r in b['records']], dtype=np.int64)
    factors = np.array([r['factor'] for b in blocks for r in b['records']], dtype=np.int64)
    matrix_offsets = np.array([0] + list(np.cumsum([len(b['records']) ** 2 for b in blocks])), dtype=np.int64)
    characters = (-1, -1, 1)
    signs = np.array([characters[pair[0]] ** ta * characters[pair[1]] ** tb for ta in range(2) for tb in range(2)], dtype=np.int64)

    @njit(parallel=True)
    def work(tables, offsets, indices, factors, matrix_offsets, signs, prime):
        output = np.zeros(matrix_offsets[-1], dtype=np.int64)
        for block in prange(len(offsets) - 1):
            lo, hi = (offsets[block], offsets[block + 1])
            size = hi - lo
            for i in range(lo, hi):
                for j in range(i, hi):
                    total = 0
                    for sector in range(tables.shape[3]):
                        value = 1
                        for group in range(tables.shape[0]):
                            x = tables[group, indices[i, group], indices[j, group], sector]
                            if x == 0:
                                value = 0
                                break
                            value = value * x % prime if prime else value * x
                        total += signs[sector % 4] * value
                    total = total % prime * factors[i] * factors[j] % prime if prime else total * factors[i] * factors[j]
                    output[matrix_offsets[block] + (i - lo) * size + j - lo] = total
                    output[matrix_offsets[block] + (j - lo) * size + i - lo] = total
        return output
    packed = work(tables, offsets, indices, factors, matrix_offsets, signs, prime)
    return [packed[matrix_offsets[i]:matrix_offsets[i + 1]].reshape(len(b['records']), -1) for i, b in enumerate(blocks)]

def standard_problem():
    rows = json.loads((resource_root() / 'support14.json').read_text())['rows']
    groups = [(0, 1), (2,), (3, 4), (5,), (6, 7), (8,)]
    return (rows, 3, groups)
