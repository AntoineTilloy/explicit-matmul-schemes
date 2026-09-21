"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
import math
import numpy as np
from flint import nmod_mat
from .quotient_invariant_gram import representation, coordinate_kernel, standard_problem
from .quotient_coordinates import shell

def mod_fraction(value, prime):
    return int(value.numerator) * pow(int(value.denominator), -1, prime) % prime

def orbital_key(width, left, right):
    if width == 2:
        seen = {}
        key = []
        for x in left + right:
            if not x:
                key.append(0)
            else:
                if x not in seen:
                    seen[x] = len(seen) + 1
                key.append(seen[x])
        return tuple(key)
    x, y = (left[0], right[0])
    return (0 if not x or not y else 1 if x == y else 2 if (x - 1) // 2 == (y - 1) // 2 else 3,)

def local_transform(instance, group, prime):
    bs = instance.bases[group]
    width = len(instance.groups[group])
    q = instance.q
    states = sorted({b['state'] for b in bs})
    channels = []
    pieces = []
    pair_index = np.full((len(bs), len(bs)), -1, dtype=np.int64)
    for left_state in states:
        for right_state in states:
            pairs = [(i, j) for i, a in enumerate(bs) for j, b in enumerate(bs) if a['state'] == left_state and b['state'] == right_state and (a['irrep'] == b['irrep'])]
            ls = next((b['labels'] for b in bs if b['state'] == left_state))
            rs = next((b['labels'] for b in bs if b['state'] == right_state))
            orbitals = {}
            for l in ls:
                for r in rs:
                    l, r = (tuple(map(int, l)), tuple(map(int, r)))
                    orbitals.setdefault(orbital_key(width, l, r), (l, r))
            orbitals = list(orbitals.values())
            assert len(pairs) == len(orbitals)
            start = len(channels)
            matrix = np.zeros((len(pairs), len(pairs)), dtype=np.int64)
            for row, (l, r) in enumerate(orbitals):
                for col, (i, j) in enumerate(pairs):
                    irrep = bs[i]['irrep']
                    dimension = (1, q - 1, q * (q - 3) // 2, (q - 1) * (q - 2) // 2)[irrep] if width == 2 else (1, q // 2 - 1, q // 2)[irrep]
                    matrix[row, col] = mod_fraction(dimension * instance.local_moment(group, i, j, l, r), prime)
            inverse = nmod_mat(matrix.tolist(), prime).inv()
            inverse = np.array([[int(inverse[i, j]) for j in range(len(pairs))] for i in range(len(pairs))], dtype=np.int64)
            pieces.append((start, len(pairs), matrix, inverse))
            for k, ((i, j), (l, r)) in enumerate(zip(pairs, orbitals)):
                pair_index[i, j] = start + k
                channels.append(dict(left_state=left_state, right_state=right_state, left=l, right=r, basis_pair=(i, j)))
    assert len(channels) == (52 if width == 2 else 11 if q % 2 else 6)
    return dict(channels=channels, pieces=pieces, pair_index=pair_index)

class HadamardEngine:

    def __init__(self, q, rows, m, groups, prime=65521, threads=1):
        from numba import set_num_threads
        set_num_threads(threads)
        self.q = q
        self.rows = rows
        self.m = m
        self.groups = groups
        self.prime = prime
        self.instance, negative, self.orbits = shell(q, rows, m, groups, -1)
        self.spec, self.bases, _, positive = representation(q, rows, m, groups, 1)
        self.quotients = {-1: negative, 1: positive}
        self.transforms = [local_transform(self.instance, i, prime) for i in range(len(groups))]
        self.shape = tuple((len(t['channels']) for t in self.transforms))
        by_irrep = {}
        for character, blocks in self.quotients.items():
            for block in blocks:
                key = tuple(block['irreps'])
                record = by_irrep.setdefault(key, dict(irreps=key, full=set(), quotient={}))
                record['quotient'][character] = block
                for r in block['records']:
                    record['full'].add(tuple(r['indices']))
                    record['full'].add(tuple(r['image']))
        self.full_blocks = []
        for key in sorted(by_irrep):
            b = by_irrep[key]
            full = sorted(b['full'])
            lookup = {indices: i for i, indices in enumerate(full)}
            indices = np.array(full, dtype=np.int64)
            size = len(full)
            positions = np.zeros((size, size), dtype=np.int64)
            for g, transform in enumerate(self.transforms):
                positions = positions * self.shape[g] + transform['pair_index'][indices[:, g, None], indices[None, :, g]]
            assert positions.min() >= 0
            state_groups = {}
            for i, ids in enumerate(full):
                states = tuple((bs[k]['state'] for bs, k in zip(self.bases, ids)))
                state_groups.setdefault(states, []).append(i)
            metrics = []
            for ids in state_groups.values():
                M = np.ones((len(ids), len(ids)), dtype=np.int64)
                for g, bs in enumerate(self.bases):
                    vectors = np.array([bs[full[i][g]]['vector'] for i in ids])
                    M = M * (vectors @ vectors.T) % prime
                metrics.append((np.array(ids), M))
            quotients = {}
            for character, block in b['quotient'].items():
                records = block['records']
                map_to = np.full(size, -1, dtype=np.int64)
                signs = np.zeros(size, dtype=np.int64)
                representatives = []
                orbit_groups = {}
                for i, r in enumerate(records):
                    x, y = (lookup[tuple(r['indices'])], lookup[tuple(r['image'])])
                    representatives.append(x)
                    map_to[x] = map_to[y] = i
                    signs[x] = 1
                    signs[y] = character * r['sign']
                    orbit_groups.setdefault(r['orbit'], []).append(i)
                self.instance.domain_character = character
                qmetrics = []
                for ids in orbit_groups.values():
                    M = np.array([[self.instance.metric(records[i], records[j]) for j in ids] for i in ids], dtype=np.int64) % prime
                    inverse = nmod_mat(M.tolist(), prime).inv()
                    inverse = np.array([[int(inverse[i, j]) for j in range(len(ids))] for i in range(len(ids))], dtype=np.int64)
                    qmetrics.append((np.array(ids), M, inverse))
                quotients[character] = dict(records=records, representatives=np.array(representatives), factors=np.array([r['factor'] for r in records]), map_to=map_to, signs=signs, metrics=qmetrics, dimension=block['dimension'])
            self.full_blocks.append(dict(irreps=key, indices=indices, positions=positions, metrics=metrics, quotients=quotients))
        self.instance.domain_character = -1
        self.deleted = set()
        widths = [piece[1] for t in self.transforms for piece in t['pieces']]
        widths += [len(ids) for b in self.full_blocks for ids, _ in b['metrics']]
        widths += [len(item[0]) for b in self.full_blocks for qb in b['quotients'].values() for item in qb['metrics']]
        self.maximum_metric_dot_width = max(widths)
        assert self.maximum_metric_dot_width * (prime - 1) ** 2 <= np.iinfo(np.int64).max, 'modular dot product can overflow int64'
        assert 2 * len(rows) * q ** (3 * m) <= np.iinfo(np.int64).max, 'initial integer kernel can overflow int64'
        print('Hadamard geometry', self.shape, math.prod(self.shape), 'full blocks', len(self.full_blocks), flush=True)

    def transform(self, values, inverse=False):
        values = values.reshape(self.shape)
        for axis, transform in enumerate(self.transforms):
            moved = np.moveaxis(values, axis, 0)
            shape = moved.shape
            flat = moved.reshape(shape[0], -1)
            out = np.zeros_like(flat)
            for start, size, forward, backward in transform['pieces']:
                matrix = backward if inverse else forward
                out[start:start + size] = matrix @ flat[start:start + size] % self.prime
            values = np.moveaxis(out.reshape(shape), 0, axis)
        return np.ascontiguousarray(values).ravel()

    @staticmethod
    def bilateral(matrix, metrics, inverse=False):
        matrix = matrix.copy()
        prime = HadamardEngine._active_prime
        for item in metrics:
            ids, M = (item[0], item[2] if inverse else item[1])
            matrix[ids, :] = M @ matrix[ids, :] % prime
        for item in metrics:
            ids, M = (item[0], item[2] if inverse else item[1])
            matrix[:, ids] = matrix[:, ids] @ M % prime
        return matrix

    def kernel_to_blocks(self, values, character):
        coefficients = self.transform(values, inverse=True)
        HadamardEngine._active_prime = self.prime
        matrices = []
        blocks = []
        for b in self.full_blocks:
            if character not in b['quotients']:
                continue
            qb = b['quotients'][character]
            full = self.bilateral(coefficients[b['positions']], b['metrics'])
            r = qb['representatives']
            f = qb['factors']
            matrix = full[np.ix_(r, r)] * f[:, None] % self.prime * f[None, :] % self.prime
            matrices.append(matrix)
            blocks.append(dict(irreps=b['irreps'], dimension=qb['dimension'], records=qb['records']))
        return (matrices, blocks)

    def blocks_to_kernel(self, matrices, character):
        coefficients = np.zeros(math.prod(self.shape), dtype=np.int64)
        HadamardEngine._active_prime = self.prime
        index = 0
        for b in self.full_blocks:
            if character not in b['quotients']:
                continue
            qb = b['quotients'][character]
            Q = self.bilateral(matrices[index], qb['metrics'], inverse=True)
            index += 1
            active = np.flatnonzero(qb['map_to'] >= 0)
            mapping = qb['map_to'][active]
            signs = qb['signs'][active]
            W = 4 * Q[np.ix_(mapping, mapping)] % self.prime * signs[:, None] % self.prime * signs[None, :] % self.prime
            coefficients[b['positions'][np.ix_(active, active)]] = W
        assert index == len(matrices)
        return self.transform(coefficients)

    def initial_kernels(self):
        from numba import njit, prange
        q, prime = (self.q, self.prime)
        dimensions = np.array(self.shape, dtype=np.int64)
        radices = np.array([5 if len(g) == 2 else 2 + self.q % 2 for g in self.groups], dtype=np.int64)
        channels = np.zeros((len(self.groups), max(self.shape), 3), dtype=np.int64)
        tables = np.zeros((3, len(self.groups), max(self.shape), len(self.rows), 2), dtype=np.int64)
        for g, (group, transform) in enumerate(zip(self.groups, self.transforms)):
            for i, channel in enumerate(transform['channels']):
                channels[g, i] = (channel['left_state'], channel['right_state'], sum((x == 0 for x in channel['left'] + channel['right'])) % 2)
                for mode in range(3):
                    for z, zero in enumerate(self.spec['families'][mode]):
                        local_zero = sum(((zero >> pos & 1) << j for j, pos in enumerate(group)))
                        for twist in range(2):
                            tables[mode, g, i, z, twist] = coordinate_kernel(q, np.array([channel['left']]), np.array([channel['right']]), self.spec['missing'][group[0]], mode, local_zero, twist, len(group))[0, 0]
        allowed = np.zeros(math.prod(radices), dtype=np.bool_)
        for orbit in self.orbits:
            for states in (orbit['states'], orbit['image']):
                code = 0
                for radix, state in zip(radices, states):
                    code = code * radix + state
                allowed[code] = True
        size = math.prod(self.shape)
        strides = np.array([math.prod(self.shape[g + 1:]) for g in range(len(self.shape))], dtype=np.int64)

        @njit(parallel=True)
        def build(size, dimensions, strides, radices, channels, tables, allowed, prime):
            out = np.zeros((3, size), dtype=np.int64)
            for index in prange(size):
                ids = np.zeros(len(dimensions), dtype=np.int64)
                for g in range(len(dimensions)):
                    ids[g] = index // strides[g] % dimensions[g]
                left = right = parity = 0
                for g in range(len(dimensions)):
                    left = left * radices[g] + channels[g, ids[g], 0]
                    right = right * radices[g] + channels[g, ids[g], 1]
                    parity += channels[g, ids[g], 2]
                if not allowed[left] or not allowed[right]:
                    continue
                for mode in range(3):
                    total = 0
                    for z in range(tables.shape[3]):
                        for twist in range(2):
                            value = 1
                            for g in range(len(dimensions)):
                                value *= tables[mode, g, ids[g], z, twist]
                                if not value:
                                    break
                            total += value * (-1 if mode != 2 and twist else 1)
                    out[mode, index] = 2 * total * (-1 if mode == 2 and parity % 2 else 1) % prime
            return out
        return build(size, dimensions, strides, radices, channels, tables, allowed, prime)

    def removal_mask(self, deleted):
        from numba import njit, prange
        radices = np.array([5 if len(g) == 2 else 2 + self.q % 2 for g in self.groups], dtype=np.int64)
        states = np.zeros((len(self.groups), max(self.shape), 2), dtype=np.int64)
        for g, transform in enumerate(self.transforms):
            for i, c in enumerate(transform['channels']):
                states[g, i] = (c['left_state'], c['right_state'])
        removed = np.zeros(math.prod(radices), dtype=np.bool_)
        for orbit in deleted:
            for word in (self.orbits[orbit]['states'], self.orbits[orbit]['image']):
                code = 0
                for radix, x in zip(radices, word):
                    code = code * radix + x
                removed[code] = True
        shape = np.array(self.shape)
        strides = np.array([math.prod(self.shape[g + 1:]) for g in range(len(self.shape))])

        @njit(parallel=True)
        def mask(size, shape, strides, radices, states, removed):
            result = np.zeros(size, dtype=np.bool_)
            for index in prange(size):
                a = b = 0
                for g in range(len(shape)):
                    c = index // strides[g] % shape[g]
                    a = a * radices[g] + states[g, c, 0]
                    b = b * radices[g] + states[g, c, 1]
                result[index] = removed[a] or removed[b]
            return result
        return mask(math.prod(self.shape), shape, strides, radices, states, removed)

def update_free(engine, kernels, free_mode, character, metadata, relations, removed):
    matrices, free_metadata = engine.kernel_to_blocks(kernels[free_mode], character)
    assert [b['irreps'] for b in free_metadata] == [b['irreps'] for b in metadata]
    geometry = [b['quotients'][character] for b in engine.full_blocks if character in b['quotients']]
    correction_entries = 0
    for i, (matrix, block, qb) in enumerate(zip(matrices, metadata, geometry)):
        records = block['records']
        if i in relations:
            relation = relations[i]
            J = relation['deleted']
            B = relation['basis']
            C = relation['coefficients']
            assert len(J) * (engine.prime - 1) ** 2 <= np.iinfo(np.int64).max, 'update dot product can overflow int64'
            T = np.zeros((len(J), len(records)), dtype=np.int64)
            T[:, B] = C.T
            for ids, M, inverse in qb['metrics']:
                T[:, ids] = T[:, ids] @ M % engine.prime
            lookup = {j: r for r, j in enumerate(J)}
            for ids, M, inverse in qb['metrics']:
                if int(ids[0]) in lookup:
                    rows = [lookup[int(j)] for j in ids]
                    T[rows, :] = inverse @ T[rows, :] % engine.prime
            cross = matrix[:, J] @ T % engine.prime
            correction = T.T @ matrix[np.ix_(J, J)] % engine.prime @ T % engine.prime
            correction_entries += int(np.count_nonzero(correction))
            matrix = (matrix + cross + cross.T + correction) % engine.prime
        zero = [j for j, r in enumerate(records) if r['orbit'] in removed or r['orbit'] in engine.deleted]
        matrix[zero, :] = 0
        matrix[:, zero] = 0
        matrices[i] = matrix
    updated = engine.blocks_to_kernel(matrices, character)
    mask = engine.removal_mask(removed)
    assert not np.any(updated[mask])
    for mode in range(3):
        if mode == free_mode:
            kernels[mode] = updated
        else:
            kernels[mode, mask] = 0
    engine.deleted.update(removed)
    return dict(correction_gram_nonzero_entries=correction_entries)
