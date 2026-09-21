"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
import json
import math
from flint import fmpq_mat
from .fixed_sector_quotient import FixedSectorQuotient
from .quotient_invariant_gram import local_bases
ROOT = resource_root()

def input_path(path):
    path = Path(path)
    return path if path.exists() or path.is_absolute() else ROOT / path

class CompressedFixedQuotient:

    def __init__(self, proof, allow_stage=False):
        self.proof_path = input_path(proof)
        self.proof = json.loads(self.proof_path.read_text())
        assert self.proof['complete']
        assert self.proof['format'] == 'cw-quotient-original-orbit-rational-v1' or (allow_stage and self.proof['format'] == 'cw-quotient-successive-stage-rational-v1')
        self.metadata = json.loads(input_path(Path(self.proof['source']).with_suffix('.json')).read_text())
        d = self.metadata
        self.groups = list(map(tuple, d['groups']))
        self.q = d['q']
        self.m = d['m']
        permutation = list(range(3 * self.m))
        flips = []
        for group in self.groups:
            if len(group) == 1:
                flips.append(group[0])
            else:
                i, j = group
                permutation[i], permutation[j] = (j, i)
        self.base = FixedSectorQuotient(self.q, d['rows'], self.m, permutation, (-1, -1, 1), flips)
        self.pair = tuple(d['pair'])
        self.free_mode = 3 - sum(self.pair)
        self.domain_character = math.prod(((-1, -1, 1)[i] for i in self.pair))
        self.bases = [local_bases(self.q, len(g)) for g in self.groups]
        self.values = [[{tuple(x): int(v) for x, v in zip(b['labels'], b['vector'])} for b in bs] for bs in self.bases]
        self.orbit_lookup = {}
        for i, orbit in enumerate(d['orbits']):
            self.orbit_lookup[tuple(orbit['states'])] = self.orbit_lookup[tuple(orbit['image'])] = i
        self.deleted_orbits = set(self.proof['deleted_orbits'])
        self._digit_code = {tuple(ds): H for H, ds in self.base.spec['digits'].items()}
        self._block_lookup = {(H, None if data is None else data[1]): i for i, (H, data, _) in enumerate(self.base.blocks)}
        deleted = set()
        for orbit in self.deleted_orbits:
            states = d['orbits'][orbit]['states']
            choices = [next((b['labels'] for b in bs if b['state'] == state)) for bs, state in zip(self.bases, states)]
            for local in product(*choices):
                labels = [0] * (3 * self.m)
                for group, values in zip(self.groups, local):
                    for i, value in zip(group, values):
                        labels[i] = int(value)
                deleted.add(self.original_index(tuple(labels))[0])
        self.deleted = sorted(deleted)
        assert len(self.deleted) == self.proof['savings']
        assert all((self.term_orbit(self.base.term(t)[1]) in self.deleted_orbits for t in self.deleted))
        self.rank = self.base.rank - len(self.deleted)
        assert self.rank == self.proof['rank']
        self.transforms = {}
        for key, relation in self.proof['relations'].items():
            block = d['blocks'][int(key)]
            records = block['records']
            J = relation['deleted']
            B = relation['basis']
            metric = fmpq_mat([[self.metric(records[a], records[b]) for b in J] for a in J])
            C = fmpq_mat([[str(x) for x in row] for row in relation['coefficients']])
            W = C * metric.inv()
            for i, b in enumerate(B):
                for j, a in enumerate(J):
                    if W[i, j]:
                        r, s = (records[b], records[a])
                        weight = F(str(W[i, j])) * block['dimension']
                        self.transforms.setdefault((r['orbit'], s['orbit']), []).append((r, s, weight))
        for name in ('n', 'rows', 'perm', 'characters', 'label_flips', 'label_fixed', 'products', 'free', 'fixed', 'fixed_terms'):
            setattr(self, name, getattr(self.base, name))
        self.updated_orbits = {a for a, b in self.transforms}

    def updates(self, original):
        return self.term_orbit(self.base.term(original)[1]) in self.updated_orbits

    def term_orbit(self, labels):
        states = []
        for group in self.groups:
            local = [labels[i] for i in group]
            if len(group) == 1:
                states.append(2 if self.q % 2 and local[0] == self.q else int(local[0] != 0))
            elif local == [0, 0]:
                states.append(0)
            elif local[1] == 0:
                states.append(1)
            elif local[0] == 0:
                states.append(2)
            elif local[0] == local[1]:
                states.append(3)
            else:
                states.append(4)
        return self.orbit_lookup[tuple(states)]

    def original_index(self, labels):
        """Invert the original unranking; also return the orientation sign."""
        labels = tuple(map(int, labels))
        ds = tuple((1 if x else 2 for x in labels))
        H = self._digit_code[ds]
        gH = self._digit_code[self.base.permute(ds)]
        sign = 1
        if H > gH:
            H = gH
            labels = self.base.variable_image(labels)
            sign = -1
        if H < max(H, gH) or (H, None) in self._block_lookup:
            b = self._block_lookup[H, None]
            _, _, ordinary = self.base.blocks[b]
            local = 0
            for i in ordinary:
                local = local * self.q + labels[i] - 1
        else:
            candidates = [(i, data) for i, (h, data, _) in enumerate(self.base.blocks) if h == H]
            components = candidates[0][1][0]
            pivot = None
            reverse = False
            for j, (kind, i, k, _, _, _) in enumerate(components):
                if kind == 'pair' and labels[i] != labels[k]:
                    pivot = j
                    reverse = labels[i] > labels[k]
                    break
                if kind == 'flip' and self.base.label_image(labels[i]) != labels[i]:
                    pivot = j
                    reverse = labels[i] > self.base.label_image(labels[i])
                    break
            assert pivot is not None, 'A fixed term has no negative-character quotient column.'
            if reverse:
                labels = self.base.variable_image(labels)
                sign = -sign
            b = self._block_lookup[H, pivot]
            _, data, radices = self.base.blocks[b]
            local = 0
            for j, ((kind, i, k, _, _, _), radix) in enumerate(zip(components, radices)):
                x, y = (labels[i] - 1, labels[k] - 1)
                if kind == 'pair':
                    value = x if j < pivot else x * (2 * self.q - x - 1) // 2 + y - x - 1 if j == pivot else x * self.q + y
                elif kind == 'flip':
                    value = labels[i] - (self.q - self.base.label_fixed) - 1 if j < pivot else x // 2 if j == pivot else x
                else:
                    value = x
                local = local * radix + value
        index = (self.base.ends[b - 1] if b else 0) + local
        assert self.base.term(index)[1] == labels
        return (index, sign)

    def metric(self, left, right):

        def dot(indices):
            value = 1
            for bs, i, j in zip(self.bases, left['indices'], indices):
                if bs[i]['state'] != bs[j]['state']:
                    return 0
                value *= int(bs[i]['vector'] @ bs[j]['vector'])
            return value
        return left['factor'] * right['factor'] * (dot(right['indices']) + self.domain_character * right['sign'] * dot(right['image']))

    @lru_cache(100000)
    def local_moment(self, group, i, j, left, right):
        bs = self.bases[group]
        a, b = (bs[i], bs[j])
        if left not in self.values[group][i] or right not in self.values[group][j]:
            return F(0)
        assert a['irrep'] == b['irrep']
        if len(self.groups[group]) == 1:
            irrep = a['irrep']
            if irrep == 0:
                return F(1)
            paired = self.q - self.q % 2
            x, y = (left[0], right[0])
            same_pair = (x - 1) // 2 == (y - 1) // 2
            if irrep == 1:
                return F(4, paired) if same_pair else F(-8, paired * (paired - 2))
            return F(2 if x == y else -2, paired) if same_pair else F(0)
        labels = sorted(set(left + right) - {0})
        canonical = {x: k + 1 for k, x in enumerate(labels)}
        cl = tuple((canonical.get(x, 0) for x in left))
        cr = tuple((canonical.get(x, 0) for x in right))
        return self.pair_moment(group, i, j, cl, cr, len(labels))

    @lru_cache(100000)
    def pair_moment(self, group, i, j, left, right, distinct):
        total = 0
        count = 0
        av, bv = (self.values[group][i], self.values[group][j])
        for assignment in permutations(range(1, self.q + 1), distinct):
            permutation = (0,) + assignment
            total += av[tuple((permutation[x] for x in left))] * bv[tuple((permutation[x] for x in right))]
            count += 1
        return F(total, count)

    def averaged_vectors(self, left_record, right_record, left_labels, right_labels):
        result = F(0)
        for a in range(2):
            sl = self.base.variable_image(left_labels) if a else left_labels
            for b in range(2):
                tl = self.base.variable_image(right_labels) if b else right_labels
                value = F(self.domain_character ** (a + b))
                for k, (group, i, j) in enumerate(zip(self.groups, left_record['indices'], right_record['indices'])):
                    value *= self.local_moment(k, i, j, tuple((sl[x] for x in group)), tuple((tl[x] for x in group)))
                    if not value:
                        break
                result += value
        return result * left_record['factor'] * right_record['factor']

    @lru_cache(20000)
    def relation_entry(self, kept_original, deleted_original):
        _, s, _ = self.base.term(kept_original)
        _, t, _ = self.base.term(deleted_original)
        so, to = (self.term_orbit(s), self.term_orbit(t))
        assert so not in self.deleted_orbits and to in self.deleted_orbits
        value = sum((weight * self.averaged_vectors(a, b, s, t) for a, b, weight in self.transforms.get((so, to), ())), F(0))
        sign = (-1) ** (sum((x == 0 for x in s)) + sum((x == 0 for x in t))) if 2 in self.pair and self.proof.get('completion_sign_omitted', True) else 1
        return sign * value

    def retained_original(self, t):
        if not 0 <= t < self.rank:
            raise IndexError(t)
        index = t
        while True:
            new = t + bisect_right(self.deleted, index)
            if new == index:
                return index
            index = new

    def factor_entry(self, mode, t, product, row, col):
        original = self.retained_original(t)
        value = self.base.factor_entry(mode, original, product, row, col)
        if mode == self.free_mode and self.updates(original):
            for deleted in self.deleted:
                coefficient = self.base.factor_entry(mode, deleted, product, row, col)
                if coefficient:
                    value += self.relation_entry(original, deleted) * coefficient
        return value

class SuccessiveCompressedFixedQuotient:

    def __init__(self, proof):
        self.proof_path = input_path(proof)
        self.proof = json.loads(self.proof_path.read_text())
        assert self.proof['format'] == 'cw-quotient-successive-rational-v1' and self.proof['complete']
        self.steps = [CompressedFixedQuotient(path, allow_stage=True) for path in self.proof['proofs']]
        self.base = self.steps[0].base
        self.deleted = sorted((t for step in self.steps for t in step.deleted))
        assert len(self.deleted) == len(set(self.deleted)) == self.proof['saved_multiplications']
        self.rank = self.base.rank - len(self.deleted)
        assert self.rank == self.proof['rank']
        for name in ('q', 'm', 'n', 'rows', 'perm', 'characters', 'label_flips', 'label_fixed', 'products', 'free', 'fixed', 'fixed_terms'):
            setattr(self, name, getattr(self.base, name))
    retained_original = CompressedFixedQuotient.retained_original

    @lru_cache(20000)
    def original_factor(self, stage, mode, original, product, row, col):
        if stage < 0:
            return self.base.factor_entry(mode, original, product, row, col)
        value = self.original_factor(stage - 1, mode, original, product, row, col)
        step = self.steps[stage]
        if mode == step.free_mode and step.updates(original):
            for deleted in step.deleted:
                coefficient = self.original_factor(stage - 1, mode, deleted, product, row, col)
                if coefficient:
                    value += step.relation_entry(original, deleted) * coefficient
        return value

    def factor_entry(self, mode, t, product, row, col):
        return self.original_factor(len(self.steps) - 1, mode, self.retained_original(t), product, row, col)

def compressed_packet(proof):
    kind = json.loads(input_path(proof).read_text())['format']
    return SuccessiveCompressedFixedQuotient(proof) if kind == 'cw-quotient-successive-rational-v1' else CompressedFixedQuotient(proof)
from .fixed_quotient_square import FixedQuotientSquare

class CompressedQuotientSquare(FixedQuotientSquare):

    def __init__(self, original_certificate, proof, ordinary_recipes=None):
        d = json.loads(input_path(original_certificate).read_text()) if isinstance(original_certificate, (str, Path)) else original_certificate
        packet = compressed_packet(proof)
        assert packet.q == d['q'] and packet.base.rank == d['packet_rank']
        super().__init__(packet, d['a'], d['b'], d['leaf_recipe'], d['n'], d.get('ordinary_recipes') if ordinary_recipes is None else ordinary_recipes)
        self.original_certificate = d
        self.compression_proof = str(proof)
        self.replacement_ordinary_recipes = ordinary_recipes

    def certificate(self):
        result = dict(format='cw-compressed-fixed-quotient-square-v1', coefficient_field='Q', n=self.n, rank=self.rank, original_square=self.original_certificate, quotient_compression_proof=self.compression_proof, packet_rank=self.packet.rank, saved_multiplications=self.original_certificate['rank'] - self.rank, coefficient_generator='compressed_fixed_quotient.load_square_certificate')
        if self.replacement_ordinary_recipes is not None:
            result['replacement_ordinary_recipes'] = self.replacement_ordinary_recipes
        return result

def load_square_certificate(value):
    d = json.loads(input_path(value).read_text()) if isinstance(value, (str, Path)) else value
    assert d['format'] == 'cw-compressed-fixed-quotient-square-v1'
    square = CompressedQuotientSquare(d['original_square'], d['quotient_compression_proof'], d.get('replacement_ordinary_recipes'))
    assert square.rank == d['rank'] and square.packet.rank == d['packet_rank']
    return square
