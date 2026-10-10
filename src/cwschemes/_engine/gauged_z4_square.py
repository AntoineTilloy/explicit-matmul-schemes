"""Gauged Strassen embedding and certified square-only column deletion.

The underlying odd-character unit remains unchanged. Deleted indices refer to
its original term order and are meaningful only after this outer embedding.
This loader checks wiring/count metadata; the cited deletion proof establishes
that each removed square term has a zero factor over Q.
"""
from bisect import bisect_right
from fractions import Fraction as F
import hashlib
import json

from .. import _resources
from .square_scheme import base_factors
from .z4_quotient import unit_from_spec, bits

FORMAT = 'cw-z4-gauged-square-v1'

# Frozen direct-output call order used by the original research certificate.
RESEARCH_STRASSEN = (
    ((1,0,0,1),(1,0,0,1),(1,0,0,1)),
    ((0,0,1,1),(1,0,0,0),(0,0,1,-1)),
    ((1,0,0,0),(0,1,0,-1),(0,1,0,1)),
    ((0,0,0,1),(-1,0,1,0),(1,0,1,0)),
    ((1,1,0,0),(0,0,0,1),(-1,1,0,0)),
    ((-1,0,1,0),(1,1,0,0),(0,0,0,1)),
    ((0,1,0,-1),(0,0,1,1),(1,0,0,0)),
)


class GaugedStrassen:
    """Seven exact outer products, indexed by packet row orbit."""
    n = 2
    rank = 7
    recipe = dict(kind='lille', n=2, rank=7)

    def __init__(self, assignment, input_gauges):
        if (not isinstance(assignment, (list, tuple)) or len(assignment) != 7
                or any(type(x) is not int for x in assignment)
                or sorted(assignment) != list(range(7))):
            raise ValueError('The assignment must be a permutation of seven calls')
        if (not isinstance(input_gauges, (list, tuple)) or len(input_gauges) != 2
                or any(len(row) != 7 for row in input_gauges)):
            raise ValueError('Two input gauges of length seven are required')
        self.assignment = tuple(assignment)
        self.input_gauges = tuple(tuple(F(x) for x in row) for row in input_gauges)
        if any(not x for row in self.input_gauges for x in row):
            raise ValueError('Input gauges must be nonzero')
        self.output_gauge = tuple(1 / (a * b) for a, b in zip(*self.input_gauges))
        gauges = (*self.input_gauges, self.output_gauge)
        old = base_factors(2)
        if any(tuple(old[mode][coordinate][call] for coordinate in range(4))
               != RESEARCH_STRASSEN[call][mode] for call in range(7) for mode in range(3)):
            raise ValueError('Public and proof Strassen call conventions differ')
        self.factors = [[[F(old[mode][coordinate][call]) * gauges[mode][h]
                          for h, call in enumerate(self.assignment)]
                         for coordinate in range(4)] for mode in range(3)]
        self.check_tensor()

    def check_tensor(self):
        """Check every one of the 64 scalar MM2 coefficients over Q."""
        for a in range(4):
            i, k = divmod(a, 2)
            for b in range(4):
                kb, j = divmod(b, 2)
                for c in range(4):
                    got = sum((self.factors[0][a][h] * self.factors[1][b][h]
                               * self.factors[2][c][h] for h in range(7)), F(0))
                    if got != int(k == kb and c == 2 * i + j):
                        raise ValueError('The gauged outer tensor is not MM2')
        return dict(coefficient_field='Q', all_coefficients=64, exact=True)

    def factor_entry(self, mode, term, row, col):
        if mode not in (0, 1, 2) or not (0 <= term < 7 and 0 <= row < 2 and 0 <= col < 2):
            raise IndexError((mode, term, row, col))
        return self.factors[mode][2 * row + col][term]


class GaugedZ4Square:
    def __init__(self, data):
        if data.get('format') != FORMAT or data.get('coefficient_field', 'Q') != 'Q':
            raise ValueError('Unsupported gauged square certificate')
        self.unit = unit_from_spec(data['unit'])
        if self.unit.input_character != 'odd' or len(self.unit.rows) != 7:
            raise ValueError('The gauged square requires the original odd seven-product unit')
        self.outer = GaugedStrassen(data['assignment'], data['input_gauges'])
        self.packets = [self.unit]
        self.inner = self.unit.q ** self.unit.m
        self.full_n = 2 * self.inner
        self.n = data['n']
        if type(self.n) is not int or not 1 <= self.n <= self.full_n or data['full_n'] != self.full_n:
            raise ValueError('Square dimension mismatch')
        self.full_rank = self.unit.rank
        deleted = data['deleted_terms']
        if (not isinstance(deleted, list) or any(type(t) is not int for t in deleted)
                or deleted != sorted(set(deleted))
                or any(not 0 <= t < self.full_rank for t in deleted)):
            raise ValueError('Deleted terms must be sorted unique original unit indices')
        self.deleted_terms = tuple(deleted)
        # d_i-i is the insertion point of original deleted term d_i in the
        # retained sequence. Repeated thresholds handle consecutive deletions.
        self._skip_thresholds = tuple(t - i for i, t in enumerate(deleted))
        self.rank = self.full_rank - len(deleted)
        if self.rank != data['rank'] or data.get('full_rank', self.full_rank) != self.full_rank:
            raise ValueError('Gauged square rank mismatch')
        self.proof_scope = data.get('proof_scope', '')
        if deleted and not self.proof_scope:
            raise ValueError('Square-only deletions require an explicit proof scope')
        self.deletion_proof = data.get('deletion_proof')
        self._deletion_audit = None
        if isinstance(self.deletion_proof, dict):
            # Bind any inline proof metadata to the exact square wiring. A
            # standalone arithmetic verifier must replay its zero-factor proof.
            for field, expected in self.deletion_binding().items():
                if field in self.deletion_proof and self.deletion_proof[field] != expected:
                    raise ValueError('Deletion proof does not match the square: ' + field)

    def deletion_binding(self):
        return dict(coefficient_field='Q', n=self.n, full_n=self.full_n,
                    full_rank=self.full_rank, assignment=list(self.outer.assignment),
                    input_gauges=[[str(x) for x in row] for row in self.outer.input_gauges],
                    deleted_terms=list(self.deleted_terms))

    def unit_ordinal(self, term):
        if type(term) is not int or not 0 <= term < self.rank:
            raise IndexError(term)
        return term + bisect_right(self._skip_thresholds, term)

    def full_factor_entry(self, mode, original_term, row, col):
        """Coefficient before deletion, in the unchanged unit's term order."""
        if (mode not in (0, 1, 2) or not 0 <= original_term < self.full_rank
                or not (0 <= row < self.n and 0 <= col < self.n)):
            raise IndexError((mode, original_term, row, col))
        br, ir = divmod(row, self.inner)
        bc, ic = divmod(col, self.inner)
        value = F(0)
        for h in range(7):
            weight = self.outer.factor_entry(mode, h, br, bc)
            if weight:
                variable = self.unit.matrix_variable(h, mode, ir, ic)
                value += weight * self.unit.factor_entry(mode, original_term, variable)
        return value

    def factor_entry(self, mode, term, row, col):
        return self.full_factor_entry(mode, self.unit_ordinal(term), row, col)

    def audit(self):
        return dict(exact_outer_tensor=self.outer.check_tensor(),
                    unit_rank_unchanged=str(self.unit.rank),
                    square_deleted_terms=len(self.deleted_terms),
                    exact_index_map=True, proof_scope=self.proof_scope,
                    deletion_proof_replayed=False)

    def _read_deletion_proof(self):
        ref = self.deletion_proof
        if not isinstance(ref, dict) or set(ref) != {'file', 'sha256'}:
            raise ValueError('A hash-bound deletion proof resource is required')
        raw = _resources.materialize(ref['file']).read_bytes()
        if hashlib.sha256(raw).hexdigest() != ref['sha256']:
            raise ValueError('Deletion proof hash mismatch')
        return json.loads(raw)

    @staticmethod
    def _normalized_labels(labels):
        names = {}; normalized = []
        for label in labels:
            if label:
                names.setdefault(label, len(names) + 1)
                normalized.append(names[label])
            else:
                normalized.append(0)
        return tuple(normalized)

    def _zero_factor_polynomial(self, mode, H, labels):
        """Rebuild the full physical factor in the centered label basis.

        A key records one of four outer blocks and all six physical digit
        slots. Coefficient cancellation in this formal basis proves zero at
        every matrix entry. Label normalization is a simultaneous S_q
        relabeling, not independent relabeling of different coordinates.
        """
        sup, pk = self.unit.sup, self.unit.packet
        polynomial = {}
        if mode in (0, 1):
            selected = [(rho, H, labels, self.unit.row_character(rho))
                        for rho in sup.rowlist]
        elif mode == 2:
            selected = [(rho, H2, lab2, 1) for rho in sup.reps
                        for H2, lab2 in self.unit.orbits.members(H, labels)]
        else:
            raise ValueError('Invalid zero-factor mode')
        for rho, pattern, lab, character in selected:
            h, g = sup.rowinfo[rho]
            a, b, c = sup.reps[h]
            first, second = ((bits(b), bits(c)), (bits(c), bits(a)), (bits(b), bits(a)))[mode]
            positions = {g[t]:j for j,t in enumerate(first + second)}
            table = pk.table(pk.perm.index(mode), rho[mode]).get(pattern, {})
            for I, coefficient in table.items():
                monomial = [0] * 6
                for t in bits(I):
                    monomial[positions[t]] = lab[t]
                for block in range(4):
                    weight = self.outer.factors[mode][block][h]
                    if weight:
                        key = (block, tuple(monomial))
                        polynomial[key] = polynomial.get(key, F(0)) + character * weight * coefficient
        # A formal zero is valid even when the normalized q-th centered
        # function is dependent (notably the two-label q=2 replay).
        return {key:value for key,value in polynomial.items() if value}

    def audit_deletions(self):
        """Exactly replay every square-only zero factor from public tables."""
        if self._deletion_audit is not None:
            return dict(self._deletion_audit)
        if not self.deleted_terms:
            result = dict(exact=True, coefficient_field='Q', deleted_terms=0,
                          distinct_normalized_factors=0, unit_rank_unchanged=str(self.unit.rank))
            self._deletion_audit = result
            return dict(result)
        proof = self._read_deletion_proof()
        if (proof.get('coefficient_field') != 'Q'
                or proof.get('assignment') != list(self.outer.assignment)
                or tuple(tuple(F(x) for x in row) for row in proof.get('scales', []))
                != (*self.outer.input_gauges, self.outer.output_gauge)):
            raise ValueError('Deletion proof does not match the exact outer gauges')
        gates = [gate for gate in proof.get('q_gates', []) if gate.get('q') == self.unit.q]
        if len(gates) != 1:
            raise ValueError('Deletion proof needs exactly one gate at this q')
        gate = gates[0]
        if (gate.get('n') != self.full_n or gate.get('original_rank') != self.full_rank
                or gate.get('deleted_terms') != list(self.deleted_terms)
                or gate.get('saving') != len(self.deleted_terms) or gate.get('rank') != self.rank):
            raise ValueError('Deletion proof dimension, ordinal or rank mismatch')
        records = gate.get('numeric_label_column_records', [])
        if (len(records) != len(self.deleted_terms)
                or sorted(r['old_term_index'] for r in records) != list(self.deleted_terms)):
            raise ValueError('Deletion proof must bind every original term once')
        checked = set()
        for record in records:
            H, labels = self.unit.orbits.unrank(record['old_term_index'])
            if record.get('H') != H or tuple(record.get('labels', [])) != tuple(labels):
                raise ValueError('Deletion proof label/ordinal mismatch')
            modes = record.get('exact_zero_modes', [])
            if (not modes or len(set(modes)) != len(modes)
                    or any(type(mode) is not int or mode not in (0,1,2) for mode in modes)):
                raise ValueError('Deletion proof needs an exact zero-factor mode')
            normalized = self._normalized_labels(labels)
            for mode in modes:
                key = (mode, H, normalized)
                if key not in checked:
                    if self._zero_factor_polynomial(mode, H, normalized):
                        raise ValueError('Deleted square factor is nonzero over Q')
                    checked.add(key)
        result = dict(exact=True, coefficient_field='Q', deleted_terms=len(records),
                      distinct_normalized_factors=len(checked),
                      unit_rank_unchanged=str(self.unit.rank),
                      proof_sha256=self.deletion_proof['sha256'],
                      method='Exact centered-basis physical polynomials; simultaneous label relabeling')
        self._deletion_audit = result
        return dict(result)

    def certificate(self):
        data = dict(format=FORMAT, unit=self.unit.spec(), rank=self.rank,
                    proof_scope=self.proof_scope,
                    coefficient_generator='general_square.load_certificate',
                    **self.deletion_binding())
        if self.deletion_proof is not None:
            data['deletion_proof'] = self.deletion_proof
        return data


def load_certificate(data):
    return GaugedZ4Square(data)
