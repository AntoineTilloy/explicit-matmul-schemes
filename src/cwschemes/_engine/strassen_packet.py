"""Exact coefficient evaluation for Strassen-tensor (S_q) packets. No discovery entry points.

Strassen's tensor Str_q = sum_{t=1..q} (x_0 y_t z_t + x_t y_0 z_t) has the rank-(q+1)
integer completion

    S(a, b, c) = sum_{t=1..q} (a_0 + a_t)(b_0 + b_t) c_t - a_0 b_0 sum_{t=1..q} c_t
               = Str_q(a, b, c) + sum_{t=1..q} a_t b_t c_t,

whose only error lies on the all-positive block. The c argument has no index 0.

A power-N packet, N = 3m, uses kind i // m on coordinate i: kind 0 is S(x, y, z), kind 1 is
S(y, z, x) and kind 2 is S(z, x, y). Rows are zero masks (A, B, C) partitioning the N
coordinates, m coordinates each. A at coordinate i means that the X leg is the scalar block
there (index 0); likewise B for Y and C for Z. Kind 0 allows only A or B, kind 1 only B or C
and kind 2 only C or A. Row matrix indices are x on B u C, y on C u A and z on A u B, with the
same wiring as CW packets, so every row is a <q^m, q^m, q^m> product.

Every local support block has at most one zero index. A block triple (X_u, Y_v, Z_w) has
|A_u| + |B_v| + |C_w| = N zero indices, so it survives only if A_u, B_v and C_w partition
the coordinates. audit() checks that no such triple exists unless u = v = w. Zeroing all
blocks not used by a row then leaves exactly the direct sum of the rows.

With compression='zero-trimmed-v1', local choice vectors whose X, Y or Z factor vanishes
identically across all packet rows are omitted: a special-choice mask triple (q_0, q_1, q_2)
on kinds 0, 1, 2 has nonzero X, Y, Z factors over Q if and only if (q_0, q_2) <= (A_u & K_0,
A_u >> 2m), (q_0, q_1) <= (B_v & K_0, B_v >> m & K_0), and (q_1, q_2) <= (C_w >> m & K_0,
C_w >> 2m) for some rows u, v, w.
"""
from bisect import bisect_right
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import accumulate, product
from .scheme import CompactPacket

LOCAL_TENSOR = 'strassen-completion-v1'
COMPRESSION = 'zero-trimmed-v1'


def role(kind, mode):
    """Argument of S read by leg `mode` (0 = x, 1 = y, 2 = z): 0 = a, 1 = b, 2 = c."""
    return (mode - kind) % 3


def local_factor(q, kind, mode, choice, x):
    """Integer coefficient of index x on leg `mode` of local term `choice`.

    Choices 0..q-1 are (e_0 + e_t)(e_0 + e_t) e_t with t = choice + 1. Choice q is
    e_0 e_0 (e_1 + ... + e_q) with sign -1, which is absorbed in mode Z.
    """
    if not 0 <= x <= q:
        raise ValueError('invalid local index')
    if role(kind, mode) == 2:
        if x == 0:
            raise ValueError('scalar index on the vector leg of S')
        value = int(choice == q or x == choice + 1)
    else:
        value = int(x == 0 or (choice < q and x == choice + 1))
    return -value if mode == 2 and choice == q else value


def verify_local(q):
    """Exact local identity of all three kinds; return the number of checked coefficients.

    One zero index: the delta on the other two (the Strassen blocks). Two zeros: 0.
    No zero: the error, 1 exactly on the diagonal. Hence no block has two zeros.
    """
    if type(q) is not int or q < 1:
        raise ValueError('q must be a positive integer')
    checked = 0
    for kind in range(3):
        ranges = [range(1, q + 1) if role(kind, mode) == 2 else range(q + 1) for mode in range(3)]
        for v in product(*ranges):
            value = sum(local_factor(q, kind, 0, s, v[0]) * local_factor(q, kind, 1, s, v[1])
                        * local_factor(q, kind, 2, s, v[2]) for s in range(q + 1))
            zeros = [mode for mode in range(3) if v[mode] == 0]
            if not zeros:
                expected = int(v[0] == v[1] == v[2])
            elif len(zeros) == 1:
                a, b = (v[mode] for mode in range(3) if mode != zeros[0])
                expected = int(a == b)
            else:
                expected = 0
            assert value == expected, (q, kind, v, value, expected)
            checked += 1
    return checked


def audit(rows, m):
    """Check the S_q support condition and return the weight-class profile."""
    rows = [tuple(r) for r in rows]
    k = 3 * m
    full = (1 << k) - 1
    kinds = [((1 << m) - 1) << (m * j) for j in range(3)]
    assert type(m) is int and m >= 1 and rows and len(set(rows)) == len(rows)
    for r in rows:
        assert len(r) == 3 and all(type(a) is int and 0 <= a <= full and a.bit_count() == m for a in r)
        a, b, c = r
        assert not (a & b or a & c or b & c) and a | b | c == full
        assert not (a & kinds[1] or b & kinds[2] or c & kinds[0]), 'letter not allowed by the coordinate kind'
    assert all(len({r[j] for r in rows}) == len(rows) for j in range(3)), 'shared input/output block'
    # No C occurs on kind 0, so a partition A_u, B_v, C_w forces B_v & kind0 = kind0 \ A_u.
    # Grouping by that part enumerates every partition triple.
    by_kind0 = {}
    for v, r in enumerate(rows):
        by_kind0.setdefault(r[1] & kinds[0], []).append(v)
    output = {r[2]: w for w, r in enumerate(rows)}
    triples = 0
    for u, r in enumerate(rows):
        for v in by_kind0.get(kinds[0] & ~r[0], ()):
            if r[0] & rows[v][1]:
                continue
            w = output.get(full ^ (r[0] | rows[v][1]))
            if w is not None:
                assert u == v == w, f'cross triple {(u, v, w)}'
                triples += 1
    assert triples == len(rows)
    weights = Counter((r[0] & kinds[0]).bit_count() for r in rows)
    return dict(m=m, k=k, products=len(rows), partition_triples=triples,
                weight_classes=[weights.get(t, 0) for t in range(m + 1)])


def _build_downsets(rows, m):
    """Build exact 2D downset tables DA[q0, q2], DB[q0, q1], DC[q1, q2] of shape (2^m, 2^m) uint8."""
    import numpy as np
    K0 = (1 << m) - 1
    size2 = 1 << (2 * m)
    D_A = np.zeros(size2, dtype=np.uint8)
    D_B = np.zeros(size2, dtype=np.uint8)
    D_C = np.zeros(size2, dtype=np.uint8)
    for a, b, c in rows:
        D_A[(a & K0) | (((a >> (2 * m)) & K0) << m)] = 1
        D_B[(b & K0) | (((b >> m) & K0) << m)] = 1
        D_C[((c >> m) & K0) | (((c >> (2 * m)) & K0) << m)] = 1
    for arr in (D_A, D_B, D_C):
        view = arr.reshape([2] * (2 * m))
        for bit in range(2 * m):
            s0 = [slice(None)] * (2 * m)
            s1 = [slice(None)] * (2 * m)
            s0[bit] = 0
            s1[bit] = 1
            view[tuple(s0)] |= view[tuple(s1)]
    M0 = 1 << m
    return D_A.reshape(M0, M0).T, D_B.reshape(M0, M0).T, D_C.reshape(M0, M0).T


def _compute_trimmed_index(DA, DB, DC, q, m):
    """Exact q2_sum[q0, q1] and q0_ends over Q/Z for zero-trimmed S_q unranking."""
    import numpy as np
    M0 = 1 << m
    pop = np.array([i.bit_count() for i in range(M0)], dtype=np.int64)
    if (q + 1) ** m < (1 << 53):
        wq = np.array([float(q ** (m - int(pop[i]))) for i in range(M0)], dtype=np.float64)
        q2_sum = np.rint((DA.astype(np.float64) * wq[None, :]) @ DC.astype(np.float64).T).astype(np.int64)
        q2_sum *= DB.astype(np.int64)
        q_pows = [q ** d for d in range(2 * m + 1)]
        q0_counts = []
        for q0 in range(M0):
            row = q2_sum[q0]
            nz = np.flatnonzero(row)
            if not len(nz):
                q0_counts.append(0)
                continue
            w0 = int(pop[q0])
            tot = 0
            for q1 in nz:
                tot += q_pows[2 * m - w0 - int(pop[q1])] * int(row[q1])
            q0_counts.append(tot)
        return q2_sum, list(accumulate(q0_counts)), pop
    # Arbitrary-precision fallback when (q+1)^m >= 2^53: group q2 by popcount.
    q2_sum_obj = np.zeros((M0, M0), dtype=object)
    for w2 in range(m + 1):
        sel = pop == w2
        if not np.any(sel):
            continue
        joint = (DA[:, sel].astype(np.int32) @ DC[:, sel].astype(np.int32).T) * DB
        factor = q ** (m - w2)
        q2_sum_obj += joint.astype(object) * factor
    q0_counts = []
    for q0 in range(M0):
        w0 = int(pop[q0])
        tot = sum(int(q2_sum_obj[q0, q1]) * (q ** (2 * m - w0 - int(pop[q1])))
                  for q1 in range(M0) if q2_sum_obj[q0, q1])
        q0_counts.append(tot)
    return q2_sum_obj, list(accumulate(q0_counts)), pop


class StrassenPacket:
    """len(rows) independent <q^m> products in (q+1)^(3m) (or zero-trimmed) multiplications."""
    compression = None

    def __init__(self, q, rows, m, compression=None):
        if type(q) is not int or q < 1:
            raise ValueError('q must be a positive integer')
        if compression not in (None, COMPRESSION):
            raise ValueError('unknown S_q packet compression')
        self.q = q
        self.m = m
        self.compression = compression
        self.rows = [tuple(r) for r in rows]
        self.witness = audit(self.rows, m)
        if compression is None:
            self.rank = (q + 1) ** (3 * m)
        else:
            self._DA, self._DB, self._DC = _build_downsets(self.rows, m)
            self._q2_sum, self._q0_ends, self._pop = _compute_trimmed_index(
                self._DA, self._DB, self._DC, q, m
            )
            self.rank = self._q0_ends[-1]

    # The row/column wiring of the CW packets: variable 0 on a row's zero mask of
    # that leg, digits + 1 elsewhere.
    matrix_variable = CompactPacket.matrix_variable

    @lru_cache(4096)
    def _unrank_trimmed(self, index):
        q, m = self.q, self.m
        M0 = 1 << m
        q0 = bisect_right(self._q0_ends, index)
        rem = index - (self._q0_ends[q0 - 1] if q0 else 0)
        w0 = int(self._pop[q0])
        row = self._q2_sum[q0]
        q1_chosen = None
        for q1 in range(M0):
            s = int(row[q1])
            if not s:
                continue
            cnt = (q ** (2 * m - w0 - int(self._pop[q1]))) * s
            if rem < cnt:
                q1_chosen = q1
                break
            rem -= cnt
        assert q1_chosen is not None
        ord01_ways = q ** (2 * m - w0 - int(self._pop[q1_chosen]))
        q2_rem, ord01_idx = divmod(rem, ord01_ways)
        da_row = self._DA[q0]
        dc_row = self._DC[q1_chosen]
        q2_chosen = None
        for q2 in range(M0):
            if da_row[q2] and dc_row[q2]:
                cnt2 = q ** (m - int(self._pop[q2]))
                if q2_rem < cnt2:
                    q2_chosen = q2
                    break
                q2_rem -= cnt2
        assert q2_chosen is not None
        ord2_idx = q2_rem
        special = q0 | (q1_chosen << m) | (q2_chosen << (2 * m))
        choices = [0] * (3 * m)
        for i in reversed(range(2 * m, 3 * m)):
            if (special >> i) & 1:
                choices[i] = q
            else:
                ord2_idx, choices[i] = divmod(ord2_idx, q)
        assert ord2_idx == 0
        for i in reversed(range(2 * m)):
            if (special >> i) & 1:
                choices[i] = q
            else:
                ord01_idx, choices[i] = divmod(ord01_idx, q)
        assert ord01_idx == 0
        return tuple(choices)

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2):
            raise ValueError('mode must be 0,1,2')
        if not 0 <= index < self.rank:
            raise IndexError(index)
        q, m = self.q, self.m
        if len(variable) != 3 * m or any(not 0 <= x <= q for x in variable):
            raise ValueError('invalid variable')
        if any(x == 0 and role(i // m, mode) == 2 for i, x in enumerate(variable)):
            raise ValueError('scalar index on the vector leg of S')
        if self.compression == COMPRESSION:
            choices = self._unrank_trimmed(index)
            value = 1
            for i in range(3 * m):
                value *= local_factor(q, i // m, mode, choices[i], variable[i])
                if not value:
                    break
            return F(value)
        value = 1
        for i in reversed(range(3 * m)):
            index, choice = divmod(index, q + 1)
            value *= local_factor(q, i // m, mode, choice, variable[i])
            if not value:
                break
        return F(value)

    def spec(self):
        out = dict(local_tensor=LOCAL_TENSOR, q=self.q, m=self.m, rows=[list(r) for r in self.rows], rank=self.rank)
        if self.compression is not None:
            out['compression'] = self.compression
        return out


def packets(engine):
    """All S_q packets reachable from a coefficient engine (MixedSquare packets and nested 'cw' recipes)."""
    from .general_square import MixedSquare, Ordinary
    out = []
    def walk(x):
        if isinstance(x, MixedSquare):
            out.extend(p for p in x.packets if isinstance(p, StrassenPacket))
            walk(x.outer)
            walk(x.leaf)
        elif isinstance(x, Ordinary):
            for name in ('child', 'left', 'right'):
                if hasattr(x, name):
                    walk(getattr(x, name))
    walk(engine)
    return out


def _sample(packet, kind, rng):
    """Draw packet variables (x, y, z) of one sample kind and their target value."""
    s = packet.q ** packet.m
    M = len(packet.rows)
    u = rng.randrange(M)
    i, k, j = (rng.randrange(s) for _ in range(3))
    if kind == 'valid':
        return (packet.matrix_variable(u, 0, i, k), packet.matrix_variable(u, 1, k, j),
                packet.matrix_variable(u, 2, i, j)), 1
    if kind == 'label':
        # Same blocks, one contracted index changed.
        k2 = (k + 1 + rng.randrange(s - 1)) % s if s > 1 else k
        return (packet.matrix_variable(u, 0, i, k), packet.matrix_variable(u, 1, k2, j),
                packet.matrix_variable(u, 2, i, j)), int(k2 == k)
    v, w = rng.randrange(M), rng.randrange(M)
    return (packet.matrix_variable(u, 0, i, k), packet.matrix_variable(v, 1, k, j),
            packet.matrix_variable(w, 2, i, j)), int(u == v == w)


def _subsets_of_mask(mask):
    out = []
    sub = mask
    while True:
        out.append(sub)
        if sub == 0:
            break
        sub = (sub - 1) & mask
    return out


def tensor_entry(packet, x, y, z, perturb=0):
    """Exact sum over all terms of U(x) V(y) W(z) for the packet over Q."""
    q, m = packet.q, packet.m
    if packet.compression == COMPRESSION:
        # Sum explicitly over surviving special-choice masks (q0, q1, q2) in DA & DB & DC!
        ord_sums = []
        spec_prods = []
        active_masks = [0, 0, 0]
        for i in range(3 * m):
            kind = i // m
            o_val = sum(local_factor(q, kind, 0, s, x[i]) * local_factor(q, kind, 1, s, y[i])
                        * local_factor(q, kind, 2, s, z[i]) for s in range(q))
            s_val = (local_factor(q, kind, 0, q, x[i]) * local_factor(q, kind, 1, q, y[i])
                     * local_factor(q, kind, 2, q, z[i]))
            ord_sums.append(o_val)
            spec_prods.append(s_val)
            if s_val:
                active_masks[kind] |= 1 << (i - kind * m)
        # Coordinates where BOTH o_val == 0 and s_val == 0 immediately kill the entire sum.
        if any(ord_sums[i] == 0 and spec_prods[i] == 0 for i in range(3 * m)):
            return F(perturb)
        # Precompute product over each kind k for each subset q_k of active_masks[k]:
        kind_vals = []
        for kind in range(3):
            d = {}
            for qk in _subsets_of_mask(active_masks[kind]):
                prod_k = 1
                base = kind * m
                for bit in range(m):
                    prod_k *= spec_prods[base + bit] if ((qk >> bit) & 1) else ord_sums[base + bit]
                    if not prod_k:
                        break
                if prod_k:
                    d[qk] = prod_k
            kind_vals.append(d)
        DA, DB, DC = packet._DA, packet._DB, packet._DC
        total = 0
        for q0, v0 in kind_vals[0].items():
            for q1, v1 in kind_vals[1].items():
                if not DB[q0, q1]:
                    continue
                v01 = v0 * v1
                for q2, v2 in kind_vals[2].items():
                    if DA[q0, q2] and DC[q1, q2]:
                        total += v01 * v2
        return F(total + perturb)

    value = 1
    for i in range(3 * m):
        kind = i // m
        value *= sum(local_factor(q, kind, 0, s, x[i]) * local_factor(q, kind, 1, s, y[i])
                     * local_factor(q, kind, 2, s, z[i]) for s in range(q + 1))
        if not value:
            break
    return F(value + perturb)


def packet_identity(packet, samples=12, seed=1, perturb=0):
    """Exact packet identity at sampled entries: valid row entries, a changed label, random block triples."""
    import random
    rng = random.Random(seed)
    kinds = ('valid', 'label', 'cross')
    records = []
    for t in range(samples):
        kind = kinds[t % 3]
        (x, y, z), target = _sample(packet, kind, rng)
        got = tensor_entry(packet, x, y, z, perturb)
        records.append(dict(kind=kind, target=target, value=str(got), ok=got == target))
    return dict(samples=samples, targets_one=sum(r['target'] for r in records),
                mismatches=sum(not r['ok'] for r in records), passed=all(r['ok'] for r in records), exact=True)


def check(packet, samples=12, seed=1):
    """Local identity, count and sampled exact packet entries of one loaded (hence audited) packet."""
    witness = packet.witness
    assert witness['products'] == witness['partition_triples'] == len(packet.rows)
    if packet.compression is None:
        assert packet.rank == (packet.q + 1) ** (3 * packet.m)
    else:
        assert packet.compression == COMPRESSION
        DA, DB, DC = _build_downsets(packet.rows, packet.m)
        _, q0_ends, _ = _compute_trimmed_index(DA, DB, DC, packet.q, packet.m)
        assert packet.rank == q0_ends[-1] and packet.rank <= (packet.q + 1) ** (3 * packet.m)
    identity = packet_identity(packet, samples=samples, seed=seed)
    assert identity['passed'], 'S_q packet identity fails'
    out = dict(q=packet.q, m=packet.m, power=3 * packet.m, products=witness['products'],
               inner_side=str(packet.q ** packet.m), packet_rank=str(packet.rank),
               weight_classes=witness['weight_classes'], local_identity_coefficients=verify_local(packet.q),
               exact_packet_entries=identity)
    if packet.compression is not None:
        out['compression'] = packet.compression
    return out
