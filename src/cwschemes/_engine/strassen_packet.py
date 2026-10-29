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
blocks not used by a row then leaves exactly the direct sum of the rows. There are no
corrections and no pure-subset counting: the packet rank is (q+1)^N.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from .scheme import CompactPacket

LOCAL_TENSOR = 'strassen-completion-v1'


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


class StrassenPacket:
    """len(rows) independent <q^m> products in (q+1)^(3m) multiplications."""

    def __init__(self, q, rows, m):
        if type(q) is not int or q < 1:
            raise ValueError('q must be a positive integer')
        self.q = q
        self.m = m
        self.rows = [tuple(r) for r in rows]
        self.witness = audit(self.rows, m)
        self.rank = (q + 1) ** (3 * m)

    # The row/column wiring of the CW packets: variable 0 on a row's zero mask of
    # that leg, digits + 1 elsewhere.
    matrix_variable = CompactPacket.matrix_variable

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
        value = 1
        for i in reversed(range(3 * m)):
            index, choice = divmod(index, q + 1)
            value *= local_factor(q, i // m, mode, choice, variable[i])
            if not value:
                break
        return F(value)

    def spec(self):
        return dict(local_tensor=LOCAL_TENSOR, q=self.q, m=self.m, rows=[list(r) for r in self.rows], rank=self.rank)


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


def tensor_entry(packet, x, y, z, perturb=0):
    """Exact sum over all (q+1)^N terms of U(x) V(y) W(z) for the packet.

    The terms are all local choices, so the sum factors coordinate-wise into local sums
    of the generator's own local factors.
    """
    q, m = packet.q, packet.m
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
    assert packet.rank == (packet.q + 1) ** (3 * packet.m)
    identity = packet_identity(packet, samples=samples, seed=seed)
    assert identity['passed'], 'S_q packet identity fails'
    return dict(q=packet.q, m=packet.m, power=3 * packet.m, products=witness['products'],
                inner_side=str(packet.q ** packet.m), packet_rank=str(packet.rank),
                weight_classes=witness['weight_classes'], local_identity_coefficients=verify_local(packet.q),
                exact_packet_entries=identity)
