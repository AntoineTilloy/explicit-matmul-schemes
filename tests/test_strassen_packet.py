"""Strassen-tensor (S_q) packets: local identity, support audit and complete small tensors."""
from fractions import Fraction
from itertools import combinations, product
import numpy as np
import pytest

import cwschemes
from cwschemes._engine import strassen_packet as sp
from cwschemes._engine.general_square import MixedSquare, load_certificate

# Tripartite Ruzsa-Szemeredi sets: RS(1) = 1, RS(2) = 2 and an optimal RS(3) = 6.
RS = {1: [(0, 0, 0)], 2: [(0, 0, 0), (1, 1, 1)],
      3: [(0, 0, 0), (0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 2), (2, 2, 0)]}


LATIN = [(a, b, -(a + b) % 3) for a in range(3) for b in range(3)]


def support(m, three=None):
    """Rows from one RS set per weight class t, labels = t-subsets of range(m)."""
    K = [((1 << m) - 1) << (m * j) for j in range(3)]
    rows = []
    for t in range(m + 1):
        labels = [sum(1 << i for i in s) for s in combinations(range(m), t)]
        for a, b, c in (three if three and len(labels) == 3 else RS[len(labels)]):
            s0, s1, s2 = labels[a], labels[b] << m, labels[c] << 2 * m
            rows.append((s0 | (K[2] & ~s2), (K[0] & ~s0) | s1, (K[1] & ~s1) | s2))
    return rows


def unaudited(q, rows, m):
    packet = object.__new__(sp.StrassenPacket)
    packet.q, packet.m, packet.rows, packet.rank = q, m, [tuple(r) for r in rows], (q + 1) ** (3 * m)
    return packet


def packet_tensor(packet):
    """Complete packet tensor on the union of the row blocks, by brute force over all terms."""
    s = packet.q ** packet.m
    keys = [[(h, i, j) for h in range(len(packet.rows)) for i in range(s) for j in range(s)] for _ in range(3)]
    factors = [np.array([[int(packet.factor_entry(mode, t, packet.matrix_variable(h, mode, i, j)))
                          for h, i, j in keys[mode]] for t in range(packet.rank)], dtype=np.int64)
               for mode in range(3)]
    U, V, W = factors
    T = np.stack([(U[:, a, None] * V).T @ W for a in range(U.shape[1])])
    expected = np.zeros_like(T)
    for a, (h, i, k) in enumerate(keys[0]):
        for j in range(s):
            expected[a, keys[1].index((h, k, j)), keys[2].index((h, i, j))] = 1
    return T, expected


def square_tensor(square):
    n = square.n
    factors = [np.array([[int(square.factor_entry(mode, t, r, c)) for r in range(n) for c in range(n)]
                         for t in range(square.rank)], dtype=np.int64) for mode in range(3)]
    U, V, W = factors
    T = np.stack([(U[:, a, None] * V).T @ W for a in range(n * n)])
    expected = np.zeros_like(T)
    for i, k, j in product(range(n), repeat=3):
        expected[i * n + k, k * n + j, i * n + j] = 1
    return T, expected


@pytest.mark.parametrize('q', [1, 2, 3, 6, 9])
def test_local_identity(q):
    assert sp.verify_local(q) == 3 * (q + 1) ** 2 * q


def test_local_identity_rejects_a_changed_sign(monkeypatch):
    original = sp.local_factor
    monkeypatch.setattr(sp, 'local_factor', lambda q, kind, mode, choice, x:
                        -original(q, kind, mode, choice, x) if (kind, mode, choice) == (1, 2, q) else original(q, kind, mode, choice, x))
    with pytest.raises(AssertionError):
        sp.verify_local(3)


def test_audit_counts_and_rejections():
    assert sp.audit(support(2), 2)['weight_classes'] == [1, 2, 1]
    assert sp.audit(support(3), 3) == dict(m=3, k=9, products=14, partition_triples=14, weight_classes=[1, 6, 6, 1])
    # {a + b + c = 0 mod 3} has injective projections but nontrivial triangles.
    with pytest.raises(AssertionError, match='cross triple'):
        sp.audit(support(3, LATIN), 3)
    # Weight class 0 row is (K2, K0, K1); swapping A and C puts C on kind 0.
    K = [7, 7 << 3, 7 << 6]
    assert sp.audit([(K[2], K[0], K[1])], 3)['products'] == 1
    with pytest.raises(AssertionError, match='letter'):
        sp.audit([(K[1], K[0], K[2])], 3)
    with pytest.raises(AssertionError):
        sp.audit(support(3) + support(3)[:1], 3)


def test_small_packet_full_tensor():
    """N = 6, q = 2: all 4 rows of <4>, every coefficient over the row blocks, 729 terms."""
    packet = sp.StrassenPacket(2, support(2), 2)
    assert packet.rank == 729 and len(packet.rows) == 4
    T, expected = packet_tensor(packet)
    assert np.array_equal(T, expected)
    rng = np.random.default_rng(0)
    s = packet.q ** packet.m
    for _ in range(40):
        h, i, k, j = (int(x) for x in rng.integers(0, [4, s, s, s]))
        v = int(rng.integers(0, 4))
        x, y, z = packet.matrix_variable(h, 0, i, k), packet.matrix_variable(v, 1, k, j), packet.matrix_variable(h, 2, i, j)
        assert sp.tensor_entry(packet, x, y, z) == int(v == h)
    assert sp.packet_identity(packet, samples=12)['passed']
    assert not sp.packet_identity(packet, samples=3, perturb=1)['passed']


def test_power9_support_at_q1_and_q2():
    rows = support(3)
    T, expected = packet_tensor(sp.StrassenPacket(1, rows, 3))
    assert T.shape == (14, 14, 14) and np.array_equal(T, expected)
    assert sp.packet_identity(sp.StrassenPacket(2, rows, 3), samples=30, seed=3)['passed']
    # Negative control: without the RS condition the q = 1 tensor has off-diagonal entries.
    T, expected = packet_tensor(unaudited(1, support(3, LATIN), 3))
    assert T.shape == (20, 20, 20) and not np.array_equal(T, expected)


@pytest.mark.parametrize('outer,leftover', [(dict(kind='lille', n=2, rank=7), 0), (dict(kind='classical', n=2), 2)])
def test_mixed_square_full_tensor(outer, leftover):
    """Size 8 = 2 * 4: a partly used final packet, unused rows with zero input, leftover leaf calls."""
    spec = sp.StrassenPacket(2, support(2), 2).spec()
    square = MixedSquare(outer, dict(kind='classical', n=4), [spec, spec], leftover_products=leftover)
    assert square.packets[0] is square.packets[1] and square.used == [4, 4 - (1 if leftover == 0 else 2)]
    assert square.rank == 2 * 729 + leftover * 64
    restored = load_certificate(square.certificate())
    assert restored.certificate() == square.certificate()
    T, expected = square_tensor(restored)
    assert np.array_equal(T, expected)


def test_unknown_local_tensor_rejected():
    spec = dict(sp.StrassenPacket(1, support(2), 2).spec(), local_tensor='strassen-completion-v0')
    with pytest.raises(AssertionError, match='local tensor'):
        MixedSquare(dict(kind='classical', n=1), dict(kind='classical', n=1), [spec], leftover_products=0)


def test_published_strassen_entries():
    e = {r['n']: r for r in cwschemes.catalogue()}
    assert int(e['8388608']['rank']) == 40 * 9 ** 18 == 6003785411879964840
    assert int(e['13947137604']['rank']) == int(e['13468840704']['rank']) == 2 * 10 ** 27
    s = cwschemes.load(8388608)
    engine = s._engine
    assert [len(p.rows) for p in engine.packets] == [358] * 40 and engine.leaves == 0
    assert sum(engine.used) == 14196 and 40 * 358 - sum(engine.used) == 124
    packet = engine.packets[-1]
    assert packet.witness['weight_classes'] == [1, 18, 90, 140, 90, 18, 1]
    # A nonzero coefficient: final call 14195 sits on row 233 of the last packet.
    h, mode = 233, 0
    br, bc = next((r, c) for r in range(32) for c in range(32) if engine.outer.factor_entry(mode, 14195, r, c))
    ir, ic = 12345, 67890
    v = packet.matrix_variable(h, mode, ir, ic)
    digits = [max(x - 1, 0) for x in v]
    term = engine.ends[-2] + sum(d * 9 ** (17 - i) for i, d in enumerate(digits))
    expected = sum((engine.outer.factor_entry(mode, engine.starts[-1] + g, br, bc)
                    * packet.factor_entry(mode, term - engine.ends[-2], packet.matrix_variable(g, mode, ir, ic))
                    for g in range(engine.used[-1])), Fraction(0))
    got = s.coefficient('U', term, br * 262144 + ir, bc * 262144 + ic)
    assert got == expected and got != 0
    padded = cwschemes.load(13468840704)._engine
    assert padded.full_n == 13947137604 == 36 * 9 ** 9 and padded.used == [9842, 9824]
    assert padded.outer.rank == 19666 and padded.rank ** 10 < padded.n ** 27
