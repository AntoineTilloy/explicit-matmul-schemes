from fractions import Fraction
from itertools import combinations, product
import numpy as np
import pytest
import cwschemes
from cwschemes._engine import strassen_packet as sp, z4_quotient as z4, z4_checks
from cwschemes._engine.gauged_z4_square import GaugedZ4Square
from cwschemes._engine.general_square import MixedSquare, load_certificate

RS = {1: [(0, 0, 0)], 2: [(0, 0, 0), (1, 1, 1)],
       3: [(0, 0, 0), (0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 2), (2, 2, 0)]}


def support(m):
    K = [((1 << m) - 1) << (m * j) for j in range(3)]
    rows = []
    for t in range(m + 1):
        labels = [sum(1 << i for i in s) for s in combinations(range(m), t)]
        for a, b, c in RS[len(labels)]:
            s0, s1, s2 = labels[a], labels[b] << m, labels[c] << 2 * m
            rows.append((s0 | (K[2] & ~s2), (K[0] & ~s0) | s1, (K[1] & ~s1) | s2))
    return rows


def packet_tensor(packet):
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


@pytest.mark.parametrize('m,q', [(2, 1), (2, 2), (3, 1)])
def test_zero_trimmed_strassen_packet_exact_full_tensor(m, q):
    rows = support(m)
    full = sp.StrassenPacket(q, rows, m)
    trimmed = sp.StrassenPacket(q, rows, m, sp.COMPRESSION)
    assert trimmed.rank < full.rank
    T_full, expected = packet_tensor(full)
    T_trim, _ = packet_tensor(trimmed)
    assert np.array_equal(T_trim, expected) and np.array_equal(T_trim, T_full)
    assert sp.check(trimmed, samples=18)['compression'] == sp.COMPRESSION
    assert not sp.packet_identity(trimmed, samples=3, perturb=1)['passed']


def test_zero_trimmed_m3_q2_identity():
    rows = support(3)
    full = sp.StrassenPacket(2, rows, 3)
    trimmed = sp.StrassenPacket(2, rows, 3, sp.COMPRESSION)
    assert trimmed.rank < full.rank
    assert sp.check(trimmed, samples=24)['compression'] == sp.COMPRESSION
    assert not sp.packet_identity(trimmed, samples=4, perturb=1)['passed']


def test_zero_trimmed_partial_row_subset_and_mixed_square():
    rows = support(2)[:3]
    trimmed = sp.StrassenPacket(2, rows, 2, sp.COMPRESSION)
    sq = MixedSquare(dict(kind='lille', n=2, rank=7), dict(kind='classical', n=4),
                     [trimmed.spec(), trimmed.spec()], leftover_products=1)
    restored = load_certificate(sq.certificate())
    T, expected = square_tensor(restored)
    assert np.array_equal(T, expected)


def test_coord8_refined_z4_q3_unit_and_gauged_square():
    rel_q3 = 'z4q/relations-split8-q3-3285ce725baaa5a6.json'
    del_file = 'z4q/deletions-split8-17aa2a1bf25f76d9.json'
    del_sha = '17aa2a1bf25f76d93236020aa33164578138b24a269e141fe8260cae18d2535e'
    sup_path = 'z4q/support-5c5e8a6df88dc849.json'
    u = z4.Z4QuotientUnit(3, sup_path, rel_q3, 'odd')
    assert u.rank == 203319
    assert z4_checks.audit(u)['kept_patterns'] == 1883
    assert z4_checks.packet_identity(u, samples=12)['passed']
    assert z4_checks.orbit_lemma(u, samples=6)['passed']
    assert not z4_checks.packet_identity(u, samples=4, perturb=1)['passed']
    sq = GaugedZ4Square(dict(
        format='cw-z4-gauged-square-v1', coefficient_field='Q', unit=u.spec(),
        assignment=[1, 0, 5, 2, 3, 4, 6],
        input_gauges=[['1', '1', '-1', '-1', '-1', '-1', '-1'], ['1', '-1', '-1', '-1', '-1', '-1', '-1']],
        n=54, full_n=54, full_rank=203319, rank=203307,
        deleted_terms=[728, 1457, 132165, 132166, 132893, 132894, 133621, 133622, 148770, 196191, 196555, 196919],
        proof_scope='q=3 exact centered-basis zero factors',
        deletion_proof=dict(file=del_file, sha256=del_sha),
    ))
    aud = sq.audit_deletions()
    assert aud['exact'] and aud['deleted_terms'] == 12
    assert z4_checks.tensor_entries(sq, u, samples=16)['passed']


def test_new_milestone_thresholds():
    e = {int(r['n']): int(r['rank']) for r in cwschemes.catalogue()}
    n27, r27 = 240619650, e[240619650]
    assert (n27 - 1) ** 27 <= r27 ** 10 < n27 ** 27
    n265 = 58481954604699954390651197457640736657662025144749186797412271422
    r265 = e[n265]
    assert (n265 - 1) ** 53 <= r265 ** 20 < n265 ** 53


def test_current_odd_lita_indexed_matches_generator():
    import hashlib
    from pathlib import Path
    from cwschemes._engine.current_odd_lita import CurrentOddLITA
    from cwschemes._engine.sources import lita_odd_current as src
    sha = hashlib.sha256(Path(src.__file__).read_bytes()).hexdigest()
    odd9 = CurrentOddLITA(9, sha)
    assert odd9.rank == 563
    ref_terms = list(src.iter_terms(9))
    assert len(ref_terms) == odd9.rank
    for idx in (0, 15, 280, 562):
        got = odd9.term(idx)
        exp = tuple(
            {(j, i) if ax == 2 else (i, j): Fraction(num, grid)
             for k, num in form.items() for i, j in [divmod(k, 9)]}
            for ax, (form, grid) in enumerate(ref_terms[idx])
        )
        assert got == exp


def test_general_uneven_and_below1024_entries():
    for n, exp_rank in ((81, 200124), (125, 675304), (330, 10010139), (488, 28101649), (686, 70079906), (844, 122067978)):
        s = cwschemes.load(n)
        assert s.rank == exp_rank
        assert isinstance(s.coefficient('U', 0, 0, 0), Fraction)
        assert isinstance(s.coefficient('W', exp_rank - 1, n - 1, n - 1), Fraction)
