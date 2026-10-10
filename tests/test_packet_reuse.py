"""Repeated packet loads retain all batch wiring and exact coefficients."""
from copy import deepcopy
from fractions import Fraction
from importlib import import_module

import pytest


@pytest.mark.parametrize('q', [1, 2])
def test_repeated_packets_preserve_partial_batch_and_fallback(q, monkeypatch):
    from cwschemes._engine import general_square as gs

    spec = dict(q=q, m=1, rows=[[1, 2, 4], [2, 4, 1]], patched=True)
    # Different mapping insertion order must not prevent structural reuse.
    reordered = {key: deepcopy(spec[key]) for key in reversed(spec)}
    outer = dict(kind='classical', n=2)
    leaf = dict(kind='classical', n=q)
    made = []
    constructor = gs.CompactPacket

    def counted(*args, **kwargs):
        packet = constructor(*args, **kwargs)
        made.append(packet)
        return packet

    monkeypatch.setattr(gs, 'CompactPacket', counted)
    cached = gs.MixedSquare(outer, leaf, [spec, reordered], leftover_products=5)
    assert len(made) == 1
    assert cached.packets[0] is cached.packets[1]
    # Ignored metadata distinguishes the cache keys and reproduces the previous
    # independent construction of otherwise identical real packet objects.
    separate = gs.MixedSquare(outer, leaf,
                              [dict(spec, regression_instance=0),
                               dict(spec, regression_instance=1)],
                              leftover_products=5)
    assert len(made) == 3
    assert separate.packets[0] is not separate.packets[1]
    assert cached.starts == separate.starts == [0, 2]
    assert cached.used == separate.used == [2, 1]
    assert cached.ends == separate.ends
    assert cached.covered == separate.covered == 3
    assert cached.leaves == separate.leaves == 5
    assert cached.rank == separate.rank
    assert cached.certificate() == separate.certificate()
    for mode in range(3):
        for term in range(cached.rank):
            for row in range(cached.n):
                for col in range(cached.n):
                    assert cached.factor_entry(mode, term, row, col) == separate.factor_entry(mode, term, row, col)
    # Check the whole resulting tensor over Q, including every zero coefficient.
    factors = [[[cached.factor_entry(mode, term, row, col)
                 for row in range(cached.n) for col in range(cached.n)]
                for term in range(cached.rank)] for mode in range(3)]
    for a in range(cached.n ** 2):
        i, k = divmod(a, cached.n)
        for b in range(cached.n ** 2):
            kb, j = divmod(b, cached.n)
            for c in range(cached.n ** 2):
                got = sum((factors[0][t][a] * factors[1][t][b] * factors[2][t][c]
                           for t in range(cached.rank)), Fraction(0))
                assert got == int(k == kb and c == i * cached.n + j)


def test_repeated_support_and_completion_are_audited_once(monkeypatch):
    verifier = import_module('cwschemes.verify')
    from cwschemes._engine import restrictions

    spec = dict(q=2, m=1, rows=[[1, 2, 4], [2, 4, 1]])
    calls = {'support': 0, 'completion': 0}
    original_audit = restrictions.audit
    original_completion = restrictions.verify_completion

    def audited(*args, **kwargs):
        calls['support'] += 1
        return original_audit(*args, **kwargs)

    def completed(*args, **kwargs):
        calls['completion'] += 1
        return original_completion(*args, **kwargs)

    monkeypatch.setattr(restrictions, 'audit', audited)
    monkeypatch.setattr(restrictions, 'verify_completion', completed)
    monkeypatch.setattr(verifier, '_AUDITED', set())
    assert verifier.audit_supports({'packets': [deepcopy(spec) for _ in range(100)]}) == 1
    assert calls == {'support': 1, 'completion': 1}
    assert verifier.audit_supports(deepcopy(spec)) == 0
    assert calls == {'support': 1, 'completion': 1}
