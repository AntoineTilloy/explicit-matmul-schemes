"""Exact outer wiring, square-only zero deletion and original ordinal replay."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json

import pytest


class ClassicalSevenUnit:
    """Seven independent scalar products plus explicit zero columns."""
    q = m = 1
    input_character = 'odd'
    rows = list(range(7))
    rank = 10
    retained = (2, 3, 4, 5, 6, 7, 8)

    def matrix_variable(self, h, mode, row, col):
        assert (row, col) == (0, 0)
        return h

    def factor_entry(self, mode, term, h):
        return F(term == self.retained[h])

    def spec(self):
        return dict(q=1, m=1, input_character='odd', rank=10)


def descriptor():
    return dict(format='cw-z4-gauged-square-v1', coefficient_field='Q',
                unit=ClassicalSevenUnit().spec(), n=2, full_n=2,
                full_rank=10, rank=7, assignment=[1, 2, 4, 0, 5, 3, 6],
                input_gauges=[['2', '1', '-1', '1/3', '1', '1', '1'],
                              ['1', '1', '1', '-1', '3/2', '-1', '1']],
                deleted_terms=[0, 1, 9],
                proof_scope='Complete zero factors of the synthetic scalar unit')


def install_unit(monkeypatch):
    from cwschemes._engine import gauged_z4_square as gz
    unit = ClassicalSevenUnit()
    monkeypatch.setattr(gz, 'unit_from_spec', lambda spec: unit)
    return gz, unit


def test_square_only_deletions_preserve_whole_rational_tensor(monkeypatch):
    gz, unit = install_unit(monkeypatch)
    from cwschemes._engine.general_square import load_certificate
    s = load_certificate(descriptor())
    assert unit.rank == s.full_rank == 10 and s.rank == 7
    assert tuple(s.unit_ordinal(t) for t in range(s.rank)) == unit.retained
    tensor = {}
    for t in range(s.rank):
        factors = [[s.factor_entry(mode, t, i, j) for i in range(2) for j in range(2)]
                   for mode in range(3)]
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    tensor[a, b, c] = tensor.get((a, b, c), F(0)) + factors[0][a] * factors[1][b] * factors[2][c]
    for a in range(4):
        i, k = divmod(a, 2)
        for b in range(4):
            kb, j = divmod(b, 2)
            for c in range(4):
                assert tensor[a, b, c] == int(k == kb and c == 2 * i + j)
    restored = load_certificate(s.certificate())
    assert restored.certificate() == s.certificate()
    for mode in range(3):
        for t in range(s.rank):
            for i in range(2):
                for j in range(2):
                    assert restored.factor_entry(mode, t, i, j) == s.factor_entry(mode, t, i, j)
    assert s.audit()['exact_outer_tensor']['all_coefficients'] == 64
    assert not s.audit()['deletion_proof_replayed']


@pytest.mark.parametrize('deleted', [[], [0], [9], [0, 1, 2], [2, 3, 4], [0, 2, 9]])
def test_ordinal_boundaries_and_consecutive_deletions(deleted, monkeypatch):
    gz, unit = install_unit(monkeypatch)
    d = descriptor();d.update(deleted_terms=deleted, rank=unit.rank-len(deleted))
    s = gz.GaugedZ4Square(d)
    expected = [t for t in range(unit.rank) if t not in deleted]
    assert [s.unit_ordinal(t) for t in range(s.rank)] == expected
    for bad in (-1, s.rank, 0.0):
        with pytest.raises(IndexError):s.unit_ordinal(bad)


@pytest.mark.parametrize('deleted', [[1, 0], [1, 1], [-1], [10], [1.0]])
def test_invalid_deletion_indices_are_rejected(deleted, monkeypatch):
    gz, unit = install_unit(monkeypatch)
    d = descriptor();d.update(deleted_terms=deleted, rank=unit.rank-len(deleted))
    with pytest.raises(ValueError, match='sorted unique'):gz.GaugedZ4Square(d)


def test_proof_metadata_and_outer_corruption_are_rejected(monkeypatch):
    gz, _ = install_unit(monkeypatch)
    d = descriptor();d['deletion_proof'] = {'n': 3}
    with pytest.raises(ValueError, match='does not match'):gz.GaugedZ4Square(d)
    d = descriptor();d['assignment'][0] = d['assignment'][1]
    with pytest.raises(ValueError, match='permutation'):gz.GaugedZ4Square(d)
    d = descriptor();d['input_gauges'][0][0] = '0'
    with pytest.raises(ValueError, match='nonzero'):gz.GaugedZ4Square(d)
    d = descriptor();d['proof_scope'] = ''
    with pytest.raises(ValueError, match='proof scope'):gz.GaugedZ4Square(d)


def test_existing_unit_traversal_reaches_gauged_square(monkeypatch):
    gz, unit = install_unit(monkeypatch)
    from cwschemes._engine import z4_quotient as z4
    # Traversal deliberately recognizes the real wrapper without changing its
    # unit. The synthetic fixture only supplies the unit interface here.
    s = gz.GaugedZ4Square(descriptor())
    assert z4.units(s) == [unit]


def real_q2_proof():
    records = [
        (0,1093,[0,0,1,1,1,1,1,1,1],1),
        (127,1093,[0,0,2,2,2,2,2,2,2],1),
        (13998,10087,[1,1,1,0,1,1,1,0,1],0),
        (14125,10087,[2,2,2,0,2,2,2,0,2],0),
        (13999,10087,[1,1,1,0,1,1,1,0,2],0),
        (14124,10087,[2,2,2,0,2,2,2,0,1],0),
    ]
    return dict(coefficient_field='Q',assignment=[1,0,5,2,3,4,6],
                scales=[['1','1','-1','-1','-1','-1','-1'],
                        ['1','-1','-1','-1','-1','-1','-1'],
                        ['1','-1','1','1','1','1','1']],
                q_gates=[dict(q=2,n=16,original_rank=15854,rank=15848,saving=6,
                              deleted_terms=sorted(r[0] for r in records),
                              numeric_label_column_records=[dict(old_term_index=t,H=h,labels=lab,exact_zero_modes=[mode])
                                                            for t,h,lab,mode in records])])


def real_square(proof, tmp_path, monkeypatch):
    from cwschemes import _resources
    from cwschemes._engine.gauged_z4_square import GaugedZ4Square
    raw=json.dumps(proof,sort_keys=True).encode();path=tmp_path/'deletion-proof.json';path.write_bytes(raw)
    original=_resources.materialize
    monkeypatch.setattr(_resources,'materialize',lambda logical,*args,**kwargs:path if logical=='test-deletion-proof.json' else original(logical,*args,**kwargs))
    gate=proof['q_gates'][0]
    d=dict(format='cw-z4-gauged-square-v1',coefficient_field='Q',
           unit=dict(q=2,m=3,input_character='odd',rank=15854,
                     support='z4q/support-5c5e8a6df88dc849.json',
                     relations='z4q/relations-304e9d51f14f180e.json'),
           assignment=proof['assignment'],input_gauges=proof['scales'][:2],
           n=16,full_n=16,full_rank=15854,rank=gate['rank'],deleted_terms=gate['deleted_terms'],
           proof_scope='Exact physical centered-basis zero factors after square embedding',
           deletion_proof=dict(file='test-deletion-proof.json',sha256=hashlib.sha256(raw).hexdigest()))
    return GaugedZ4Square(d),path


def test_real_q2_deleted_columns_replay_exactly(tmp_path,monkeypatch):
    s,_=real_square(real_q2_proof(),tmp_path,monkeypatch)
    result=s.audit_deletions()
    assert result['exact'] and result['deleted_terms']==6 and result['distinct_normalized_factors']==3
    assert s.unit.rank==15854 and s.rank==15848
    assert s.audit_deletions()==result


def test_real_deletion_wrong_gauge_fails_polynomial_replay(tmp_path,monkeypatch):
    proof=real_q2_proof();proof['scales'][0][1]='-1'
    proof['scales'][2][1]='1'  # Preserve the complete outer MM2 identity.
    s,_=real_square(proof,tmp_path,monkeypatch)
    with pytest.raises(ValueError,match='nonzero over Q'):s.audit_deletions()


def test_real_bad_deletion_fails_polynomial_replay(tmp_path,monkeypatch):
    proof=real_q2_proof();gate=proof['q_gates'][0]
    gate['deleted_terms'][0]=1
    record=gate['numeric_label_column_records'][0];record.update(old_term_index=1,labels=[0,0,1,1,1,1,1,1,2])
    s,_=real_square(proof,tmp_path,monkeypatch)
    with pytest.raises(ValueError,match='nonzero over Q'):s.audit_deletions()


def test_deletion_resource_hash_is_checked(tmp_path,monkeypatch):
    s,path=real_square(real_q2_proof(),tmp_path,monkeypatch)
    path.write_text('{}')
    with pytest.raises(ValueError,match='hash mismatch'):s.audit_deletions()
