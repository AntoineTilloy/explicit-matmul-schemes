from fractions import Fraction
import gzip
import json
from pathlib import Path
import pytest
import cwschemes
from cwschemes import _resources


def test_catalogue_and_exact_threshold():
    entries=cwschemes.catalogue()
    assert len(entries)==41 and len({r['n'] for r in entries})==41
    d=next(r for r in entries if r['n']=='5873299070')
    n,r=int(d['n']),int(d['rank'])
    assert (n-1)**27<=r**10<n**27


def test_corrupt_object_is_rejected(tmp_path,monkeypatch):
    h=cwschemes.catalogue()[0]['root_object']
    (tmp_path/'objects').mkdir()
    (tmp_path/'objects'/f'{h}.gz').write_bytes(gzip.compress(b'{}'))
    _resources.index()
    monkeypatch.setattr(_resources,'DATA',tmp_path)
    with pytest.raises(ValueError,match='integrity'): _resources.raw_object(h)


def test_bad_support_is_rejected():
    from cwschemes._engine.restrictions import audit
    with pytest.raises(AssertionError):audit([[1,2,4],[1,4,2]],1)


def test_exact_small_multiplication():
    for p in _resources.index()['paths']:
        if p.startswith('seeds/'): _resources.materialize(p)
    from cwschemes._engine.general_square import Ordinary
    s=Ordinary({'kind':'lille','n':2,'rank':7})
    A=[[Fraction(2,3),-1],[4,5]];B=[[7,8],[9,Fraction(1,2)]]
    got=[[Fraction(0) for _ in range(2)] for _ in range(2)]
    for t in range(s.rank):
        a=sum(s.factor_entry(0,t,i,j)*A[i][j] for i in range(2) for j in range(2))
        b=sum(s.factor_entry(1,t,i,j)*B[i][j] for i in range(2) for j in range(2))
        for i in range(2):
            for j in range(2):got[i][j]+=s.factor_entry(2,t,i,j)*a*b
    assert got==[[sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def test_public_index_validation():
    from cwschemes import Scheme
    class Tiny:
        n=2;rank=7
        def factor_entry(self,*args):return Fraction(1)
    s=Scheme({},Tiny())
    with pytest.raises(IndexError):s.coefficient('U',7,0,0)
    with pytest.raises(TypeError):s.coefficient('U',0.0,0,0)
    with pytest.raises(ValueError):s.coefficients('U',[(0,0,0)]*2,max_queries=1)


def test_literal_data_is_packaged():
    assert all((_resources.DATA/'objects'/f'{h}.gz').is_file() for h in _resources.index()['objects'])


def test_explicit_small_scheme_is_exact():
    s=cwschemes.load(14)
    assert s.rank==1593
    assert s._engine.random_check(trials=1)
    assert s._engine.exact_check()['exact']


Z4_SUPPORT='z4q/support-5c5e8a6df88dc849.json'
Z4_RELATIONS='z4q/relations-304e9d51f14f180e.json'


def test_z4_quotient_catalogue_counts():
    from cwschemes._engine.z4_checks import closed_form_rank
    e={r['n']:int(r['rank']) for r in cwschemes.catalogue()}
    assert e['9826']==closed_form_rank(17)==74352484826 and e['1024']==closed_form_rank(8)==198683936
    assert e['16384']==7*closed_form_rank(16) and e['32768']==48*closed_form_rank(16)
    assert e['65536']==319*closed_form_rank(16)+3*closed_form_rank(13)
    s=cwschemes.load(9826)
    assert s.rank==74352484826 and s._engine.packets[0].rank==s.rank


def test_z4_quotient_small_exact():
    """q = 2 instance (n = 16) of the bundled Z4 support and relation objects."""
    from cwschemes._engine import z4_quotient, z4_checks
    S=z4_quotient.square(2,Z4_SUPPORT,Z4_RELATIONS);u=S.packets[0]
    assert (S.n,S.rank)==(16,15866)
    assert z4_checks.audit(u)['unquotiented_rank']=='63360'
    # coefficients recorded from the research generator (checks/z4_quotient.json compares all of them)
    assert S.factor_entry(0,2201,2,8)==Fraction(5,128)
    assert S.factor_entry(1,1931,15,14)==Fraction(81,256)
    assert S.factor_entry(2,7737,12,6)==Fraction(-66210351742441,22773104640)
    assert z4_checks.packet_identity(u,samples=8)['passed']
    assert z4_checks.orbit_lemma(u,samples=2)['passed']
    assert not z4_checks.packet_identity(u,samples=4,perturb=1)['passed']
