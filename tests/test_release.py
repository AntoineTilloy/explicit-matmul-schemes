from fractions import Fraction
import gzip
import json
from pathlib import Path
import pytest
import cwschemes
from cwschemes import _resources


def test_catalogue_and_exact_threshold():
    entries=cwschemes.catalogue()
    assert len(entries)==97 and len({r['n'] for r in entries})==97
    d=next(r for r in entries if r['n']=='528108169')
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
    assert e['9826']==closed_form_rank(17)==74352484826
    gauged=cwschemes.load(1024)
    relation=gauged._engine.unit.packet.proof
    assert relation['certified_q']==8
    assert e['1024']==closed_form_rank(8,'odd',relation['rank_polynomials'])-72==198679192
    assert gauged._engine.unit.rank==198679264
    assert gauged._engine.audit_deletions()['deleted_terms']==72
    assert e['16384']==7*closed_form_rank(16) and e['32768']==48*closed_form_rank(16)
    assert e['65536']==929*closed_form_rank(14)
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


def test_odd_character_survives_mixed_certificate():
    """The real character changes coefficients and must survive serialization."""
    from cwschemes._engine import z4_quotient,z4_checks
    from cwschemes._engine.general_square import load_certificate
    square=z4_quotient.square(2,Z4_SUPPORT,Z4_RELATIONS,input_character='odd')
    unit=square.packets[0]
    assert square.rank==15854 and unit.orbits.removed==12
    assert z4_checks.audit(unit)['input_character']=='odd'
    assert z4_checks.orbit_lemma(unit,samples=8)['passed']
    restored=load_certificate(square.certificate())
    assert restored.packets[0].spec()==unit.spec()
    for mode in range(3):
        for t,row,col in ((0,0,0),(213,8,7),(15853,15,15)):
            assert restored.factor_entry(mode,t,row,col)==square.factor_entry(mode,t,row,col)
    with pytest.raises(ValueError,match='character'):
        z4_quotient.Z4QuotientUnit(2,Z4_SUPPORT,Z4_RELATIONS,'invalid')


def test_uneven_strassen_full_tensor():
    """Unequal 3+2 blocks: independently check every size-five coefficient."""
    import numpy as np
    from cwschemes._engine.rectangular_square import Recipe
    h='93f7e9ad5997975bf9eab2005f9a2ed5cca29fa0ae92c70d713d1a5085755c5d'
    nodes={'outer':dict(operation='qcsr',tensor=[2]*3,rank=7,children=[],
                        file=f'rectbases/{h}.npz',sha256=h)}
    dimensions=[(3,3,3),(2,3,3),(3,3,2),(2,2,3),(3,2,2),(2,3,2),(3,2,3)]
    for j,d in enumerate(dimensions):
        nodes[str(j)]=dict(operation='classical',tensor=d,rank=int(np.prod(d)),children=[])
    rank=sum(nodes[str(j)]['rank'] for j in range(7))
    nodes['root']=dict(operation='uneven_strassen',tensor=[5]*3,rank=rank,
                       children=list(map(str,range(7))),outer='outer',p=3,q=2)
    s=Recipe(dict(format='exact-rectangular-recipe-v1',coefficient_field='Q',
                  tensor=[5]*3,rank=rank,root='root',nodes=nodes,
                  primitive_identity_basis='Exact Strassen identity'))
    got=np.zeros((25,25,25),dtype=np.int64)
    for t in range(rank):
        u,v,w=s.root.term(t)
        for (i,k),a in u.items():
            for (kk,j),b in v.items():
                for (ii,jj),c in w.items():
                    value=a*b*c;assert value.denominator==1
                    got[i*5+k,kk*5+j,ii*5+jj]+=int(value)
    expected=np.zeros_like(got)
    for i in range(5):
        for k in range(5):
            for j in range(5):expected[i*5+k,k*5+j,i*5+j]=1
    assert np.array_equal(got,expected)


def test_split_with_shared_child_full_tensor():
    """A DAG can reuse the same rectangular algorithm for both split halves."""
    from cwschemes._engine.rectangular_square import Recipe
    nodes={'leaf':dict(operation='classical',tensor=[1,2,2],rank=4,children=[]),
           'root':dict(operation='split',tensor=[2,2,2],rank=8,axis=0,
                       children=['leaf','leaf'])}
    s=Recipe(dict(format='exact-rectangular-recipe-v1',coefficient_field='Q',
                  tensor=[2]*3,rank=8,root='root',nodes=nodes,
                  primitive_identity_basis='Classical matrix multiplication'))
    got={}
    for t in range(s.rank):
        u,v,w=s.root.term(t)
        for (i,k),a in u.items():
            for (kk,j),b in v.items():
                for (ii,jj),c in w.items():
                    key=(i*2+k,kk*2+j,ii*2+jj)
                    got[key]=got.get(key,Fraction(0))+a*b*c
    expected={(i*2+k,k*2+j,i*2+j):1
              for i in range(2) for k in range(2) for j in range(2)}
    assert got==expected


def test_new_rectangular_catalogue_entries():
    for n,rank in ((106,424690),(146,1040656),(432,19665128)):
        s=cwschemes.load(n)
        assert s.rank==rank and s._engine.audit()['exact_dimensions_and_counts']
        assert isinstance(s.coefficient('W',rank-1,n-1,n-1),Fraction)
