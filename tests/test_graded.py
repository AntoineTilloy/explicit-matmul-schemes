"""Finite proof and package integration checks for graded CW constructions."""
from fractions import Fraction
from itertools import product
import pytest
import cwschemes
from cwschemes._engine.graded_square import load_certificate


def tiny_certificate():
    return dict(format='graded-cw-square-v1',coefficient_field='Q',
        packet=dict(rows=['1'],potentials=[[0],[0],[0]],q=2,t=1,
            family='full_balanced_interpolation',
            expected=dict(inner_side=2,rank=27,capacity=1,power=3,interpolation_factor=1)),
        outer=dict(bases=[],strassen_copies=0,expected_side=1,expected_rank=1),
        n=2,natural_n=2,rank=27,crop=None)


def test_graded_package_full_tiny_tensor():
    engine=load_certificate(tiny_certificate())
    scheme=cwschemes.Scheme({},engine)
    variables=list(product(range(2),repeat=2))
    for (i,k),(l,j),(a,b) in product(variables,repeat=3):
        value=sum((scheme.coefficient('U',t,i,k)*scheme.coefficient('V',t,l,j)*
                   scheme.coefficient('W',t,a,b) for t in range(27)),Fraction())
        assert value==int(k==l and a==i and b==j)
    assert list(scheme.coefficient_summands('U',0,0,0))==[scheme.coefficient('U',0,0,0)]
    assert engine.audit()['coefficient_field']=='Q'


def test_graded_false_count_and_grade_rejected():
    data=tiny_certificate();data['packet']['expected']['capacity']=2
    with pytest.raises(ValueError,match='counts'):load_certificate(data)
    data=tiny_certificate();data['packet']['potentials'][0][0]=1
    with pytest.raises(ValueError,match='desired grade'):load_certificate(data)


def test_published_graded_threshold_and_query_guard():
    entries=[r for r in cwschemes.catalogue() if 'graded' in r['provenance']['construction'].lower()]
    assert len(entries)==2
    for record in entries:
        scheme=cwschemes.load(int(record['n']))
        audit=scheme._engine.audit()
        assert audit['strict_sub265'] and scheme.rank**20<scheme.n**53
        assert all(not item['engine_loaded'] for item in audit['base_references'])
        with pytest.raises(OverflowError,match='slot sum'):
            scheme.coefficient('U',0,0,0)
        # Creating the stream does not evaluate the astronomically large sum.
        stream=scheme.coefficient_summands('U',0,0,0)
        assert iter(stream) is stream
        if scheme._engine.data.get('crop') is not None:
            assert (scheme.n-1)**53<=scheme.rank**20
