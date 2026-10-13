"""Coefficient-level checks for independently split rectangular dimensions."""
import itertools
import math
import numpy as np
import pytest

from cwschemes._engine.rectangular_square import Node


@pytest.mark.parametrize('permutation',list(itertools.permutations(range(3))))
def test_asymmetric_full_tensor(permutation):
    dims=tuple((3,4,5)[i] for i in permutation)
    p=tuple((2,3,3)[i] for i in permutation)
    q=tuple(d-pi for d,pi in zip(dims,p))
    h='93f7e9ad5997975bf9eab2005f9a2ed5cca29fa0ae92c70d713d1a5085755c5d'
    specs={'outer':dict(operation='qcsr',tensor=[2]*3,rank=7,children=[],
                       file=f'rectbases/{h}.npz',sha256=h)}
    # Actual Strassen product order: M1, M2, ..., M7.
    corners=((0,0,0),(1,0,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(0,1,0))
    for i,corner in enumerate(corners):
        d=tuple(q[k] if corner[k] else p[k] for k in range(3))
        specs[str(i)]=dict(operation='classical',tensor=d,rank=math.prod(d),children=[])
    rank=sum(specs[str(i)]['rank'] for i in range(7))
    specs['root']=dict(operation='uneven_strassen',tensor=dims,rank=rank,
                       children=list(map(str,range(7))),outer='outer',p=p,q=q)
    cache={}
    def load(key):
        if key not in cache:cache[key]=Node(specs[key],load,None)
        return cache[key]
    node=load('root'); a,b,c=dims
    got=np.zeros((a*b,b*c,a*c),dtype=np.int64)
    for t in range(rank):
        U,V,W=node.term(t)
        for (i,k),u in U.items():
            for (kk,j),v in V.items():
                for (ii,jj),w in W.items():
                    value=u*v*w
                    assert value.denominator==1
                    got[i*b+k,kk*c+j,ii*c+jj]+=int(value)
    expected=np.zeros_like(got)
    for i,k,j in itertools.product(range(a),range(b),range(c)):
        expected[i*b+k,k*c+j,i*c+j]=1
    assert np.array_equal(got,expected)
