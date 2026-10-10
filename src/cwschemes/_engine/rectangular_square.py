#!/usr/bin/env python3
"""Exact rational coefficient access for explicit rectangular recipe DAGs."""
from bisect import bisect_right
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from .. import _resources

EDGES=((0,1),(1,2),(0,2))


class Node:
    def __init__(self,spec,load,root):
        self.s=spec;self.d=tuple(spec['tensor']);self.rank=spec['rank'];self.op=spec['operation'];self.children=[load(x) for x in spec.get('children',[])];self.arrays=None
        if self.op=='qcsr':
            p=_resources.materialize(spec['file']);assert hashlib.sha256(p.read_bytes()).hexdigest()==spec['sha256']
            with np.load(p,allow_pickle=False) as z:self.arrays={k:z[k] for k in z.files}
            m=json.loads(self.arrays['metadata_json'].item());assert tuple(m['tensor'])==self.d and m['rank']==self.rank and m['coefficient_field']=='Q'
        elif self.op=='classical':assert self.rank==np.prod(self.d)
        elif self.op=='product':
            a,b=self.children;assert self.d==tuple(x*y for x,y in zip(a.d,b.d)) and self.rank==a.rank*b.rank
        elif self.op=='split':
            a,b=self.children;ax=spec['axis'];assert self.rank==a.rank+b.rank
            assert all(self.d[j]==a.d[j]==b.d[j] for j in range(3) if j!=ax);assert self.d[ax]==a.d[ax]+b.d[ax]
        elif self.op=='pad':
            a,=self.children;assert self.rank==a.rank and all(x<=y for x,y in zip(self.d,a.d))
        elif self.op=='permute':
            a,=self.children;self.perm=spec['permutation'];assert sorted(self.perm)==[0,1,2] and self.d==tuple(a.d[j] for j in self.perm) and self.rank==a.rank
        elif self.op=='uneven_strassen':
            self.outer=load(spec['outer']);pp,q=spec['p'],spec['q'];assert self.d==(pp+q,)*3 and pp>=q
            self.prefix=[0]
            for t,c in enumerate(self.children):
                forms=self.outer.term(t);need=[]
                for vertex in range(3):
                    maxima=[]
                    for mode,edge in enumerate(EDGES):
                        if vertex in edge:
                            k=edge.index(vertex);maxima.append(max(pp if coordinate[k]==0 else q for coordinate in forms[mode]))
                    need.append(min(maxima))
                assert c.d==tuple(need),(t,c.d,need);self.prefix.append(self.prefix[-1]+c.rank)
            assert len(self.children)==7 and self.outer.d==(2,)*3 and self.outer.rank==7 and self.prefix[-1]==self.rank
        else:raise ValueError(self.op)

    def term(self,t):
        assert 0<=t<self.rank
        if self.op=='qcsr':
            out=[]
            for axis,edge in zip('uvw',EDGES):
                ptr=self.arrays[axis+'_indptr'];start,end=int(ptr[t]),int(ptr[t+1]);out.append({divmod(int(self.arrays[axis+'_indices'][j]),self.d[edge[1]]):F(int(self.arrays[axis+'_numerators'][j]),int(self.arrays[axis+'_denominators'][j])) for j in range(start,end)})
            return out
        if self.op=='classical':
            i,rest=divmod(t,self.d[1]*self.d[2]);k,j=divmod(rest,self.d[2]);return [{(i,k):F(1)},{(k,j):F(1)},{(i,j):F(1)}]
        if self.op=='product':
            a,b=self.children;ta,tb=divmod(t,b.rank);fa,fb=a.term(ta),b.term(tb)
            return [{(ra*b.d[edge[0]]+rb,ca*b.d[edge[1]]+cb):va*vb for (ra,ca),va in fa[mode].items() for (rb,cb),vb in fb[mode].items()} for mode,edge in enumerate(EDGES)]
        if self.op=='split':
            a,b=self.children;first=t<a.rank;child=a if first else b;forms=child.term(t if first else t-a.rank);offset=0 if first else a.d[self.s['axis']]
            return [{(r+(offset if edge[0]==self.s['axis'] else 0),c+(offset if edge[1]==self.s['axis'] else 0)):v for (r,c),v in forms[mode].items()} for mode,edge in enumerate(EDGES)]
        if self.op=='pad':
            forms=self.children[0].term(t)
            return [{(r,c):v for (r,c),v in forms[mode].items() if r<self.d[edge[0]] and c<self.d[edge[1]]} for mode,edge in enumerate(EDGES)]
        if self.op=='permute':
            forms=self.children[0].term(t);out=[]
            for edge in EDGES:
                old=(self.perm[edge[0]],self.perm[edge[1]]);mode=next(j for j,e in enumerate(EDGES) if set(e)==set(old));out.append(forms[mode] if old==EDGES[mode] else {(c,r):v for (r,c),v in forms[mode].items()})
            return out
        if self.op=='uneven_strassen':
            ot=bisect_right(self.prefix,t)-1;forms=self.children[ot].term(t-self.prefix[ot]);outer=self.outer.term(ot);pp,q=self.s['p'],self.s['q'];out=[]
            for mode in range(3):
                entries={}
                for (br,bc),ov in outer[mode].items():
                    nr,nc=(pp if br==0 else q),(pp if bc==0 else q)
                    for (r,c),v in forms[mode].items():
                        if r<nr and c<nc:
                            key=(r+br*pp,c+bc*pp);entries[key]=entries.get(key,F(0))+ov*v
                out.append({k:v for k,v in entries.items() if v})
            return out
        raise ValueError(self.op)

    def coefficient(self,mode,t,row,col):return self.term(t)[mode].get((row,col),F(0))


class Recipe:
    def __init__(self,value):
        self.spec=json.loads(Path(value).read_text()) if isinstance(value,(str,Path)) else value
        assert self.spec['format']=='exact-rectangular-recipe-v1' and self.spec['coefficient_field']=='Q';self.cache={}
        def load(key):
            if key not in self.cache:self.cache[key]=Node(self.spec['nodes'][key],load,_resources.root())
            return self.cache[key]
        self.root=load(self.spec['root']);self.n=self.root.d[0];self.rank=self.root.rank
        assert self.root.d==(self.n,)*3==tuple(self.spec['tensor']) and self.rank==self.spec['rank']

    def factor_entry(self,mode,t,row,col):return self.root.coefficient(mode,t,row,col)

    def certificate(self):return self.spec

    def audit(self):
        return dict(nodes=len(self.cache),primitive_objects=sum(x.op=='qcsr' for x in self.cache.values()),
                    exact_dimensions_and_counts=True,hash_checked_coefficients=True,
                    primitive_identity_basis=self.spec['primitive_identity_basis'])


def load_certificate(value):return Recipe(value)
