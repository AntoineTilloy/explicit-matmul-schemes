#!/usr/bin/env python3
"""Exact square restriction after three independent tensor-vertex permutations."""
from bisect import bisect_right
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
from crop_accessor import FrozenBasis, select_bit

EDGES=((0,1),(1,2),(0,2))


def allocations(outer,inner,n,orders,inner_orders=None):
    positions=[np.argsort(order) for order in orders]
    inner_positions=[np.argsort(order) for order in inner_orders] if inner_orders is not None else [np.arange(inner.n)]*3
    blocks,delta=divmod(n,inner.n);all_inner=(1<<inner.rank)-1;masks=[]
    for mode in range(3):
        partial=[0,0,0]
        for t in range(inner.rank):
            support=inner.support(mode,t);rv,cv=EDGES[mode];rr=inner_positions[rv][support//inner.n];cc=inner_positions[cv][support%inner.n]
            for j,alive in enumerate((np.any(rr<delta),np.any(cc<delta),np.any((rr<delta)&(cc<delta)))):
                if alive:partial[j]|=1<<t
        masks.append([((partial[0] if c&1 else 0)|(partial[1] if c&2 else 0)|(partial[2] if c&4 else 0)) for c in range(8)]+[all_inner])
    keys=[];kept={};prefix=[0]
    for t in range(outer.rank):
        codes=[]
        for mode,(row_vertex,col_vertex) in enumerate(EDGES):
            code=0
            for k in outer.support(mode,t):
                row,col=divmod(int(k),outer.n);row=int(positions[row_vertex][row]);col=int(positions[col_vertex][col])
                if row<blocks and col<blocks:code=8;break
                if row==blocks and col<blocks:code|=1
                elif row<blocks and col==blocks:code|=2
                elif row==blocks and col==blocks:code|=4
            codes.append(code)
        key=tuple(codes)
        if key not in kept:kept[key]=masks[0][key[0]]&masks[1][key[1]]&masks[2][key[2]]
        keys.append(key);prefix.append(prefix[-1]+kept[key].bit_count())
    return keys,kept,prefix


class VertexCrop:
    def __init__(self,path):
        self.path=Path(path);s=json.loads(self.path.read_text());self.spec=s
        assert s['format']=='exact-three-vertex-crop-v1';self.n=s['tensor'][0];assert s['tensor']==[self.n]*3
        self.field=s['coefficient_field'];self.orders=s['vertex_orders'];assert len(self.orders)==3
        self.outer=FrozenBasis(self.path.parent/s['outer']['file'],s['outer']['sha256']);self.inner=FrozenBasis(self.path.parent/s['inner']['file'],s['inner']['sha256'])
        assert self.field==self.outer.field==self.inner.field
        assert all(sorted(order)==list(range(self.outer.n)) for order in self.orders)
        assert 0<self.n<=self.outer.n*self.inner.n
        self.positions=[np.argsort(order) for order in self.orders]
        self.inner_orders=s.get('inner_vertex_orders',[list(range(self.inner.n)) for _ in range(3)])
        assert len(self.inner_orders)==3 and all(sorted(order)==list(range(self.inner.n)) for order in self.inner_orders)
        self.inner_positions=[np.argsort(order) for order in self.inner_orders]
        self.keys,self.kept,self.prefix=allocations(self.outer,self.inner,self.n,self.orders,self.inner_orders)
        self.rank=self.prefix[-1];assert self.rank==s['rank']
        bp=self.path.parent/s['bookkeeping'];assert hashlib.sha256(bp.read_bytes()).hexdigest()==s['bookkeeping_sha256']
        with np.load(bp,allow_pickle=False) as z:
            assert np.array_equal(z['outer_indptr'],np.asarray(self.prefix,np.int64))
            signature=[tuple(map(int,row)) for row in z['class_signature']]
            assert len(signature)==len(set(signature))==len(self.kept) and set(signature)==set(self.kept)
            assert z['outer_class'].shape==(self.outer.rank,)
            assert [signature[int(j)] for j in z['outer_class']]==self.keys
            assert z['class_inner_kept'].shape==(len(signature),(self.inner.rank+7)//8)
            for j,key in enumerate(signature):assert int.from_bytes(z['class_inner_kept'][j].tobytes(),'little')==self.kept[key]

    def locate(self,t):
        assert 0<=t<self.rank
        ot=bisect_right(self.prefix,t)-1
        return ot,select_bit(self.kept[self.keys[ot]],t-self.prefix[ot])

    def coefficient(self,mode,t,row,col):
        assert mode in (0,1,2)
        if not(0<=row<self.n and 0<=col<self.n):return 0
        ot,it=self.locate(t);ro,ri=divmod(row,self.inner.n);co,ci=divmod(col,self.inner.n);rv,cv=EDGES[mode]
        a=self.outer.coefficient(mode,ot,self.orders[rv][ro],self.orders[cv][co]);b=self.inner.coefficient(mode,it,self.inner_orders[rv][ri],self.inner_orders[cv][ci])
        return a*b if self.field=='Q' else (a*b)%2

    factor_entry=coefficient

    def term(self,t):
        ot,it=self.locate(t);out=[]
        for mode,(rv,cv) in enumerate(EDGES):
            entries=[]
            for ia,na,da in self.outer.sparse(mode,ot):
                ra,ca=divmod(int(ia),self.outer.n);ra=int(self.positions[rv][ra]);ca=int(self.positions[cv][ca])
                for ib,nb,db in self.inner.sparse(mode,it):
                    rb,cb=divmod(int(ib),self.inner.n);rb=int(self.inner_positions[rv][rb]);cb=int(self.inner_positions[cv][cb]);row=ra*self.inner.n+rb;col=ca*self.inner.n+cb
                    if row>=self.n or col>=self.n:continue
                    value=Fraction(int(na)*int(nb),int(da)*int(db))
                    if self.field=='F2':assert value.denominator%2;value=Fraction(value.numerator%2)
                    if value:entries.append((row,col,value.numerator,value.denominator))
            assert entries;out.append(sorted(entries))
        return out
