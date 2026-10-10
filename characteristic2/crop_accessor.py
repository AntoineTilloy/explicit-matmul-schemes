#!/usr/bin/env python3
"""Exact coefficient access to frozen two-level zero-padded restrictions.

The proof rule is Kronecker composition followed by simultaneous matrix-index
permutation and restriction. A product is removed precisely when at least one
of its three restricted linear forms vanishes. Seed correctness is recorded
separately; this module never uses a probabilistic tensor test as a proof.
"""
import argparse
from bisect import bisect_right
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np


class FrozenBasis:
    def __init__(self,path,sha256):
        self.path=Path(path);assert hashlib.sha256(self.path.read_bytes()).hexdigest()==sha256
        with np.load(self.path,allow_pickle=False) as z:self.data={k:z[k] for k in z.files}
        m=json.loads(self.data['metadata_json'].item());self.n=m['tensor'][0];self.rank=m['rank'];assert m['tensor']==[self.n]*3
        self.field=m['coefficient_field']

    def support(self,mode,t):
        axis='uvw'[mode];p=self.data[axis+'_indptr'];s=slice(int(p[t]),int(p[t+1]))
        return self.data[axis+'_indices'][s]

    def sparse(self,mode,t):
        axis='uvw'[mode];p=self.data[axis+'_indptr'];s=slice(int(p[t]),int(p[t+1]))
        return zip(self.data[axis+'_indices'][s],self.data[axis+'_numerators'][s],self.data[axis+'_denominators'][s])

    def coefficient(self,mode,t,row,col):
        for k,num,den in self.sparse(mode,t):
            if int(k)==row*self.n+col:
                if self.field=='F2':assert int(den)%2;return int(num)%2
                return Fraction(int(num),int(den))
        return 0


def select_bit(mask,k):
    """Index of the kth set bit, without scanning k preceding products."""
    assert 0<=k<mask.bit_count();offset=0
    while mask.bit_length()>256:
        cut=((mask.bit_length()+1)//2+63)//64*64
        low=mask&((1<<cut)-1);count=low.bit_count()
        if k<count:mask=low
        else:k-=count;offset+=cut;mask>>=cut
    while k:mask&=mask-1;k-=1
    return offset+(mask&-mask).bit_length()-1


class Crop:
    def __init__(self,path):
        self.path=Path(path);self.spec=json.loads(self.path.read_text());s=self.spec
        assert s['format']=='exact-two-level-crop-v1';self.n=s['tensor'][0];assert s['tensor']==[self.n]*3
        self.field=s['coefficient_field'];self.order=s['outer_order'];self.position=np.argsort(self.order)
        self.outer=FrozenBasis(self.path.parent/s['outer']['file'],s['outer']['sha256']);self.inner=FrozenBasis(self.path.parent/s['inner']['file'],s['inner']['sha256'])
        assert self.field==self.outer.field==self.inner.field
        assert sorted(self.order)==list(range(self.outer.n));assert self.n<=self.outer.n*self.inner.n
        blocks,delta=divmod(self.n,self.inner.n);all_inner=(1<<self.inner.rank)-1
        masks=[]
        for mode in range(3):
            partial=[0,0,0]
            for t in range(self.inner.rank):
                support=self.inner.support(mode,t);rr=support//self.inner.n;cc=support%self.inner.n
                for j,alive in enumerate((np.any(rr<delta),np.any(cc<delta),np.any((rr<delta)&(cc<delta)))):
                    if alive:partial[j]|=1<<t
            masks.append([((partial[0] if c&1 else 0)|(partial[1] if c&2 else 0)|(partial[2] if c&4 else 0)) for c in range(8)]+[all_inner])
        keys=[];kept={};prefix=[0]
        for t in range(self.outer.rank):
            codes=[]
            for mode in range(3):
                code=0
                for k in self.outer.support(mode,t):
                    row,col=divmod(int(k),self.outer.n);row=int(self.position[row]);col=int(self.position[col])
                    if row<blocks and col<blocks:code=8;break
                    if row==blocks and col<blocks:code|=1
                    elif row<blocks and col==blocks:code|=2
                    elif row==blocks and col==blocks:code|=4
                codes.append(code)
            key=tuple(codes)
            if key not in kept:kept[key]=masks[0][key[0]]&masks[1][key[1]]&masks[2][key[2]]
            keys.append(key);prefix.append(prefix[-1]+kept[key].bit_count())
        self.keys=keys;self.kept=kept;self.prefix=prefix;self.rank=prefix[-1]
        assert self.rank==s['rank'],(self.rank,s['rank'])
        if s.get('bookkeeping_sha256'):
            book=self.path.parent/s['bookkeeping'];assert hashlib.sha256(book.read_bytes()).hexdigest()==s['bookkeeping_sha256']
            with np.load(book,allow_pickle=False) as z:
                assert np.array_equal(z['outer_indptr'],np.asarray(self.prefix,np.int64))
                signature=[tuple(map(int,row)) for row in z['class_signature']]
                assert [signature[int(j)] for j in z['outer_class']]==self.keys
                for j,key in enumerate(signature):assert int.from_bytes(z['class_inner_kept'][j].tobytes(),'little')==self.kept[key]

    def locate(self,t):
        assert 0<=t<self.rank
        outer=bisect_right(self.prefix,t)-1;inner=select_bit(self.kept[self.keys[outer]],t-self.prefix[outer])
        return outer,inner

    def coefficient(self,mode,t,row,col):
        assert mode in (0,1,2)
        if not (0<=row<self.n and 0<=col<self.n):return 0
        ot,it=self.locate(t);ro,ri=divmod(row,self.inner.n);co,ci=divmod(col,self.inner.n)
        a=self.outer.coefficient(mode,ot,self.order[ro],self.order[co]);b=self.inner.coefficient(mode,it,ri,ci)
        return a*b if self.field=='Q' else (a*b)%2

    def term(self,t):
        ot,it=self.locate(t);result=[]
        for mode in range(3):
            entries=[]
            for ia,na,da in self.outer.sparse(mode,ot):
                ra,ca=divmod(int(ia),self.outer.n);ra=int(self.position[ra]);ca=int(self.position[ca])
                for ib,nb,db in self.inner.sparse(mode,it):
                    rb,cb=divmod(int(ib),self.inner.n);row=ra*self.inner.n+rb;col=ca*self.inner.n+cb
                    if row>=self.n or col>=self.n:continue
                    value=Fraction(int(na)*int(nb),int(da)*int(db))
                    if self.field=='F2':assert value.denominator%2;value=Fraction(value.numerator%2)
                    if value:entries.append((row,col,value.numerator,value.denominator))
            assert entries;result.append(sorted(entries))
        return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('recipe',type=Path);p.add_argument('--term',type=int,default=0);args=p.parse_args();c=Crop(args.recipe)
    print(json.dumps(dict(field=c.field,n=c.n,rank=c.rank,index=args.term,parents=c.locate(args.term),term=c.term(args.term)),indent=2))
