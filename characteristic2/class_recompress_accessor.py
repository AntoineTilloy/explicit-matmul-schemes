#!/usr/bin/env python3
"""Exact compact GF(2) recompression, grouped by inner restriction profile."""
from bisect import bisect_right
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
from crop_accessor import Crop


def inner_profiles(inner,delta,positions=None):
    key=np.zeros(inner.rank,np.uint16)
    for mode,axis in enumerate('uvw'):
        p=inner.data[axis+'_indptr'];idx=inner.data[axis+'_indices'];rows=np.repeat(np.arange(inner.rank),np.diff(p));rr,cc=idx//inner.n,idx%inner.n
        if positions is not None:
            rv,cv=((0,1),(1,2),(0,2))[mode];rr,cc=positions[rv][rr],positions[cv][cc]
        flags=np.zeros(inner.rank,np.uint16)
        for bit,condition in ((1,rr<delta),(2,cc<delta),(4,(rr<delta)&(cc<delta))):
            alive=np.zeros(inner.rank,np.uint8);alive[rows[condition]]=1;flags|=alive.astype(np.uint16)*bit
        key|=flags<<(3*mode)
    return key


def outer_profile(c,key):
    blocks,delta=divmod(c.n,c.inner.n);forms=[]
    for t in range(c.outer.rank):
        triple=[]
        for mode in range(3):
            flags=(key>>(3*mode))&7;mask=0
            for idx in c.outer.support(mode,t):
                r,co=divmod(int(idx),c.outer.n)
                if hasattr(c,'positions'):
                    rv,cv=((0,1),(1,2),(0,2))[mode];r=int(c.positions[rv][r]);co=int(c.positions[cv][co])
                else:r=int(c.position[r]);co=int(c.position[co])
                alive=False
                if r<blocks and co<blocks:alive=True
                elif r==blocks and co<blocks:alive=bool(flags&1)
                elif r<blocks and co==blocks:alive=bool(flags&2)
                elif r==blocks and co==blocks:alive=bool(flags&4)
                if alive:mask|=1<<(r*c.outer.n+co)
            triple.append(mask)
        forms.append(tuple(triple))
    return forms


class ClassRecompressed:
    def __init__(self,path):
        path=Path(path);s=json.loads(path.read_text());assert s['format']=='exact-f2-class-recompression-v1' and s['coefficient_field']=='F2';self.spec=s
        bp=path.parent/s['base'];assert hashlib.sha256(bp.read_bytes()).hexdigest()==s['base_sha256']
        if json.loads(bp.read_text())['format']=='exact-three-vertex-crop-v1':
            from vertex_crop_accessor import VertexCrop
            self.base=VertexCrop(bp)
        else:self.base=Crop(bp)
        self.n=self.base.n;self.field='F2';self.rank=s['rank'];assert self.base.field=='F2' and s['tensor']==[self.n]*3
        ap=path.parent/s['allocation'];assert hashlib.sha256(ap.read_bytes()).hexdigest()==s['allocation_sha256']
        with np.load(ap,allow_pickle=False) as z:self.a={k:z[k] for k in z.files}
        self.profile=inner_profiles(self.base.inner,self.n%self.base.inner.n,getattr(self.base,'inner_positions',None));assert np.array_equal(self.profile,self.a['inner_profile'])
        keys=sorted(set(map(int,self.profile)));assert keys==list(map(int,self.a['class_keys'])) and len(s['classes'])==len(keys)
        ptr=self.a['outer_indptr'];assert ptr.shape==(len(keys)+1,) and ptr[0]==0 and np.all(np.diff(ptr)>=0) and ptr[-1]==len(self.a['outer_ids'])
        self.keys=keys;self.outer_forms=[];self.inner_ids=[];self.outer_ids=[];self.prefix=[0];original_rank=0
        for j,key in enumerate(keys):
            info=s['classes'][j];multiplicity=int(np.sum(self.profile==key));forms=outer_profile(self.base,key);live={t for t,f in enumerate(forms) if all(f)};initial=len(live)
            assert info['key']==key and info['inner_multiplicity']==multiplicity and info['old_outer_calls']==initial
            original_rank+=initial*multiplicity;merged={};next_id=self.base.outer.rank
            for event in s['classes'][j]['merges']:
                members=event['members'];assert len(set(members))==len(members)>=2 and all(t in live for t in members);pair=event['pair'];assert len(set(pair))==2 and set(pair)<={0,1,2};other=next(k for k in range(3) if k not in pair)
                ff=[forms[t] if t<self.base.outer.rank else merged[t] for t in members];assert all(all(f[p]==ff[0][p] for p in pair) for f in ff)
                out=list(ff[0]);v=0
                for f in ff:v^=f[other]
                out[other]=v;live.difference_update(members)
                if v:assert event['new_id']==next_id;merged[next_id]=tuple(out);live.add(next_id);next_id+=1
                else:assert event['new_id'] is None
            p=self.a['outer_indptr'];ids=self.a['outer_ids'][int(p[j]):int(p[j+1])];assert len(set(map(int,ids)))==len(ids) and set(map(int,ids))==live
            assert info['new_outer_calls']==len(ids) and info['weighted_saving']==multiplicity*(initial-len(ids))
            inner=np.flatnonzero(self.profile==key);self.inner_ids.append(inner);self.outer_ids.append(ids);self.outer_forms.append((forms,merged));self.prefix.append(self.prefix[-1]+len(ids)*len(inner))
        assert original_rank==self.base.rank and self.prefix[-1]==self.rank

    def locate(self,t):
        assert 0<=t<self.rank;c=bisect_right(self.prefix,t)-1;inner_pos,outer_pos=divmod(t-self.prefix[c],len(self.outer_ids[c]));return c,int(self.inner_ids[c][inner_pos]),int(self.outer_ids[c][outer_pos])

    def _outer(self,c,t):
        forms,merged=self.outer_forms[c];return forms[t] if t<self.base.outer.rank else merged[t]

    def coefficient(self,mode,t,row,col):
        if not (0<=row<self.n and 0<=col<self.n):return 0
        c,it,ot=self.locate(t);ro,ri=divmod(row,self.base.inner.n);co,ci=divmod(col,self.base.inner.n);mask=self._outer(c,ot)[mode]
        if hasattr(self.base,'inner_orders'):
            rv,cv=((0,1),(1,2),(0,2))[mode];ri,ci=self.base.inner_orders[rv][ri],self.base.inner_orders[cv][ci]
        return self.base.inner.coefficient(mode,it,ri,ci) if mask>>(ro*self.base.outer.n+co)&1 else 0

    factor_entry=coefficient

    def term(self,t):
        c,it,ot=self.locate(t);outer=self._outer(c,ot);result=[]
        for mode in range(3):
            mask=outer[mode];entries=[]
            while mask:
                low=mask&-mask;r,co=divmod(low.bit_length()-1,self.base.outer.n);mask-=low
                for idx,num,den in self.base.inner.sparse(mode,it):
                    ri,ci=divmod(int(idx),self.base.inner.n)
                    if hasattr(self.base,'inner_positions'):
                        rv,cv=((0,1),(1,2),(0,2))[mode];ri=int(self.base.inner_positions[rv][ri]);ci=int(self.base.inner_positions[cv][ci])
                    row=r*self.base.inner.n+ri;col=co*self.base.inner.n+ci
                    if row<self.n and col<self.n and int(num)%2:assert int(den)%2;entries.append((row,col,1,1))
            assert entries;result.append(sorted(entries))
        return result
