#!/usr/bin/env python3
"""SOURCE ONLY: exact tiny tensor checks plus large recipe binding, Atlas only."""
import argparse,hashlib,json,resource,time
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
from graded_lazy_scheme import GradedPacket,FrequencySequences,AggregateProfileSelection,FullBalancedSelection,GradedSquare,UnitOuter,StrassenOuter,TensorOuter
START=time.process_time()
LITA32_ROOT='a7fb19d69d8e078fb0d6d614ade10417225a186f4fbf94a8d4292f3935995e05'

def save(path,obj):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)

def check_packet(packet):
    variables=[(slot,i,j) for slot in range(packet.capacity) for i in range(packet.n) for j in range(packet.n)]
    assert len(variables)<=16 and packet.rank<=4096
    factors=[[[packet.coefficient(h,term,*v) for term in range(packet.rank)] for v in variables] for h in range(3)]
    checked=0
    for ia,A in enumerate(variables):
      for ib,B in enumerate(variables):
        pair=[u*v for u,v in zip(factors[0][ia],factors[1][ib])]
        for ic,C in enumerate(variables):
            got=sum((x*y for x,y in zip(pair,factors[2][ic])),Fraction(0))
            expected=int(A[0]==B[0]==C[0] and A[1]==C[1] and A[2]==B[1] and B[2]==C[2])
            assert got==expected,(A,B,C,got,expected)
            checked+=1
    return dict(side=packet.n,products=packet.capacity,rank=packet.rank,coefficients=checked,status='PASS exact Q')

def check_square(scheme):
    variables=list(product(range(scheme.n),repeat=2));assert scheme.n<=4 and scheme.rank<=2048
    factors=[[[scheme.coefficient(h,term,*v) for term in range(scheme.rank)] for v in variables] for h in range(3)]
    checked=0
    for ia,(i,j) in enumerate(variables):
      for ib,(k,l) in enumerate(variables):
        pair=[u*v for u,v in zip(factors[0][ia],factors[1][ib])]
        for ic,(o,p) in enumerate(variables):
            got=sum((x*y for x,y in zip(pair,factors[2][ic])),Fraction(0))
            expected=int(j==k and o==i and p==l)
            assert got==expected,(i,j,k,l,o,p,got,expected)
            checked+=1
    return dict(side=scheme.n,rank=scheme.rank,coefficients=checked,status='PASS exact Q direct-output convention')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--grading',type=Path,required=True)
    ap.add_argument('--candidate',type=Path,required=True);ap.add_argument('--grading-sha',required=True)
    ap.add_argument('--candidate-sha',required=True);ap.add_argument('--generator-sha',required=True)
    ap.add_argument('--catalogue',type=Path,required=True);ap.add_argument('--catalogue-sha',required=True)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(512<<20,512<<20))
    args.out.mkdir(parents=True,exist_ok=True)
    generator=Path(__file__).with_name('graded_lazy_scheme.py')
    assert hashlib.sha256(generator.read_bytes()).hexdigest()==args.generator_sha
    gradingraw=args.grading.read_bytes();candidateraw=args.candidate.read_bytes()
    assert hashlib.sha256(gradingraw).hexdigest()==args.grading_sha
    assert hashlib.sha256(candidateraw).hexdigest()==args.candidate_sha
    catalogue_raw=args.catalogue.read_bytes();assert hashlib.sha256(catalogue_raw).hexdigest()==args.catalogue_sha
    catalogue={x['root_object']:x for x in json.loads(catalogue_raw)}
    assert (int(catalogue[LITA32_ROOT]['n']),int(catalogue[LITA32_ROOT]['rank']))==(32,14196)
    result=dict(status='RUNNING',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        generator_sha256=args.generator_sha,grading_sha256=args.grading_sha,candidate_sha256=args.candidate_sha,catalogue_sha256=args.catalogue_sha,
        controls=[],large_recipe=None)
    def checkpoint():result['CPU_seconds']=time.process_time()-START;save(args.out/'result.json',result)
    checkpoint()
    rows=('11','23');pot=((0,2),(0,-1),(0,-1));profiles=[(2,0,0),(0,1,1)]
    sequences=FrequencySequences(profiles,(1,1))
    assert sequences.count==2 and {sequences.unrank(i) for i in range(2)}=={(0,1),(1,0)}
    aggregate=AggregateProfileSelection(rows,2)
    brute=[]
    for seqs in product(tuple(product(range(2),repeat=2)),repeat=3):
        totals=[tuple(sum(profiles[row][h] for row in seq) for h in range(3)) for seq in seqs]
        if totals[0]==totals[1]==totals[2]:brute.append(tuple((row,o) for o,seq in enumerate(seqs) for row in seq))
    assert len(brute)==aggregate.capacity==10
    assert {aggregate.slot_rows(i) for i in range(aggregate.capacity)}==set(brute)
    result['controls'].append(dict(kind='complete_sequence_unranking',frequency_count=2,aggregate_count=10,status='PASS'));checkpoint()
    full_rows=('112','123','233');full=FullBalancedSelection(full_rows,2);brute_full=[]
    for sequence in product(range(len(full_rows)),repeat=6):
        selected=tuple((row,pos//2) for pos,row in enumerate(sequence))
        counts=[0,0,0]
        for row,orientation in selected:
            for ch in full_rows[row]:counts[(int(ch)-1-orientation)%3]+=1
        if counts==[6,6,6]:brute_full.append(selected)
    assert full.capacity==len(brute_full)
    assert [full.slot_rows(i) for i in range(full.capacity)]==brute_full
    result['controls'].append(dict(kind='full_balanced_complete_lexicographic_unranking',count=full.capacity,enumerated_row_tuples=3**6,status='PASS'));checkpoint()
    ordinary=GradedPacket(('1',),((0,),(0,),(0,)),2,(1,),family='profile_interpolation')
    result['controls'].append(dict(kind='q2_full_packet',**check_packet(ordinary)));checkpoint()
    multi=GradedPacket(rows,pot,1,family='aggregate_profile_interpolation',t=1,allow_q1_control=True)
    result['controls'].append(dict(kind='q1_complete_multislot_interpolation_packet',**check_packet(multi)));checkpoint()
    full_multi=GradedPacket(rows,pot,1,family='full_balanced_interpolation',t=1,allow_q1_control=True)
    result['controls'].append(dict(kind='q1_complete_full_balanced_packet',**check_packet(full_multi)));checkpoint()
    zero=GradedPacket(rows,pot,1,(0,1),family='actual_row_grade_zero',allow_q1_control=True)
    result['controls'].append(dict(kind='q1_complete_grade_zero_packet',**check_packet(zero)));checkpoint()
    result['controls'].append(dict(kind='strassen_base',**check_square(StrassenOuter())));checkpoint()
    result['controls'].append(dict(kind='outer_tensor_product',**check_square(TensorOuter((StrassenOuter(),StrassenOuter())))));checkpoint()
    multi.coefficient=lru_cache(maxsize=None)(multi.coefficient)
    square=GradedSquare(multi,StrassenOuter())
    result['controls'].append(dict(kind='full_batched_square_composition',**check_square(square)));checkpoint()
    crop=GradedSquare(multi,StrassenOuter(),crop=1)
    result['controls'].append(dict(kind='full_zero_padding_crop',**check_square(crop)));checkpoint()
    candidate_data=json.loads(candidateraw);candidate=candidate_data['best_natural']
    grading_data=json.loads(gradingraw)
    puzzle=candidate.get('puzzle',candidate_data.get('puzzle'))
    target=next(x for x in grading_data['targets'] if x['puzzle']==puzzle)
    assert target['status']=='CERTIFIED_POSITIVE_GRADING'
    if 'family' in candidate and 'frequencies' in candidate:
        packet=GradedPacket(target['rows'],target['certificate']['potentials'],candidate['q'],candidate['frequencies'],family=candidate['family'])
        assert packet.sequences.count==candidate['sequence_count']
    else:
        family=('full_balanced_interpolation' if candidate['construction'].startswith('All balanced row tuples')
                else 'aggregate_profile_interpolation')
        packet=GradedPacket(target['rows'],target['certificate']['potentials'],candidate['q'],family=family,t=candidate['t'])
    assert packet.rank==candidate['rank'] and packet.capacity==candidate['square_capacity'] and packet.n==candidate['inner_side']
    i,j=candidate['outer']['lita32_copies'],candidate['outer']['strassen2_copies']
    h=candidate['outer'].get('public_base_copies',0)
    bases=[dict(n=32,rank=14196,copies=i,root_object=LITA32_ROOT)]
    outer_side=32**i*2**j;outer_rank=14196**i*7**j
    if h:
        root=candidate['outer']['public_root'];base=catalogue[root]
        assert base['field']=='Q'
        if 'public_outer' in candidate_data:assert base==candidate_data['public_outer']
        bn,br=int(base['n']),int(base['rank'])
        outer_side*=bn**h;outer_rank*=br**h
        bases.append(dict(n=bn,rank=br,copies=h,root_object=root))
    assert outer_side==candidate['outer']['side'] and outer_rank==candidate['outer']['rank']<=packet.capacity
    assert outer_side*packet.n==candidate['natural_n']
    assert packet.rank**20<candidate['natural_n']**53
    crop_n=candidate['padded_n'];assert (crop_n-1)**53<=packet.rank**20<crop_n**53
    for slot in (0,outer_rank-1,packet.capacity-1):
        zeros,grades=packet.slot_geometry(slot)
        assert all(zeros.count(h)==packet.axis for h in range(3)) and sum(grades)==0
        for mode in range(3):
            word,grade=packet.word(mode,slot,packet.n-1,0)
            assert word.count(0)==packet.axis and set(word)<=set(range(packet.q+1))
    # One computable large packet coefficient per mode, NOT a large composed
    # coefficient: the latter may require an astronomical finite slot sum.
    packet_samples=[str(packet.coefficient(mode,0,0,0,0)) for mode in range(3)]
    result['large_recipe']=dict(status='EXACT_RECIPE_BOUND',packet=packet.binding(),
        outer=dict(kind='tensor_product',catalogue_bases=bases,strassen2_copies=j,side=outer_side,rank=outer_rank,
            dependency='Hash-bound existing catalogue certificates; this control does not re-prove those bases'),
        natural_n=candidate['natural_n'],padded_n=crop_n,rank=packet.rank,
        strict_sub265='PASS integer twentieth/fifty-third powers',large_packet_coefficient_samples=packet_samples,
        limits='Large square tensor not expanded; large composed coefficients not evaluated; each is defined by a potentially enormous finite slot sum')
    result['status']='COMPLETE';checkpoint()

if __name__=='__main__':main()
