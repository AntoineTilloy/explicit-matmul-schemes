"""Exact rational coefficient DEFINITION for finite graded CW squares.

Pure Python/Fraction core. No arrays proportional to the target rank/side.
A composed coefficient is a finite sum over outer calls: this can be enormous
and is NOT claimed to be fast random access. The default query guard refuses
large slot sums; max_slot_terms=None explicitly requests the full finite sum.
Scientific checks must run on Atlas after review.
"""
from collections import Counter
from fractions import Fraction
from math import comb, factorial, prod
import hashlib,json


def zero_masks(rows):
    return [[sum(1<<j for j,ch in enumerate(row) if ch==str(h+1)) for row in rows] for h in range(3)]


def replay_grading(rows,potentials):
    s=len(rows);w=len(rows[0])
    if not all(len(row)==w and set(row)<=set('123') for row in rows):raise ValueError('invalid puzzle')
    masks=zero_masks(rows)
    if any(len(set(mode))!=s for mode in masks):raise ValueError('mode masks not injective')
    if len(potentials)!=3 or any(len(mode)!=s or any(type(x)is not int for x in mode) for mode in potentials):raise ValueError('invalid potentials')
    degree=0;edges=0
    for i,A in enumerate(masks[0]):
      for j,B in enumerate(masks[1]):
        for k,C in enumerate(masks[2]):
          if A&B==A&C==B&C:
            edges+=1;grade=potentials[0][i]+potentials[1][j]+potentials[2][k]
            if i==j==k:
                if grade!=0:raise ValueError('nonzero desired grade')
            elif grade<=0:raise ValueError('unwanted nonpositive grade')
            degree=max(degree,grade)
    return dict(degree=degree,structural_edges=edges,mode_injective=True)


class FrequencySequences:
    """Lexicographic unranking of actual row sequences with type frequencies."""
    def __init__(self,row_types,frequencies):
        self.types=sorted(set(row_types));self.groups=[[i for i,t in enumerate(row_types) if t==key] for key in self.types]
        self.row_group={row:g for g,group in enumerate(self.groups) for row in group}
        self.frequencies=tuple(frequencies);self.t=sum(frequencies)
        if len(frequencies)!=len(self.types) or any(type(x)is not int or x<0 for x in frequencies):raise ValueError('bad frequencies')
        numerator=factorial(self.t)*prod(len(group)**v for group,v in zip(self.groups,frequencies))
        denominator=prod(factorial(v) for v in frequencies)
        if numerator%denominator:raise AssertionError('nonintegral count')
        self.count=numerator//denominator;self.row_count=len(row_types)
    def unrank(self,index):
        if not 0<=index<self.count:raise IndexError('sequence rank')
        remaining=list(self.frequencies);N=self.t;count=self.count;out=[]
        while N:
            for row in range(self.row_count):
                g=self.row_group[row]
                if not remaining[g]:continue
                numerator=count*remaining[g];denominator=N*len(self.groups[g])
                if numerator%denominator:raise AssertionError('nonintegral branch count')
                block=numerator//denominator
                if index>=block:index-=block;continue
                out.append(row);remaining[g]-=1;N-=1;count=block;break
            else:raise AssertionError('unrank failed')
        if index:raise AssertionError('unrank remainder')
        return tuple(out)


class AggregateProfileSelection:
    """Exact C_t(profile)^3 selection; DP counts actual row sequences.

    All three orientations use the same aggregate profile, not necessarily
    the same type frequencies. Capacity is the SUM of cubes, not the cube
    of the total number of rectangular constituents.
    """
    def __init__(self,rows,t):
        if type(t)is not int or t<1:raise ValueError('positive t')
        self.t=t;self.width=len(rows[0]);self.profiles=[tuple(row.count(c) for c in '123') for row in rows]
        histogram=Counter(self.profiles);self.dp=[{(0,0):1}]
        for k in range(1,t+1):
            table=Counter()
            for (a,b),count in self.dp[-1].items():
                for (da,db,dc),multiplicity in histogram.items():table[a+da,b+db]+=count*multiplicity
            self.dp.append(dict(table))
        self.states=sorted(self.dp[-1]);self.capacity=sum(count**3 for count in self.dp[-1].values())
        self.histogram=histogram
    def unrank_sequence(self,profile,index):
        a,b=profile;count=self.dp[self.t][profile]
        if not 0<=index<count:raise IndexError('aggregate sequence')
        out=[]
        for remaining in range(self.t-1,-1,-1):
            for row,(da,db,dc) in enumerate(self.profiles):
                block=self.dp[remaining].get((a-da,b-db),0)
                if index>=block:index-=block;continue
                out.append(row);a-=da;b-=db;break
            else:raise AssertionError('aggregate unrank')
        if a or b or index:raise AssertionError('aggregate remainder')
        return tuple(out)
    def slot_rows(self,slot):
        if not 0<=slot<self.capacity:raise IndexError('aggregate slot')
        for profile in self.states:
            L=self.dp[self.t][profile];block=L**3
            if slot>=block:slot-=block;continue
            out=[]
            for orientation in range(3):
                slot,index=divmod(slot,L)
                out.extend((row,orientation) for row in self.unrank_sequence(profile,index))
            if slot:raise AssertionError('aggregate slot remainder')
            return tuple(out)
        raise AssertionError('aggregate selection')
    def binding(self):
        return dict(kind='equal_aggregate_profile_three_orientations',t=self.t,capacity=self.capacity,
            profile_histogram=[dict(profile=list(p),multiplicity=h) for p,h in sorted(self.histogram.items())],
            aggregate_profile_counts=[dict(first_two=list(p),count=self.dp[-1][p]) for p in self.states])


class FullBalancedSelection:
    """All balanced row tuples in fixed cyclic orientation order.

    Suffix DP counts the first two zero-mask sizes; the third follows from
    total width. Each orientation has t consecutive copies. No row tuple
    enumeration or tensor coefficient arrays are used.
    """
    def __init__(self,rows,t):
        if type(t)is not int or t<1:raise ValueError('positive t')
        self.t=t;self.width=len(rows[0]);self.target=self.width*t
        profiles=[tuple(row.count(c) for c in '123') for row in rows]
        self.profiles=[[(p[o],p[(o+1)%3]) for p in profiles] for o in range(3)]
        self.orientations=tuple(o for o in range(3) for _ in range(t))
        histograms=[Counter(p) for p in self.profiles]
        self.dp=[None]*(3*t+1);self.dp[-1]={(0,0):1}
        for pos in range(3*t-1,-1,-1):
            table=Counter()
            for (a,b),count in self.dp[pos+1].items():
                for (da,db),multiplicity in histograms[self.orientations[pos]].items():
                    na,nb=a+da,b+db
                    # The third mask must also fit the target. All increments
                    # are nonnegative, so this loses no balanced completion.
                    nc=(3*t-pos)*self.width-na-nb
                    if na<=self.target and nb<=self.target and nc<=self.target:
                        table[na,nb]+=count*multiplicity
            self.dp[pos]=dict(table)
        self.capacity=self.dp[0].get((self.target,self.target),0)
        if self.capacity<1:raise ValueError('no balanced tuples')
    def slot_rows(self,slot):
        if not 0<=slot<self.capacity:raise IndexError('balanced slot')
        a=b=self.target;out=[]
        for pos,orientation in enumerate(self.orientations):
            for row,(da,db) in enumerate(self.profiles[orientation]):
                block=self.dp[pos+1].get((a-da,b-db),0)
                if slot>=block:slot-=block;continue
                out.append((row,orientation));a-=da;b-=db;break
            else:raise AssertionError('balanced unrank')
        if a or b or slot:raise AssertionError('balanced unrank remainder')
        return tuple(out)
    def binding(self):
        return dict(kind='all_balanced_three_orientations',t=self.t,
            target_zero_count=self.target,capacity=self.capacity,
            orientation_order='t copies each of orientations 0, 1, 2',
            suffix_state_counts=[len(table) for table in self.dp])


class GradedPacket:
    """Direct sum of explicitly selected square products; direct-output W."""
    def __init__(self,rows,potentials,q,frequencies=None,family='profile_interpolation',*,t=None,allow_q1_control=False):
        if type(q)is not int or q<(1 if allow_q1_control else 2):raise ValueError('q>=2 required')
        self.rows=tuple(rows);self.potentials=tuple(tuple(x) for x in potentials);self.q=q
        self.replay=replay_grading(self.rows,self.potentials);self.degree=self.replay['degree'];self.width=len(rows[0])
        profiles=[tuple(row.count(ch) for ch in '123') for row in rows]
        self.aggregate=None
        if family=='full_balanced_interpolation':
            self.aggregate=FullBalancedSelection(rows,t);types=None;self.interpolate=True
        elif family=='aggregate_profile_interpolation':
            self.aggregate=AggregateProfileSelection(rows,t);types=None;self.interpolate=True
        elif family=='profile_interpolation':types=profiles;self.interpolate=True
        elif family=='joint_type_grade_zero':types=[profiles[i]+tuple(potentials[h][i] for h in range(3)) for i in range(len(rows))];self.interpolate=False
        elif family=='actual_row_grade_zero':types=list(range(len(rows)));self.interpolate=False
        else:raise ValueError('unknown selected tuple family')
        self.family=family
        self.sequences=None if self.aggregate is not None else FrequencySequences(types,frequencies)
        self.t=self.aggregate.t if self.aggregate is not None else self.sequences.t
        if self.t<1:raise ValueError('positive t required')
        self.copies=3*self.t;self.power=self.width*self.copies;self.axis=self.width*self.t
        self.n=q**self.axis;self.capacity=self.aggregate.capacity if self.aggregate is not None else self.sequences.count**3
        self.interpolation_degree=self.copies*self.degree if self.interpolate else 0
        self.evaluations=self.interpolation_degree+1;self.column_count=(q+1)**self.power
        self.rank=self.evaluations*self.column_count
    def binding(self):
        return dict(format='graded-cw-packet-frequency-v1',rows=list(self.rows),potentials=[list(x) for x in self.potentials],
            q=self.q,family=self.family,selection=self.aggregate.binding() if self.aggregate is not None else
                dict(kind='fixed_type_frequencies',frequencies=list(self.sequences.frequencies),types=self.sequences.types,sequence_count=self.sequences.count),t=self.t,power=self.power,inner_side=self.n,
            capacity=self.capacity,rank=self.rank,interpolation_factor=self.evaluations,grading_replay=self.replay)
    def slot_rows(self,slot):
        if not 0<=slot<self.capacity:raise IndexError('slot')
        if self.aggregate is not None:return self.aggregate.slot_rows(slot)
        L=self.sequences.count;out=[]
        for orientation in range(3):
            slot,index=divmod(slot,L)
            out.extend((row,orientation) for row in self.sequences.unrank(index))
        if slot:raise AssertionError('slot remainder')
        return tuple(out)
    def slot_geometry(self,slot):
        zeros=[];grade=[0,0,0]
        for row,orientation in self.slot_rows(slot):
            zeros.extend((int(ch)-1-orientation)%3 for ch in self.rows[row])
            for h in range(3):grade[h]+=self.potentials[(h+orientation)%3][row]
        if any(zeros.count(h)!=self.axis for h in range(3)):raise AssertionError('unbalanced selected tuple')
        if sum(grade):raise AssertionError('desired total grade')
        if not self.interpolate and any(grade):raise AssertionError('grade-zero subset failed')
        return tuple(zeros),tuple(grade)
    def word(self,mode,slot,row,col):
        if mode not in (0,1,2) or not(0<=row<self.n and 0<=col<self.n):raise IndexError('variable')
        zeros,grades=self.slot_geometry(slot)
        row_role,col_role=((1,2),(2,0),(1,0))[mode]
        out=[]
        for z in zeros:
            if z==mode:out.append(0)
            elif z==row_role:row,digit=divmod(row,self.q);out.append(digit+1)
            elif z==col_role:col,digit=divmod(col,self.q);out.append(digit+1)
            else:raise AssertionError('role assignment')
        if row or col:raise AssertionError('base-q remainder')
        return tuple(out),grades[mode]
    def coefficient(self,mode,term,slot,row,col):
        if not 0<=term<self.rank:raise IndexError('term')
        evaluation,column=divmod(term,self.column_count);node=evaluation+1
        word,grade=self.word(mode,slot,row,col)
        coefficient=Fraction(node)**grade if self.interpolate else Fraction(1)
        if mode==2 and self.interpolate:
            coefficient*=(-1)**evaluation*comb(self.evaluations,node)
        for value in word:
            column,label=divmod(column,self.q+1)
            if label==self.q:
                factor=Fraction(1) if value==0 else Fraction(3,2*self.q)
                if mode==2:factor*=Fraction(-4*self.q,3)
            else:
                factor=Fraction(1) if value==0 else Fraction(label==value-1)+Fraction(1,self.q)
            coefficient*=factor
        if column:raise AssertionError('column digit remainder')
        return coefficient


class UnitOuter:
    n=1;rank=1
    def coefficient(self,mode,term,row,col):
        if (mode not in (0,1,2)) or(term,row,col)!=(0,0,0):raise IndexError
        return Fraction(1)
    def binding(self):return dict(kind='unit',side=1,rank=1)


class StrassenOuter:
    n=2;rank=7
    U=((1,0,0,1),(0,0,1,1),(1,0,0,0),(0,0,0,1),(1,1,0,0),(-1,0,1,0),(0,1,0,-1))
    V=((1,0,0,1),(1,0,0,0),(0,1,0,-1),(-1,0,1,0),(0,0,0,1),(1,1,0,0),(0,0,1,1))
    W=((1,0,0,1),(0,0,1,-1),(0,1,0,1),(1,0,1,0),(-1,1,0,0),(0,0,0,1),(1,0,0,0))
    def coefficient(self,mode,term,row,col):
        if mode not in (0,1,2) or not 0<=term<7 or not(0<=row<2 and 0<=col<2):raise IndexError
        return Fraction((self.U,self.V,self.W)[mode][term][2*row+col])
    def binding(self):return dict(kind='explicit_strassen_7',side=2,rank=7)


class CatalogueOuter:
    """Generic pinned-catalogue adapter, with exact n/rank/root validation."""
    def __init__(self,scheme,expected_n,expected_rank,expected_root_object):
        if(scheme.n,scheme.rank)!=(expected_n,expected_rank):raise ValueError('catalogue dimensions/rank mismatch')
        if scheme.record['root_object']!=expected_root_object:raise ValueError('catalogue root mismatch')
        self.scheme=scheme;self.n=expected_n;self.rank=expected_rank;self.root=expected_root_object
    def coefficient(self,mode,term,row,col):return Fraction(self.scheme.coefficient(mode,term,row,col))
    def binding(self):return dict(kind='catalogue_outer',side=self.n,rank=self.rank,root_object=self.root)


class CatalogueLita32:
    """Adapter to the existing verified base; only the small base is loaded.

    The arithmetic core remains Fraction based. The existing catalogue loader
    may internally use its own small arrays/dependencies. Callers must supply
    and hash-bind the intended catalogue root; this is not an unpinned download.
    """
    n=32;rank=14196
    def __init__(self,scheme,expected_root_object):
        if(scheme.n,scheme.rank)!=(32,14196):raise ValueError('wrong LITA base')
        if scheme.record['root_object']!=expected_root_object:raise ValueError('LITA root mismatch')
        self.scheme=scheme;self.root=expected_root_object
    def coefficient(self,mode,term,row,col):return Fraction(self.scheme.coefficient(mode,term,row,col))
    def binding(self):return dict(kind='catalogue_lita32',side=32,rank=14196,root_object=self.root)


class TensorOuter:
    def __init__(self,factors):
        self.factors=tuple(factors);self.n=prod(f.n for f in factors);self.rank=prod(f.rank for f in factors)
    def coefficient(self,mode,term,row,col):
        if mode not in(0,1,2)or not 0<=term<self.rank or not(0<=row<self.n and 0<=col<self.n):raise IndexError
        value=Fraction(1)
        for factor in self.factors:
            term,t=divmod(term,factor.rank);row,r=divmod(row,factor.n);col,c=divmod(col,factor.n)
            value*=factor.coefficient(mode,t,r,c)
            if not value:return value
        return value
    def binding(self):return dict(kind='tensor_product',side=self.n,rank=self.rank,factors=[f.binding() for f in self.factors])


class GradedSquare:
    def __init__(self,packet,outer,*,crop=None):
        self.packet=packet;self.outer=outer;self.natural_n=packet.n*outer.n
        self.n=self.natural_n if crop is None else crop
        if type(self.n)is not int or not 1<=self.n<=self.natural_n:raise ValueError('invalid crop')
        self.batches=(outer.rank+packet.capacity-1)//packet.capacity;self.rank=self.batches*packet.rank
    def coefficient_summands(self,mode,term,row,col):
        if mode not in(0,1,2)or not 0<=term<self.rank or not(0<=row<self.n and 0<=col<self.n):raise IndexError
        batch,local_term=divmod(term,self.packet.rank)
        block_row,inner_row=divmod(row,self.packet.n);block_col,inner_col=divmod(col,self.packet.n)
        start=batch*self.packet.capacity;stop=min(start+self.packet.capacity,self.outer.rank)
        for call in range(start,stop):
            scalar=self.outer.coefficient(mode,call,block_row,block_col)
            if scalar:yield scalar*self.packet.coefficient(mode,local_term,call-start,inner_row,inner_col)
    def coefficient(self,mode,term,row,col,*,max_slot_terms=10000):
        batch=term//self.packet.rank
        calls=max(0,min(self.packet.capacity,self.outer.rank-batch*self.packet.capacity))
        if max_slot_terms is not None and calls>max_slot_terms:
            raise OverflowError('This exact coefficient requires an enormous finite slot sum; pass max_slot_terms=None to request it explicitly.')
        return sum(self.coefficient_summands(mode,term,row,col),Fraction(0))
    def binding(self):
        return dict(format='graded-cw-square-frequency-v1',packet=self.packet.binding(),outer=self.outer.binding(),
            natural_n=self.natural_n,n=self.n,rank=self.rank,batches=self.batches,
            sub265_exact=self.rank**20<self.n**53,
            coefficient_access='Exact finite slot sum; generally computationally prohibitive at target size; not a fast random-access claim')
