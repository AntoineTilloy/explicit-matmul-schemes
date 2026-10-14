"""Package adapter for exact finite graded CW square constructions.

The coefficient core is .graded_packet_core. Existing
proof objects remain immutable dependencies. Base engines load only when a
coefficient is actually requested; audit reads hash-checked proof metadata.
"""
from fractions import Fraction
import json
from pathlib import Path
import re
from .. import _resources
from .graded_packet_core import GradedPacket,GradedSquare,StrassenOuter,TensorOuter

FORMAT='graded-cw-square-v1'


def _positive_int(value,label):
    if isinstance(value,bool):raise ValueError(label+' must be an integer')
    if isinstance(value,str):
        if not value.isdecimal():raise ValueError(label+' must be positive decimal')
        value=int(value)
    if type(value)is not int or value<1:raise ValueError(label+' must be positive')
    return value


def _nonnegative_int(value,label):
    if type(value)is not int or value<0:raise ValueError(label+' must be nonnegative integer')
    return value


def _proof_digest(proof):
    if not isinstance(proof,str) or not re.fullmatch(r'proofs/[0-9a-f]{64}\.json',proof):
        raise ValueError('Base proof must be proofs/<sha256>.json')
    digest=proof[len('proofs/'):-len('.json')]
    if _resources.index()['paths'].get(proof)!=digest:raise ValueError('Base proof path/hash mismatch')
    return digest


def _load_base_engine(data):
    kind=data.get('format')
    if kind==FORMAT:return load_certificate(data)
    if kind in ('cw-single-square-v1','ordinary-square-product-v1'):
        from .legacy_square import load_certificate as loader
    elif kind=='explicit-qcsr-square-v1':
        from .explicit_square import load_certificate as loader
    else:
        from .general_square import load_certificate as loader
    return loader(data)


class DeferredProofOuter:
    def __init__(self,proof,n,rank):
        self.proof=proof;self.digest=_proof_digest(proof)
        self.n=_positive_int(n,'base side');self.rank=_positive_int(rank,'base rank')
        self._engine=None;self._metadata=None
    def audit(self):
        if self._metadata is None:
            # raw_object verifies size and SHA without instantiating an engine.
            data=json.loads(_resources.raw_object(self.digest))
            if(_positive_int(data['n'],'proof side'),_positive_int(data['rank'],'proof rank'))!=(self.n,self.rank):
                raise ValueError('Referenced base proof counts differ from certificate')
            if data.get('coefficient_field')!='Q':
                raise ValueError('Referenced base proof must declare rational coefficients')
            self._metadata=data
        return dict(proof=self.proof,root_object=self.digest,n=self.n,rank=self.rank,
                    proof_format=self._metadata['format'],integrity='PASS hash-checked resource',
                    engine_loaded=self._engine is not None)
    def coefficient(self,mode,term,row,col):
        if self._engine is None:
            self.audit()
            path=_resources.prepare(self.digest)
            self._engine=_load_base_engine(json.loads(path.read_text()))
            if(self._engine.n,self._engine.rank)!=(self.n,self.rank):raise ValueError('Loaded base counts mismatch')
        return Fraction(self._engine.factor_entry(mode,term,row,col))
    def binding(self):return dict(kind='deferred_proof_outer',proof=self.proof,side=self.n,rank=self.rank)


class GradedPackageEngine:
    finite_slot_sum=True
    def __init__(self,data):
        if data.get('format')!=FORMAT:raise ValueError('Unsupported graded certificate')
        if data.get('coefficient_field')!='Q':raise ValueError('Rational coefficient field required')
        self.data=data;p=data['packet']
        self.packet=GradedPacket(p['rows'],p['potentials'],_positive_int(p['q'],'q'),
            p.get('frequencies'),family=p['family'],t=p.get('t'))
        expected=p['expected']
        if(self.packet.n,self.packet.rank,self.packet.capacity)!=(
                _positive_int(expected['inner_side'],'packet side'),
                _positive_int(expected['rank'],'packet rank'),
                _positive_int(expected['capacity'],'packet capacity')):
            raise ValueError('Packet counts do not replay')
        if 'power' in expected and self.packet.power!=expected['power']:raise ValueError('Packet power mismatch')
        if 'interpolation_factor' in expected and self.packet.evaluations!=expected['interpolation_factor']:
            raise ValueError('Interpolation factor mismatch')
        outer=data['outer'];factors=[];self.base_proofs=[];memo={}
        for entry in outer.get('bases',[]):
            copies=_positive_int(entry['copies'],'base copies')
            key=(entry['proof'],_positive_int(entry['n'],'base side'),_positive_int(entry['rank'],'base rank'))
            if key not in memo:
                memo[key]=DeferredProofOuter(*key);self.base_proofs.append(memo[key])
                memo[key].audit()
            factors.extend([memo[key]]*copies)
        strassen_copies=_nonnegative_int(outer.get('strassen_copies',0),'Strassen copies')
        factors.extend([StrassenOuter()]*strassen_copies)
        self.outer=TensorOuter(factors)
        if(self.outer.n,self.outer.rank)!=(
                _positive_int(outer['expected_side'],'outer side'),
                _positive_int(outer['expected_rank'],'outer rank')):
            raise ValueError('Outer product counts do not replay')
        crop=None if data.get('crop') is None else _positive_int(data['crop'],'crop')
        self.square=GradedSquare(self.packet,self.outer,crop=crop)
        self.n,self.rank=self.square.n,self.square.rank;self.field='Q'
        if(self.n,self.rank)!=(
                _positive_int(data['n'],'square side'),_positive_int(data['rank'],'square rank')):
            raise ValueError('Square counts do not replay')
        if self.square.natural_n!=_positive_int(data['natural_n'],'natural side'):
            raise ValueError('Natural square size mismatch')
        claim=data.get('exponent_claim')
        if claim is not None:
            if claim!={'numerator':53,'denominator':20,'strict':True}:
                raise ValueError('Only the explicit strict sub-2.65 claim is supported here')
            if not self.rank**20<self.n**53:raise ValueError('Claimed exponent threshold fails')
    def factor_entry(self,mode,term,row,col):
        """Standard engine API; refuses >10000 active slot terms by default.

        The finite definition exists for all valid indices, but evaluation at
        large catalogue dimensions can be computationally prohibitive. Use
        coefficient(..., max_slot_terms=None) only to explicitly request the
        full sum, or inspect coefficient_summands lazily.
        """
        return self.square.coefficient(mode,term,row,col,max_slot_terms=10000)
    def coefficient(self,mode,term,row,col,*,max_slot_terms=10000):
        return self.square.coefficient(mode,term,row,col,max_slot_terms=max_slot_terms)
    def coefficient_summands(self,mode,term,row,col):
        return self.square.coefficient_summands(mode,term,row,col)
    def audit(self):
        references=[base.audit() for base in self.base_proofs]
        return dict(format='graded-cw-square-audit-v1',n=self.n,rank=self.rank,
            natural_n=self.square.natural_n,packet=self.packet.binding(),
            outer=self.outer.binding(),base_references=references,
            batches=self.square.batches,strict_sub265=self.rank**20<self.n**53,
            coefficient_field='Q',block_recursive=True,
            limits='Exact construction replay and referenced base metadata; no expansion of the large square tensor; coefficient evaluation may require an enormous finite slot sum')


def load_certificate(path_or_data):
    data=path_or_data if isinstance(path_or_data,dict) else json.loads(Path(path_or_data).read_text())
    return GradedPackageEngine(data)
