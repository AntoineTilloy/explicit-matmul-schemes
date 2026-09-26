"""Integrity checks and exact replay of the supplied finite reduction witnesses.

The proof basis is the documented CW identities and published primitive schemes;
this is an arithmetic checker, not a formal proof-assistant implementation.
"""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import time
from . import _resources, load

_AUDITED=set()
_EXACT=set()


def read(path):return json.loads(Path(path).read_text())


def proof_path(s):
    p=Path(s)
    if not p.is_absolute():p=_resources.root()/p
    return p if p.suffix else p.with_suffix('.json')


def audit_supports(value):
    from ._engine.restrictions import audit,verify_completion
    checked=0
    if isinstance(value,dict):
        if 'rows' in value and 'm' in value and isinstance(value['rows'],list) and value['rows']:
            rows=value['rows'];m=value['m']
            if isinstance(rows[0],list) and len(rows[0])==3 and all(type(x) is int for x in rows[0]):
                key=(m,tuple(map(tuple,rows)))
                if key not in _AUDITED:audit(rows,m);_AUDITED.add(key);checked+=1
        if 'q' in value and type(value['q']) is int:
            key=('completion',value['q'])
            if key not in _AUDITED:verify_completion(value['q']);_AUDITED.add(key)
        for v in value.values():checked+=audit_supports(v)
    elif isinstance(value,list):
        for v in value:checked+=audit_supports(v)
    return checked


def replay_orbit(d,threads):
    import numpy as np
    from ._engine.rational_state import ExactState
    from ._engine.rational_refinement import ExactMultiSplitState
    from ._engine.rational_arrays import qm,from_qm,Q
    from ._engine.tailored_completion import make_state
    if d['format']=='cw-refined-orbit-relations-v1':
        state=replay_orbit(read(proof_path(d['base_proof'])),threads) if 'base_proof' in d else make_state(d['q'],d['rows'],d['m'],exact=True)
        for coordinate in d['coordinates']:state=ExactMultiSplitState.split(state,coordinate)
    else:state=ExactState(d['q'],d['rows'],d['m'])
    for step in d['steps']:
        pair=step['pair'];assert len(pair)==2 and len(set(pair))==2 and set(pair)<={0,1,2}
        zeros=set(step.get('zero_patterns',[]));assert zeros<=set(state.zeros())
        keep=set(step['kept']);deleted=set(step['deleted'])
        assert not(keep&deleted or keep&zeros or deleted&zeros)
        assert keep|deleted|zeros==set(state.kept)
        if zeros:
            state.eliminate_exact([],dict(kept=sorted(set(state.kept)-zeros),deleted=[],bases={}),pair)
        G=state.pair_gram(pair) if deleted else []
        new=[[] for _ in range(3)];free=3-sum(pair);used=set()
        for I,patterns in enumerate(state.patterns):
            K=[j for j,H in enumerate(patterns) if H in keep]
            J=[j for j,H in enumerate(patterns) if H in deleted]
            for mode in range(3):new[mode].append(state.coefs[mode][I][:,K].copy())
            if not J:continue
            b=step['blocks'][str(I)];used.add(str(I))
            assert b['deleted']==[patterns[j] for j in J]
            assert set(b['basis'])<=keep and len(set(b['basis']))==len(b['basis'])
            B=[state.lookups[I][H] for H in b['basis']]
            C=np.array([[Q(x) for x in row] for row in b['coefficients']],dtype=object).reshape(len(B),len(J))
            assert qm(G[I][:,B])*qm(C)==qm(G[I][:,J]),('false rational relation',I)
            pos={old:j for j,old in enumerate(K)}
            dest=[pos[j] for j in B]
            new[free][I][:,dest]+=from_qm(qm(state.coefs[free][I][:,J])*qm(C).transpose())
        assert used==set(step['blocks'])
        state.kept=step['kept'];state.coefs=new;state.reindex()
        assert state.rank()==step['rank']
    assert state.rank()==d['rank']
    return state


def replay_quotient(path,d,threads):
    import numpy as np
    from flint import fmpq_mat
    from ._engine.quotient_rational_replay import reconstruct
    from ._engine.quotient_successive_replay import residue,coefficient_bound
    from ._engine.quotient_invariant_gram import representation
    if d.get('certificate_method')=='exact-single-factor-zero-norm':
        for h in d['history']:exact_proof(proof_path(h),threads)
        from ._engine.quotient_zero_factor_replay import verify_frozen
        verify_frozen(path);return
    metadata=read(proof_path(d['source']));q=metadata['q'];pair=tuple(d['pair'])
    _,_,orbits,blocks=representation(q,metadata['rows'],metadata['m'],metadata['groups'],math.prod((-1,-1,1)[i] for i in pair))
    assert json.loads(json.dumps(blocks))==metadata['blocks'],'Incorrect representation metadata'
    removed=set(d['deleted_orbits'])
    assert len(removed)==len(d['deleted_orbits'])
    assert sum(orbits[i]['size'] for i in removed)==d['savings']
    ids=[i for i,b in enumerate(blocks) if any(r['orbit'] in removed for r in b['records'])]
    assert set(map(str,ids))==set(d['relations'])
    history=[proof_path(h) for h in d.get('history',[])]
    for h in history:exact_proof(h,threads)
    if not history:
        matrices,_,_,_=reconstruct(metadata,ids,threads)
    else:
        from sympy import nextprime
        denominators=[1,1,1];magnitudes=[2,2,2]
        for h in history:
            prev=read(h);den,mag=coefficient_bound(h);free=3-sum(prev['pair'])
            denominators[free]*=den;magnitudes[free]*=1+prev['savings']*mag
        denominator=(denominators[pair[0]]*denominators[pair[1]])**2
        bound=denominator*4*len(metadata['rows'])**2*q**(10*metadata['m'])*(magnitudes[pair[0]]*magnitudes[pair[1]])**2
        modulus=1;prime=100000000;matrices=None
        while modulus<=2*bound:
            prime=int(nextprime(prime));current=residue(q,prime,history,pair,ids,denominator,threads)
            if matrices is None:matrices=[g.astype(object) for g in current]
            else:
                inv=pow(modulus,-1,prime)
                for g,r in zip(matrices,current):g+=modulus*(((r-g)%prime)*inv%prime)
            modulus*=prime
        for g in matrices:
            g[g>modulus//2]-=modulus
            assert all(abs(int(v))<=bound for v in g.flat)
    for i,G in zip(ids,matrices):
        relation=d['relations'][str(i)];B=relation['basis'];J=relation['deleted']
        assert J==[j for j,r in enumerate(blocks[i]['records']) if r['orbit'] in removed]
        assert not set(B)&set(J) and len(B)==len(set(B))
        C=fmpq_mat(relation['coefficients'])
        assert fmpq_mat(G[:,B].tolist())*C==fmpq_mat(G[:,J].tolist()),('false quotient relation',i)


def verify_z4(engine,exact,threads):
    """Z4-quotient units: structure and count audit, exact packet identity at sampled entries, orbit lemma;
    with exact=True the bundled q-independent relation certificate is replayed exactly over Q at each unit's q."""
    from ._engine import z4_checks
    from ._engine.z4_quotient import units,replay_document
    distinct={}
    for u in units(engine):distinct.setdefault((u.q,u.support_path,u.relations_path),u)
    out=[]
    for (q,_,relations),u in sorted(distinct.items()):
        rec=dict(q=q,unit_rank=str(u.rank),audit=z4_checks.audit(u))
        packet=z4_checks.packet_identity(u,samples=4,seed=q);packet.pop('records')
        assert packet['passed'],'Z4 packet identity fails'
        rec['exact_packet_entries']=packet
        rec['orbit_lemma']=z4_checks.orbit_lemma(u,samples=2,seed=q)
        if exact:
            start=time.monotonic();state=replay_orbit(replay_document(relations,q),threads)
            assert sorted(state.kept)==sorted(u.packet.kept) and state.rank()==u.packet.rank
            rec['relation_replay']=dict(exact=True,unquotiented_rank=str(state.rank()),seconds=time.monotonic()-start)
        out.append(rec)
    return out


def exact_proof(path,threads=1):
    path=Path(path);digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest in _EXACT:return
    d=read(path);fmt=d.get('format','')
    if fmt in ('cw-compact-orbit-relations-v1','cw-refined-orbit-relations-v1'):
        replay_orbit(d,threads)
    elif fmt in ('cw-quotient-original-orbit-rational-v1','cw-quotient-successive-stage-rational-v1'):
        replay_quotient(path,d,threads)
    elif fmt=='cw-quotient-successive-rational-v1':
        for h in d['proofs']:exact_proof(proof_path(h),threads)
    elif 'relations' in fmt or 'rational' in fmt:
        raise NotImplementedError(f'Unsupported exact proof format: {fmt}')
    _EXACT.add(digest)


def verify(n,*,exact=False,threads=1):
    if not __debug__:raise RuntimeError('Verification requires Python without -O')
    if not 1<=threads<=12:raise ValueError('threads must be between 1 and 12')
    start=time.monotonic();s=load(n)
    path=_resources.root()/f"proofs/{s.record['root_object']}.json"
    visited=set();supports=0;proofs=0
    def visit(p):
        nonlocal supports,proofs
        if p in visited:return
        visited.add(p);d=read(p);supports+=audit_supports(d)
        for ref in set(_resources.references(d)):visit(_resources.root()/ref)
        fmt=d.get('format','')
        if 'relations' in fmt or 'rational' in fmt:
            proofs+=1
            if exact:exact_proof(p,threads)
    visit(path)
    explicit=None
    if read(path).get('format')=='explicit-qcsr-square-v1':
        explicit=dict(random_evaluation_mod_2_31_minus_1=s._engine.random_check(trials=2))
        assert explicit['random_evaluation_mod_2_31_minus_1'],'Explicit scheme fails random evaluation'
        if exact:
            explicit['all_coefficients']=s._engine.exact_check()
            assert explicit['all_coefficients']['exact'],'Explicit scheme fails the exact coefficient check'
    z4=None
    from ._engine.z4_quotient import units as z4_units
    if z4_units(s._engine):
        z4=verify_z4(s._engine,exact,threads)
    queries=[]
    for factor in ('U','V','W'):
        for term,i,j in [(0,s.n-1,s.n-1),(s.rank-1,s.n-1,s.n-1)]:
            queries.append(str(s.coefficient(factor,term,i,j)))
    return dict(n=str(s.n),rank=str(s.rank),mode='exact reduction replay' if exact else 'integrity, support, count and coefficient checks',
                dependencies=len(visited),new_support_audits=supports,reduction_proofs=proofs,coefficient_queries=queries,
                sub_2_7=s.rank**10<s.n**27,**({'explicit_checks':explicit} if explicit else {}),**({'z4_quotient_checks':z4} if z4 else {}),seconds=time.monotonic()-start,
                proof_basis='Documented base identities and published primitive schemes; see docs/verification.md')
