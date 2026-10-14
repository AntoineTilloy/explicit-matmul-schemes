#!/usr/bin/env python3
"""SOURCE ONLY until review: register exact graded rational square certificates.

Run on Atlas. All certificate/count checks finish before any repository write.
No large coefficient queries or new scientific search occur in this exporter.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import resource
import sys

LITA32_ROOT='a7fb19d69d8e078fb0d6d614ade10417225a186f4fbf94a8d4292f3935995e05'


def sha(raw):return hashlib.sha256(raw).hexdigest()
def canonical(data):return json.dumps(data,sort_keys=True,separators=(',',':')).encode()
def pretty(data):return (json.dumps(data,indent=2)+'\n').encode()
def read_bound(path,expected):
    raw=path.read_bytes()
    if sha(raw)!=expected:raise ValueError('Input hash mismatch: '+str(path))
    return json.loads(raw),raw


def integer(value,minimum=1):
    if isinstance(value,bool):raise ValueError('Boolean count')
    if isinstance(value,str) and value.isdecimal():value=int(value)
    if type(value)is not int or value<minimum:raise ValueError('Invalid integer count')
    return value


def lita32_comparison(n,byroot):
    """A specific valid tensor-power/padding baseline, not a best bound."""
    base=byroot[LITA32_ROOT]
    if (integer(base['n']),integer(base['rank']),base['field'])!=(32,14196,'Q'):
        raise ValueError('Pinned rational size32 comparison base mismatch')
    leaf=dict(kind='catalogue',n=32,rank=14196,root_object=LITA32_ROOT)
    recipe=leaf
    while recipe['n']<n:
        recipe=dict(kind='product',n=32*recipe['n'],rank=14196*recipe['rank'],left=leaf,right=recipe)
    if recipe['n']!=n:recipe=dict(kind='pad',n=n,rank=recipe['rank'],child=recipe)
    def replay(node):
        if node['kind']=='catalogue':
            if node!=leaf:raise ValueError('Unexpected comparison primitive')
            return 32,14196
        if node['kind']=='product':
            ln,lr=replay(node['left']);rn,rr=replay(node['right'])
            if (node['n'],node['rank'])!=(ln*rn,lr*rr):raise ValueError('Comparison product mismatch')
        elif node['kind']=='pad':
            cn,cr=replay(node['child'])
            if not 1<=node['n']<=cn or node['rank']!=cr:raise ValueError('Comparison padding mismatch')
        else:raise ValueError('Comparison kind')
        return node['n'],node['rank']
    if replay(recipe)!=(n,recipe['rank']):raise ValueError('Comparison replay failed')
    return recipe


def outer_spec(item,byroot):
    old=item['outer'];bases=[]
    def add(root,copies,n=None,rank=None):
        copies=integer(copies,0)
        if not copies:return
        record=byroot[root]
        if record['field']!='Q':raise ValueError('Non-rational outer base')
        bn,br=integer(record['n']),integer(record['rank'])
        if n is not None and integer(n)!=bn:raise ValueError('Base side mismatch')
        if rank is not None and integer(rank)!=br:raise ValueError('Base rank mismatch')
        bases.append(dict(proof='proofs/'+root+'.json',n=str(bn),rank=str(br),copies=copies))
    if 'bases' in old:
        for entry in old['bases']:
            root=entry.get('root_object')
            if root is None:
                proof=entry['proof']
                if not proof.startswith('proofs/') or not proof.endswith('.json'):raise ValueError('Invalid proof path')
                root=proof[7:-5]
            add(root,entry['copies'],entry.get('n'),entry.get('rank'))
    else:
        if old.get('public_base_copies',0):
            add(old['public_root'],old['public_base_copies'],old.get('public_base_n'),old.get('public_base_rank'))
        add(LITA32_ROOT,old.get('lita32_copies',0))
    j=integer(old.get('strassen2_copies',old.get('strassen_copies',0)),0)
    return dict(bases=bases,strassen_copies=j,
        expected_side=str(integer(old['side'])),expected_rank=str(integer(old['rank'])))


def certificate(item,label,grading,byroot,bindings):
    puzzle=item['puzzle'];target=next(x for x in grading['targets'] if x['puzzle']==puzzle)
    if target['status']!='CERTIFIED_POSITIVE_GRADING':raise ValueError('Uncertified grading')
    t,q=integer(item['t']),integer(item['q'],2)
    power=3*puzzle[1]*t;degree=integer(item['degree'],0)
    if integer(item['power'])!=power or integer(item['base_copies'])!=3*t:raise ValueError('Power mismatch')
    if target['certificate']['degree']!=degree:raise ValueError('Grading degree mismatch')
    evaluations=3*t*degree+1
    if integer(item['interpolation_factor'])!=evaluations:raise ValueError('Interpolation price mismatch')
    packet_rank=evaluations*(q+1)**power;capacity=integer(item['square_capacity'])
    outer=outer_spec(item,byroot);calls=int(outer['expected_rank'])
    batches=(calls+capacity-1)//capacity
    if 'batches' in item and integer(item['batches'])!=batches:raise ValueError('Batch mismatch')
    rank=integer(item['rank'])
    if rank!=batches*packet_rank:raise ValueError('Total rank mismatch')
    natural=integer(item['natural_n']);crop=integer(item['padded_n'])
    if not (crop<=natural and (crop-1)**53<=rank**20<crop**53):raise ValueError('Strict crop threshold mismatch')
    n=natural if label=='best_natural' else crop
    family=item.get('family')
    if family is None:
        if not item['construction'].startswith(('All balanced row tuples','Full balanced cyclic packet')):raise ValueError('Unrecognized tuple selection')
        family='full_balanced_interpolation'
    packet=dict(rows=target['rows'],potentials=target['certificate']['potentials'],q=q,t=t,family=family,
        expected=dict(inner_side=str(integer(item['inner_side'])),rank=str(packet_rank),capacity=str(capacity),
            power=power,interpolation_factor=evaluations))
    if 'frequencies' in item:packet['frequencies']=item['frequencies']
    return dict(format='graded-cw-square-v1',coefficient_field='Q',n=str(n),rank=str(rank),natural_n=str(natural),
        crop=None if n==natural else str(n),packet=packet,outer=outer,
        exponent_claim=dict(numerator=53,denominator=20,strict=True),
        provenance=dict(selection=label,input_hashes=bindings,construction=item['construction'],
            scope='Finite rational bilinear construction; bounded search, no minimum-size claim'))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--repo',type=Path,required=True)
    for name in ('candidate','grading','audit','controls'):
        p.add_argument('--'+name,type=Path,required=True)
        p.add_argument('--'+name+'-sha',required=True)
    for name in ('core','adapter','catalogue','index','comparison'):
        p.add_argument('--'+name+'-sha',required=True)
    a=p.parse_args();sys.set_int_max_str_digits(0)
    resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(512<<20,512<<20))
    repo=a.repo.resolve();data_dir=repo/'src/cwschemes/data'
    core=repo/'src/cwschemes/_engine/graded_packet_core.py'
    adapter=repo/'src/cwschemes/_engine/graded_square.py'
    for path,expected in ((core,a.core_sha),(adapter,a.adapter_sha)):
        if sha(path.read_bytes())!=expected:raise ValueError('Package source hash mismatch: '+str(path))
    catalogue,cat_raw=read_bound(repo/'catalogue.json',a.catalogue_sha)
    if (data_dir/'catalogue.json').read_bytes()!=cat_raw:raise ValueError('Catalogue mirrors differ')
    index,index_raw=read_bound(data_dir/'index.json',a.index_sha)
    comparison_path=repo/'comparisons/large_recipes.json'
    comparisons,_=read_bound(comparison_path,a.comparison_sha)
    if len(catalogue)!=97:raise ValueError('This reviewed export expects 97 prior catalogue entries')
    inputs={};raw_inputs={};bindings={}
    for name in ('candidate','grading','audit','controls'):
        expected=getattr(a,name+'_sha');bindings[name]=expected
        inputs[name],raw_inputs[name]=read_bound(getattr(a,name),expected)
    candidate,grading,audit,controls=(inputs[k] for k in ('candidate','grading','audit','controls'))
    if candidate['grading_sha256']!=a.grading_sha or candidate['catalogue_sha256']!=a.catalogue_sha:
        raise ValueError('Candidate provenance mismatch')
    if audit['status']!='PASS':raise ValueError('Independent audit not PASS')
    for key,expected in (('result',a.candidate_sha),('grading',a.grading_sha),('catalogue',a.catalogue_sha)):
        if audit['input_hashes'][key]!=expected:raise ValueError('Independent audit binding mismatch')
    audited={x['item'] for x in audit['items'] if all(x.get(k)=='PASS' for k in ('exact_capacity','exact_price','exact_crop'))}
    if not {'best_natural','best_padded'}<=audited:raise ValueError('Both exported candidates must be independently audited')
    if controls['status']!='COMPLETE' or controls['generator_sha256']!=a.core_sha:raise ValueError('Controls/core mismatch')
    if len(controls['controls'])<10 or not all(x['status'].startswith('PASS') for x in controls['controls']):
        raise ValueError('Incomplete exact controls')
    byroot={x['root_object']:x for x in catalogue};by_n={int(x['n']):x for x in catalogue}
    sys.path.insert(0,str(repo/'src'))
    from cwschemes._engine.graded_square import load_certificate
    from cwschemes import _resources
    if _resources.DATA.resolve()!=data_dir.resolve():raise ValueError('Wrong package import')
    staged={};records=[];replays=[]
    def stage(path,raw):
        path=Path(path)
        if path in staged and staged[path]!=raw:raise ValueError('Conflicting staged output')
        staged[path]=raw
    for label in ('best_natural','best_padded'):
        cert=certificate(candidate[label],label,grading,byroot,bindings)
        engine=load_certificate(cert);replay=engine.audit()
        if not replay['strict_sub265'] or any(ref['engine_loaded'] for ref in replay['base_references']):
            raise ValueError('Audit must pass without loading a base coefficient engine')
        raw=canonical(cert);root=sha(raw);n,rank=int(cert['n']),int(cert['rank'])
        if n in by_n:raise ValueError('Refusing to replace an existing catalogue dimension')
        if any(int(x['n'])==n for x in records):raise ValueError('Two export labels have the same dimension')
        packed=gzip.compress(raw,compresslevel=9,mtime=0)
        if sha(gzip.decompress(packed))!=root:raise AssertionError('Compressed root mismatch')
        if root in index['objects']:raise ValueError('Root already registered')
        index['objects'][root]=dict(bytes=len(raw),compressed_bytes=len(packed))
        index['paths']['proofs/'+root+'.json']=root
        stage(data_dir/'objects'/f'{root}.gz',packed)
        record=dict(format='cwschemes-catalogue-entry-v1',n=str(n),rank=str(rank),field='Q',root_object=root,
            effective_exponent=math.log(rank)/math.log(n),
            provenance=dict(construction='Explicit finite positively graded CW tensor power with all balanced row tuples, rational interpolation, and shared outer-call batches',selection=label),
            verification='Exact grading/capacity/rank/crop replay and small end-to-end rational tensor controls in checks/graded. Large coefficients are explicit finite slot sums, generally computationally prohibitive; the large square tensor was not expanded.')
        records.append(record);stage(repo/'schemes'/f'{n}.json',pretty(record))
        if str(n) in comparisons:raise ValueError('Comparison size already registered')
        comparisons[str(n)]=lita32_comparison(n,byroot)
        replays.append(dict(label=label,root_object=root,audit=replay))
    catalogue.extend(records);catalogue.sort(key=lambda x:int(x['n']))
    stage(repo/'catalogue.json',pretty(catalogue));stage(data_dir/'catalogue.json',pretty(catalogue))
    stage(data_dir/'index.json',pretty(index))
    stage(comparison_path,pretty(comparisons))
    receipt_manifest={}
    receipt_names={'source.py','helper.py','full.py','graded_lazy_scheme.py','run.sh','controller.pid',
        'controller.log','resources.txt','terminal.json','stdout.log','stderr.log','launch_hashes.txt','bindings.sha256','runtime.txt',
        'catalogue.json','grading.json','cached.json','al_puzzles.py'}
    for name in ('candidate','grading','audit','controls'):
        receipt_dir=repo/'checks/graded'/name
        stage(receipt_dir/'result.json',raw_inputs[name])
        receipt_manifest[str((receipt_dir/'result.json').relative_to(repo))]=sha(raw_inputs[name])
        parent=getattr(a,name).parent
        for child in sorted(parent.iterdir()):
            if child.name not in receipt_names or not child.is_file():continue
            raw=child.read_bytes()
            source_sha=inputs[name].get('source_sha256',inputs[name].get('bindings',{}).get('source_sha256'))
            if child.name=='source.py' and source_sha!=sha(raw):
                raise ValueError('Receipt source mismatch: '+str(child))
            output=receipt_dir/child.name;stage(output,raw)
            receipt_manifest[str(output.relative_to(repo))]=sha(raw)
    stage(repo/'checks/graded/export_graded.py',Path(__file__).read_bytes())
    stage(repo/'checks/graded/export.json',pretty(dict(status='PASS',source_sha256=sha(Path(__file__).read_bytes()),
        input_hashes=bindings,package_core_sha256=a.core_sha,package_adapter_sha256=a.adapter_sha,
        previous_catalogue_sha256=a.catalogue_sha,previous_index_sha256=a.index_sha,
        previous_comparison_sha256=a.comparison_sha,
        catalogue_entries_before=len(by_n),catalogue_entries_after=len(catalogue),
        independent_replays=replays,receipt_sha256=receipt_manifest,
        scope='Exact construction audit before export; no new search or enormous coefficient expansion')))
    # No public files have been changed before this point. Verify every existing
    # destination either is one of the three replaceable manifests or absent.
    replaceable={repo/'catalogue.json',data_dir/'catalogue.json',data_dir/'index.json',comparison_path}
    for path,raw in staged.items():
        if path.exists() and path not in replaceable and path.read_bytes()!=raw:
            raise ValueError('Refusing to overwrite existing output: '+str(path))
    for path,raw in staged.items():
        path.parent.mkdir(parents=True,exist_ok=True)
        temporary=path.with_name(path.name+'.graded-tmp');temporary.write_bytes(raw);temporary.replace(path)
    print(json.dumps(dict(status='EXPORTED',records=records,files=len(staged),catalogue_entries=len(catalogue))))


if __name__=='__main__':main()
