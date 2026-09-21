"""Hash-checked, offline object storage with on-demand cache materialization."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

DATA=Path(__file__).parent/'data'
_INDEX=None
_VERIFIED=set()


def index():
    global _INDEX
    if _INDEX is None:
        _INDEX=json.loads((DATA/'index.json').read_text())
        if _INDEX['format']!='cwschemes-object-index-v1':raise ValueError('Unsupported object index')
    return _INDEX


def root():
    version=hashlib.sha256((DATA/'index.json').read_bytes()).hexdigest()[:16]
    base=Path(os.environ.get('CWSCHEMES_CACHE',Path.home()/'.cache/cwschemes'))
    path=base/version;path.mkdir(parents=True,exist_ok=True)
    return path


def raw_object(digest):
    if not re.fullmatch('[0-9a-f]{64}',digest):raise ValueError('Invalid object digest')
    meta=index()['objects'][digest]
    with gzip.open(DATA/'objects'/f'{digest}.gz','rb') as f:
        raw=f.read(meta['bytes']+1)
    if len(raw)!=meta['bytes'] or hashlib.sha256(raw).hexdigest()!=digest:
        raise ValueError(f'Object integrity failure: {digest}')
    return raw


def references(value):
    if isinstance(value,dict):
        for v in value.values():yield from references(v)
    elif isinstance(value,list):
        for v in value:yield from references(v)
    elif isinstance(value,str) and value.startswith('proofs/'):
        logical=value if value.endswith('.json') else value+'.json'
        if logical not in index()['paths']:raise ValueError(f'Unresolved proof: {logical}')
        yield logical


def materialize(logical, visiting=None):
    if logical not in index()['paths']:raise ValueError('Unknown resource')
    p=Path(logical)
    if p.is_absolute() or '..' in p.parts:raise ValueError('Unsafe resource path')
    digest=index()['paths'][logical];target=root()/logical
    if (logical,digest) in _VERIFIED and target.is_file():return target
    visiting=set() if visiting is None else visiting
    if digest in visiting:raise ValueError('Cyclic proof dependency')
    visiting.add(digest);raw=raw_object(digest)
    if logical.endswith('.json'):
        for ref in set(references(json.loads(raw))):materialize(ref,visiting)
    visiting.remove(digest)
    target.parent.mkdir(parents=True,exist_ok=True)
    if not target.exists() or hashlib.sha256(target.read_bytes()).hexdigest()!=digest:
        with tempfile.NamedTemporaryFile(dir=target.parent,delete=False) as f:f.write(raw);tmp=Path(f.name)
        os.replace(tmp,target)
    _VERIFIED.add((logical,digest));return target


def prepare(digest):
    for logical in index()['paths']:
        if logical.startswith('seeds/') or logical=='support14.json':materialize(logical)
    return materialize(f'proofs/{digest}.json')
