import json
from fractions import Fraction
from pathlib import Path
import cwschemes
from cwschemes.verify import verify

assert '/wheel_install/' in str(Path(cwschemes.__file__).resolve())
records=cwschemes.catalogue()
assert len(records)==99
graded=[r for r in records if 'graded' in r['provenance']['construction'].lower()]
assert len(graded)==2
results=[]
for record in graded:
    n=int(record['n']);scheme=cwschemes.load(n)
    assert scheme.rank==int(record['rank']) and scheme.rank**20<n**53
    try:
        scheme.coefficient('U',0,0,0)
    except OverflowError:
        pass
    else:
        raise AssertionError('Large coefficient query was not guarded')
    audit=verify(n,exact=True)
    assert audit['sub_2_65'] and all(not r['engine_loaded'] for r in audit['graded_construction']['base_references'])
    results.append(dict(n=str(n),rank=str(scheme.rank),strict_sub265=True,verification=audit['mode']))
old=cwschemes.load(32)
assert old.rank==14196 and isinstance(old.coefficient('U',0,0,0),Fraction)
assert next(int(r['rank']) for r in records if int(r['n'])==1024)==198679192
print(json.dumps(dict(status='PASS',installed_module=cwschemes.__file__,catalogue_entries=len(records),
                     graded=results,old_size32_coefficient='PASS',size1024_rank=198679192),indent=2))
