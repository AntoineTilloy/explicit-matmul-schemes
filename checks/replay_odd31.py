"""Replay every coefficient of the bundled odd LITA31 primitive; requires g++."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import tempfile
import numpy as np
from cwschemes import _resources


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    digest='9effecdc238c8425993e7d73e3cefd45a40e67cd3cc08caf2825493f884e7955'
    primitive=_resources.materialize(f'rectbases/{digest}.npz')
    assert hashlib.sha256(primitive.read_bytes()).hexdigest()==digest
    with np.load(primitive,allow_pickle=False) as z:arrays={key:z[key] for key in z.files}
    meta=json.loads(arrays['metadata_json'].item())
    assert meta['tensor']==[31]*3 and meta['rank']==13433
    factors=[];scales=[]
    for axis in 'uvw':
        scale=math.lcm(*map(int,arrays[axis+'_denominators']));scales.append(scale)
        nums=[int(x)*(scale//int(d)) for x,d in zip(arrays[axis+'_numerators'],arrays[axis+'_denominators'])]
        assert all(-(2**63)<x<2**63 for x in nums)
        factors.append((np.asarray(arrays[axis+'_indptr'],dtype='<i8'),
                        np.asarray(arrays[axis+'_indices'],dtype='<i8'),np.asarray(nums,dtype='<i8')))
    cpp=root/'sources/verify_odd31.cpp'
    with tempfile.TemporaryDirectory(prefix='odd31-exact-') as directory:
        folder=Path(directory);binary=folder/'scaled.bin';exe=folder/'check'
        with binary.open('wb') as f:
            f.write(b'ODD12801');f.write(struct.pack('<5q',31,13433,*scales))
            for ptr,indices,nums in factors:
                f.write(struct.pack('<q',len(indices)));ptr.tofile(f);indices.tofile(f);nums.tofile(f)
        subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
        run=subprocess.run([str(exe),str(binary),str(folder)],check=True,capture_output=True,text=True)
        result=json.loads(run.stdout)
    assert result['exact'] and result['status']=='EXACT_ODD31_FULL_TENSOR_PASS'
    result.update(primitive_sha256=digest,native_source_sha256=hashlib.sha256(cpp.read_bytes()).hexdigest())
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
