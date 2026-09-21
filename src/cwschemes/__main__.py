import argparse
import json
from . import catalogue, load


def main():
    p=argparse.ArgumentParser(description='Explicit rational square matmul catalogue')
    sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('list')
    info=sub.add_parser('info');info.add_argument('n',type=int)
    q=sub.add_parser('coefficient');q.add_argument('n',type=int);q.add_argument('factor',choices=['U','V','W'])
    for name in ('term','row','col'):q.add_argument(name,type=int)
    v=sub.add_parser('verify');v.add_argument('n',type=int,nargs='?');v.add_argument('--all',action='store_true')
    v.add_argument('--exact',action='store_true');v.add_argument('--threads',type=int,default=1)
    a=p.parse_args()
    if a.command=='list':print(json.dumps(catalogue(),indent=2))
    elif a.command=='info':
        r=next(r for r in catalogue() if int(r['n'])==a.n);print(json.dumps(r,indent=2))
    elif a.command=='coefficient':print(load(a.n).coefficient(a.factor,a.term,a.row,a.col))
    else:
        from .verify import verify
        if not a.all and a.n is None:p.error('Specify a size or --all')
        sizes=[int(r['n']) for r in catalogue()] if a.all else [a.n]
        for n in sizes:print(json.dumps(verify(n,exact=a.exact,threads=a.threads)),flush=True)


if __name__=='__main__':main()
