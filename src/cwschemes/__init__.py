"""Exact, lazy coefficient access to the published catalogue."""
from fractions import Fraction
import json
from . import _resources

__version__='0.1.0'


def catalogue():
    """Return catalogue records; dimensions and ranks are decimal strings."""
    return json.loads((_resources.DATA/'catalogue.json').read_text())


def load(n):
    """Load one catalogue size without allocating its full coefficient arrays."""
    record=next((r for r in catalogue() if int(r['n'])==int(n)),None)
    if record is None:raise KeyError(f'No saved scheme at size {n}')
    path=_resources.prepare(record['root_object'])
    data=json.loads(path.read_text())
    if data['format'] in ('cw-single-square-v1','ordinary-square-product-v1'):
        from ._engine.legacy_square import load_certificate
    else:
        from ._engine.general_square import load_certificate
    engine=load_certificate(data)
    if (engine.n,engine.rank)!=(int(record['n']),int(record['rank'])):raise ValueError('Catalogue count mismatch')
    return Scheme(record,engine)


class Scheme:
    def __init__(self,record,engine):
        self.record=record;self._engine=engine;self.n=engine.n;self.rank=engine.rank;self.field='Q'

    def coefficient(self,factor,term,row,col):
        """Return a Fraction in the direct-output U/V/W convention."""
        mode={'U':0,'V':1,'W':2}.get(factor,factor)
        if mode not in (0,1,2):raise ValueError('Factor must be U, V or W')
        if any(type(x) is not int for x in (term,row,col)):raise TypeError('Indices must be integers')
        if not (0<=term<self.rank and 0<=row<self.n and 0<=col<self.n):raise IndexError('Coefficient index out of range')
        return Fraction(self._engine.factor_entry(mode,term,row,col))

    def coefficients(self,factor,queries,*,max_queries=1000):
        """Evaluate a bounded batch; never allocate an entire factor matrix."""
        from itertools import islice
        if type(max_queries) is not int or max_queries<0:raise ValueError('max_queries must be a nonnegative integer')
        queries=list(islice(iter(queries),max_queries+1))
        if len(queries)>max_queries:raise ValueError('Requested batch exceeds max_queries')
        return [self.coefficient(factor,*q) for q in queries]
