"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
import json, math
from pathlib import Path
from .square_scheme import SquareCW, product_entry
from .published_seeds import optimal_bases, RANKS

class OrdinaryScheme:

    def __init__(self, levels, use_lita=True):
        self.n = 2 ** levels
        self.bases = optimal_bases(levels, use_lita)
        self.rank = math.prod((RANKS[b] for b in self.bases))

    def factor_entry(self, mode, index, row, col):
        return product_entry(mode, self.bases, index, row, col)

    def certificate(self):
        return dict(format='ordinary-square-product-v1', coefficient_field='Q', n=self.n, rank=self.rank, bases=self.bases, coefficient_generator='optimise_hierarchy.OrdinaryScheme.factor_entry')

def load_certificate(path_or_data):
    """Reconstruct a delivered certificate, checking its counts and base order."""
    data = path_or_data if isinstance(path_or_data, dict) else json.loads(Path(path_or_data).read_text())
    if data['format'] == 'ordinary-square-product-v1':
        n = data['n']
        assert n > 0 and n & n - 1 == 0
        scheme = OrdinaryScheme(n.bit_length() - 1, any((b in (16, 32) for b in data['bases'])))
        assert list(scheme.bases) == data['bases']
    else:
        assert data['format'] == 'cw-single-square-v1'
        child = load_certificate(data['remainder_scheme']) if data.get('remainder_scheme') else None
        bases = data['outer_bases'] + data['inner_bases']
        w = data['packet']['witness']
        scheme = SquareCW(data['q'], w['rows'], data['m'], data['outer_levels'], 4 in bases, data.get('padded_tail', False), any((b in (16, 32) for b in bases)), child)
        assert list(scheme.outer_bases) == data['outer_bases'] and list(scheme.inner_bases) == data['inner_bases']
        assert scheme.packet.rank == data['packet']['rank']
    assert (scheme.n, scheme.rank) == (data['n'], data['rank'])
    return scheme
