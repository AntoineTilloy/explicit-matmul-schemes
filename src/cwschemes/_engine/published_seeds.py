"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
import json
RANKS = {2: 7, 4: 48, 16: 2247, 32: 14215}
PUBLISHED_RANKS = {13: 1420, 14: 1603, 16: 2247, 18: 3043, 19: 3981, 20: 4007, 21: 5160, 22: 5155, 23: 6545, 24: 6503, 25: 8152, 26: 8067, 27: 9997, 28: 9863, 29: 12096, 30: 11907, 31: 14465, 32: 14215, 34: 16803, 86: 234787}

@lru_cache(None)
def load(n):
    import numpy as np
    if n not in PUBLISHED_RANKS:
        raise ValueError(n)
    path = resource_root() / f'seeds/lita{n}.npz'
    with np.load(path, allow_pickle=False) as data:
        out = {k: data[k] for k in data.files}
    meta = json.loads(out['metadata_json'].item())
    assert meta['rank'] == PUBLISHED_RANKS[n] and meta['tensor'] == [n, n, n] and (meta['coefficient_field'] == 'Q')
    return out

def entry(n, mode, term, variable):
    import numpy as np
    data = load(n)
    axis = 'uvw'[mode]
    ptr = data[axis + '_indptr']
    start, end = (int(ptr[term]), int(ptr[term + 1]))
    indices = data[axis + '_indices']
    offset = start + int(np.searchsorted(indices[start:end], variable))
    if offset == end or indices[offset] != variable:
        return F(0)
    return F(int(data[axis + '_numerators'][offset]), int(data[axis + '_denominators'][offset]))

def optimal_bases(levels, use_lita=False):
    choices = [2, 4, 16, 32] if use_lita else [2, 4]
    best = [(1, ())]
    for k in range(1, levels + 1):
        options = []
        for n in choices:
            l = n.bit_length() - 1
            if l <= k:
                options.append((best[k - l][0] * RANKS[n], tuple(sorted(best[k - l][1] + (n,), reverse=True))))
        best.append(min(options))
    return best[levels][1]
