"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from fractions import Fraction as F
from functools import lru_cache
from math import prod
import json
from pathlib import Path
from .scheme import CompactPacket
from .published_seeds import RANKS, optimal_bases, entry as published_entry
STRASSEN = [[[1, 0, 1, 0, 1, -1, 0], [0, 0, 0, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1, 0], [1, 1, 0, 1, 0, 0, -1]], [[1, 1, 0, -1, 0, 1, 0], [0, 0, 1, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 1], [1, 0, -1, 0, 1, 0, 1]], [[1, 0, 0, 1, -1, 0, 1], [0, 0, 1, 0, 1, 0, 0], [0, 1, 0, 1, 0, 0, 0], [1, -1, 1, 0, 0, 1, 0]]]

class SquareCW:

    def __init__(self, q, rows, m, outer_levels=3, base4=False, pad_tail=False, use_lita=False, remainder_scheme=None):
        if q < 1 or q & q - 1:
            raise ValueError('self-contained wrapper requires q a power of two')
        if outer_levels < 0:
            raise ValueError('outer_levels must be nonnegative')
        self.packet = CompactPacket(q, rows, m, patched=False)
        self.q = q
        self.m = m
        self.outer_levels = outer_levels
        self.inner_levels = (q.bit_length() - 1) * m
        self.inner = q ** m
        self.outer = 2 ** outer_levels
        self.n = self.inner * self.outer
        bases = lambda levels: optimal_bases(levels, True) if use_lita else tuple(([2] if levels % 2 else []) + [4] * (levels // 2)) if base4 else (2,) * levels
        self.outer_bases = bases(outer_levels)
        self.inner_bases = bases(self.inner_levels)
        self.calls = prod((RANKS[b] for b in self.outer_bases))
        self.batch = len(rows)
        self.packets, self.leftovers = divmod(self.calls, self.batch)
        self.remainder_scheme = remainder_scheme
        if remainder_scheme is not None and remainder_scheme.n != self.inner:
            raise ValueError('remainder size mismatch')
        self.leaf_rank = remainder_scheme.rank if remainder_scheme else prod((RANKS[b] for b in self.inner_bases))
        self.padded_tail = bool(pad_tail and self.leftovers and (self.packet.rank < self.leftovers * self.leaf_rank))
        if self.padded_tail:
            self.packets += 1
            self.leftovers = 0
        self.rank = self.packets * self.packet.rank + self.leftovers * self.leaf_rank

    def outer_entry(self, mode, term, row, col):
        return product_entry(mode, self.outer_bases, term, row, col)

    def inner_entry(self, mode, term, row, col):
        if self.remainder_scheme is not None:
            return self.remainder_scheme.factor_entry(mode, term, row, col)
        return product_entry(mode, self.inner_bases, term, row, col)

    def factor_entry(self, mode, index, row, col):
        if mode not in (0, 1, 2) or not (0 <= index < self.rank and 0 <= row < self.n and (0 <= col < self.n)):
            raise IndexError((mode, index, row, col))
        br, ir = divmod(row, self.inner)
        bc, ic = divmod(col, self.inner)
        if index < self.packets * self.packet.rank:
            batch, term = divmod(index, self.packet.rank)
            value = F(0)
            for h in range(self.batch):
                if batch * self.batch + h >= self.calls:
                    break
                outer = self.outer_entry(mode, batch * self.batch + h, br, bc)
                if outer:
                    variable = self.packet.matrix_variable(h, mode, ir, ic)
                    value += outer * self.packet.factor_entry(mode, term, variable)
            return value
        remainder = index - self.packets * self.packet.rank
        call, term = divmod(remainder, self.leaf_rank)
        return F(self.outer_entry(mode, self.packets * self.batch + call, br, bc) * self.inner_entry(mode, term, ir, ic))

    def certificate(self):
        return dict(format='cw-single-square-v1', coefficient_field='Q', n=self.n, rank=self.rank, q=self.q, m=self.m, outer_levels=self.outer_levels, inner_levels=self.inner_levels, outer_bases=self.outer_bases, inner_bases=self.inner_bases, outer_calls=self.calls, packet_count=self.packets, products_per_packet=self.batch, padded_tail=self.padded_tail, zero_filled_products=max(0, self.packets * self.batch - self.calls), remainder_scheme=self.remainder_scheme.certificate() if self.remainder_scheme else None, leftover_products=self.leftovers, leftover_rank=self.leaf_rank, packet=self.packet.certificate(), coefficient_generator='square_scheme.SquareCW.factor_entry')

@lru_cache(None)
def base_factors(n):
    if n == 2:
        return STRASSEN
    if n == 4:
        data = json.loads((resource_root() / 'seeds/dps48.json').read_text())
        return [[[F(v) for v in row] for row in matrix] for matrix in data['factors']]
    raise ValueError(n)

def product_entry(mode, bases, term, row, col):
    value = F(1)
    for base in reversed(bases):
        rank = RANKS[base]
        term, t = divmod(term, rank)
        row, i = divmod(row, base)
        col, j = divmod(col, base)
        value *= published_entry(base, mode, t, base * i + j) if base in (16, 32) else base_factors(base)[mode][base * i + j][t]
        if not value:
            return F(0)
    return value
