"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
from .scheme import CompactPacket
from .square_scheme import base_factors
ROOT = resource_root()

@lru_cache(None)
def primitive(n, rank=None):
    if n in (2, 4):
        return base_factors(n)
    filename = 'simple7_249.json' if n == 7 and rank == 249 else f'simple{n}.json'
    data = json.loads((ROOT / 'seeds' / filename).read_text())
    return [[[F(x) for x in row] for row in matrix] for matrix in data['factors']]

class Ordinary:

    def __init__(self, recipe):
        self.recipe = recipe
        self.n = recipe['n']
        self.kind = recipe['kind']
        if self.kind == 'classical':
            self.rank = self.n ** 3
        elif self.kind == 'lille':
            self.rank = recipe['rank']
        elif self.kind == 'lita':
            from .sources.lita_round6 import lita_rank
            self.rank = lita_rank(self.n)
        elif self.kind == 'cw':
            self.child = load_certificate(recipe['certificate'])
            assert self.child.n == self.n
            self.rank = self.child.rank
        elif self.kind == 'product':
            self.left = Ordinary(recipe['left'])
            self.right = Ordinary(recipe['right'])
            assert self.n == self.left.n * self.right.n
            self.rank = self.left.rank * self.right.rank
        elif self.kind in ('pad', 'peel'):
            self.child = Ordinary(recipe['child'])
            assert self.child.n > self.n if self.kind == 'pad' else self.child.n == self.n - 1
            self.rank = self.child.rank
            if self.kind == 'peel':
                self.rank += 3 * self.n * self.n - 3 * self.n + 1
        else:
            raise ValueError(self.kind)

    def factor_entry(self, mode, term, row, col):
        if mode not in (0, 1, 2) or not (0 <= term < self.rank and 0 <= row < self.n and (0 <= col < self.n)):
            raise IndexError((mode, term, row, col))
        n = self.n
        if self.kind == 'cw':
            return self.child.factor_entry(mode, term, row, col)
        if self.kind == 'lita':
            from .lita_coefficients import entry
            return entry(n, mode, term, row, col)
        if self.kind == 'lille':
            if n in (2, 3, 4, 5, 6, 7):
                factors = primitive(n, self.rank)
                assert len(factors[0][0]) == self.rank
                return F(factors[mode][n * row + col][term])
            if n >= 8 and n % 2 == 0 and (not (ROOT / f'seeds/lita{n}.npz').exists()):
                from .sources.lita_round6 import lita_rank
                from .lita_coefficients import entry as formula_entry
                assert self.rank == lita_rank(n)
                return formula_entry(n, mode, term, row, col)
            from .published_seeds import entry, load
            assert len(load(n)['u_indptr']) - 1 == self.rank
            return entry(n, mode, term, n * row + col)
        if self.kind == 'product':
            a, b = divmod(term, self.right.rank)
            i, x = divmod(row, self.right.n)
            j, y = divmod(col, self.right.n)
            return self.left.factor_entry(mode, a, i, j) * self.right.factor_entry(mode, b, x, y)
        if self.kind == 'pad':
            return self.child.factor_entry(mode, term, row, col)
        if self.kind == 'peel':
            if term < self.child.rank:
                return self.child.factor_entry(mode, term, row, col) if row < n - 1 and col < n - 1 else F(0)
            term -= self.child.rank
            if term < n * n:
                i = n - 1
                k, j = divmod(term, n)
            elif term < n * n + (n - 1) * n:
                i, j = divmod(term - n * n, n)
                k = n - 1
            else:
                i, k = divmod(term - n * n - (n - 1) * n, n - 1)
                j = n - 1
        else:
            i, rest = divmod(term, n * n)
            k, j = divmod(rest, n)
        return F((row, col) == ((i, k), (k, j), (i, j))[mode])

class MixedSquare:

    def __init__(self, outer_recipe, leaf_recipe, packets, target=None, leftover_products=None):
        self.outer = Ordinary(outer_recipe)
        self.leaf = Ordinary(leaf_recipe)
        self.inner = self.leaf.n
        self.full_n = self.outer.n * self.inner
        self.n = self.full_n if target is None else target
        assert 1 <= self.n <= self.full_n
        self.specs = packets
        self.packets = []
        self.ends = []
        self.starts = []
        self.used = []
        budget = self.outer.rank if leftover_products is None else self.outer.rank - leftover_products
        assert 0 <= budget <= self.outer.rank
        offset = calls = 0
        for spec in packets:
            if spec.get('compression') == 'z4-quotient-v1':
                from .z4_quotient import unit_from_spec
                p = unit_from_spec(spec)
            elif spec.get('compression') == 'tailored-completion-v1':
                from .tailored_completion import TailoredPacket
                p = TailoredPacket(spec['q'], spec['rows'], spec['m'])
            elif spec.get('compression') == 'refined-orbit-v1':
                from .refined_orbit_packet import RefinedOrbitPacket
                p = RefinedOrbitPacket(spec['q'], spec['rows'], spec['m'], spec['relations'])
            elif spec.get('compression') == 'joint-orbit-v1':
                from .compact_orbit_packet import CompactOrbitPacket
                p = CompactOrbitPacket(spec['q'], spec['rows'], spec['m'], spec['relations'])
            elif spec.get('compression') == 'orbit-iterated-v1':
                from .iterated_orbit_packet import IteratedOrbitPacket
                p = IteratedOrbitPacket(spec['q'], spec['rows'], spec['m'], spec['relations'])
            elif spec.get('compression') == 'pair-schur-multi-v1':
                from .second_schur_packet import MultiSchurPacket
                p = MultiSchurPacket(spec['q'], spec['rows'], spec['m'], spec['relations'])
            elif spec.get('compression') == 'pair-schur-second-v1':
                from .second_schur_packet import SecondSchurPacket
                p = SecondSchurPacket(spec['q'], spec['rows'], spec['m'], spec['relations'])
            elif spec.get('compression') == 'pair-schur-v1':
                from .schur_compression import SchurPacket
                p = SchurPacket(spec['q'], spec['rows'], spec['m'])
            else:
                assert not spec.get('compression'), 'unknown compression format'
                p = CompactPacket(spec['q'], spec['rows'], spec['m'], patched=spec.get('patched', True))
            assert p.q ** p.m == self.inner
            if 'rank' in spec:
                assert p.rank == spec['rank']
            assert calls < budget, 'Only the final packet may overfill'
            self.starts.append(calls)
            used = min(len(p.rows), budget - calls)
            self.used.append(used)
            calls += used
            offset += p.rank
            self.ends.append(offset)
            self.packets.append(p)
        self.packet_rank = offset
        self.covered = calls
        self.leaves = self.outer.rank - calls
        if leftover_products is not None:
            assert self.leaves == leftover_products
        self.rank = offset + self.leaves * self.leaf.rank

    def factor_entry(self, mode, term, row, col):
        if mode not in (0, 1, 2) or not (0 <= term < self.rank and 0 <= row < self.n and (0 <= col < self.n)):
            raise IndexError((mode, term, row, col))
        br, ir = divmod(row, self.inner)
        bc, ic = divmod(col, self.inner)
        if term >= self.packet_rank:
            call, t = divmod(term - self.packet_rank, self.leaf.rank)
            return self.outer.factor_entry(mode, self.covered + call, br, bc) * self.leaf.factor_entry(mode, t, ir, ic)
        batch = bisect_right(self.ends, term)
        t = term - (self.ends[batch - 1] if batch else 0)
        p = self.packets[batch]
        value = F(0)
        for h in range(self.used[batch]):
            a = self.outer.factor_entry(mode, self.starts[batch] + h, br, bc)
            if a:
                v = p.matrix_variable(h, mode, ir, ic)
                value += a * p.factor_entry(mode, t, v)
        return value

    def certificate(self):
        return dict(format='cw-mixed-square-v1', coefficient_field='Q', n=self.n, full_n=self.full_n, rank=self.rank, outer_recipe=self.outer.recipe, leaf_recipe=self.leaf.recipe, packets=[dict(q=p.q, m=p.m, rows=p.rows, patched=p.patched, rank=p.rank, **{'compression': p.compression} if hasattr(p, 'compression') else {}, **{'relations': p.relations_path} if hasattr(p, 'relations_path') else {}) for p in self.packets], packet_products_used=self.used, leftover_products=self.leaves, coefficient_generator='general_square.load_certificate')

def load_certificate(value):
    data = json.loads(Path(value).read_text()) if isinstance(value, (str, Path)) else value
    if data['format'] == 'cw-z4-quotient-square-v1':
        from .z4_quotient import load_certificate as load_z4
        return load_z4(data)
    if data['format'] == 'explicit-qcsr-square-v1':
        from .explicit_square import load_certificate as load_explicit
        return load_explicit(data)
    if data['format'] == 'cw-compressed-fixed-quotient-square-v1':
        from .compressed_fixed_quotient import load_square_certificate
        return load_square_certificate(data)
    if data['format'] == 'cw-fixed-quotient-square-v1':
        from .fixed_quotient_square import load_certificate as load_fixed
        return load_fixed(data)
    assert data['format'] == 'cw-mixed-square-v1'
    s = MixedSquare(data['outer_recipe'], data['leaf_recipe'], data['packets'], data['n'], data['leftover_products'])
    assert s.rank == data['rank'] and s.full_n == data['full_n']
    assert s.used == data['packet_products_used']
    assert s.leaves == data['leftover_products']
    return s
