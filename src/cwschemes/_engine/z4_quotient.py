"""Z4-quotient CW packets fed to an outer scheme (formats ``cw-z4-quotient-square-v1`` and packet
compression ``z4-quotient-v1`` inside ``cw-mixed-square-v1``).

Construction (see docs/construction.md).  A 28-row CW support at power N = 9 (m = 3) is invariant under a
coordinate permutation group Z4 and splits into 7 free Z4-orbits of rows.  The packet tensor
    T_pkt = G_q^{(x)9}|supp + sum_K beta_K (-D)^{(x)K} (x) G_q^{(x)(9-K)}|supp  (= 28 copies of <q^3>)
is written in the joint triangular representation (pure / ordinary / special columns (H, labels),
H in {0,1,2}^9) of the mode-permuted rows, then compressed by two certified relation steps:
    step 1  zero-pattern removal (patterns whose factor vanishes identically on the support);
    step 2  a Z4-invariant pair-(0,2) Gram elimination whose free mode is the OUTPUT C.
The coordinate permutation acts on columns; both input factors are equivariant, so on Z4-invariant inputs the
columns of one Z4-orbit give the same product.  One multiplication per column orbit (output factor = sum over
the orbit) computes the 7 row-orbit products independently: a scheme for 7 disjoint <q^3> products ("unit").
Fed to Strassen's <2,2,2;7> it is a <2q^3> scheme; inside a larger outer scheme, groups of 7 outer products
use one unit each.

Coefficients are exact Fractions computed lazily; per-(mode, zero-sector) coefficient tables are built from the
bundled q-independent relation object on first use and cached (``coefficient_cache/joint``).
Conventions: mode 0 = A (row-major), mode 1 = B, mode 2 = C = AB (direct output).
"""
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from itertools import accumulate, product
import hashlib
import json
import math
from .. import _resources
from .compact_orbit_packet import CompactOrbitPacket
from .triangular_full_compression import allowed

FORMAT = 'cw-z4-quotient-square-v1'
COMPRESSION = 'z4-quotient-v1'
RELATIONS_FORMAT = 'cw-z4-quotient-relations-v1'
K9 = 9


# ------------------------------------------------------------------------------------------ bundled objects
@lru_cache(None)
def _object_bytes(logical):
    return _resources.materialize(logical).read_bytes()


@lru_cache(None)
def support_object(logical):
    return Support(json.loads(_object_bytes(logical)), hashlib.sha256(_object_bytes(logical)).hexdigest())


@lru_cache(None)
def relations_object(logical):
    raw = _object_bytes(logical)
    d = json.loads(raw)
    if d.get('format') != RELATIONS_FORMAT or d.get('coefficient_field') != 'Q':
        raise ValueError('Unsupported Z4 relation object')
    steps = [dict(free=3 - sum(s['pair']), kept=frozenset(s['kept']),
                  blocks={int(I): dict(basis=b['basis'], deleted=b['deleted'],
                                       coefficients=[[F(x) for x in row] for row in b['coefficients']])
                          for I, b in s['blocks'].items()}) for s in d['steps']]
    return d, steps, hashlib.sha256(raw).hexdigest()


def replay_document(logical, q):
    """The relation object instantiated at q in the ``cw-compact-orbit-relations-v1`` layout (with the
    q-dependent step ranks recomputed), for the exact replay in ``verify.replay_orbit``."""
    d, _, _ = relations_object(logical)
    ordinary = [sum((v == 1) << i for i, v in enumerate(ds)) for ds in product(range(3), repeat=3 * d['m'])]
    rank = lambda kept: sum(q ** ordinary[H].bit_count() for H in kept)
    steps = [dict(pair=s['pair'], kept=s['kept'], deleted=s['deleted'], zero_patterns=s['zero_patterns'],
                  blocks=s['blocks'], rank=rank(s['kept'])) for s in d['steps']]
    return dict(format='cw-compact-orbit-relations-v1', coefficient_field='Q', q=q, m=d['m'], rows=d['rows'],
                initial_pair=d['initial_pair'], steps=steps, rank=steps[-1]['rank'], complete=True)


# ------------------------------------------------------------------------------------------ support / group
def closure(gens, N):
    idt = tuple(range(N))
    el = {idt}
    fr = [idt]
    while fr:
        nf = []
        for x in fr:
            for g in gens:
                y = tuple(g[x[i]] for i in range(N))
                if y not in el:
                    el.add(y)
                    nf.append(y)
        fr = nf
    return sorted(el)


def app(g, mask, N=K9):
    return sum(1 << g[i] for i in range(N) if mask >> i & 1)


def bits(mask, N=K9):
    return [i for i in range(N) if mask >> i & 1]


class Support:
    """28 rows (A, B, C zero masks), the group Z4 on the 9 coordinates, row orbits in stored order
    (so outer product r <-> row orbit r)."""

    def __init__(self, d, sha256):
        self.sha256 = sha256
        self.N = N = d['N']
        self.m = N // 3
        if N != K9:
            raise ValueError('Z4 support must have N = 9')
        self.G = closure([tuple(g) for g in d['gens']], N)
        self.rows = [tuple(r) for r in d['rows']]
        rs = set(self.rows)
        full = (1 << N) - 1
        for r in self.rows:
            assert all(bin(x).count('1') == self.m for x in r) and (r[0] | r[1] | r[2]) == full
            assert not (r[0] & r[1]) and not (r[0] & r[2]) and not (r[1] & r[2])
        assert len(rs) == len(self.rows) == 28
        self.orbs = []
        seen = set()
        for r in self.rows:
            if r in seen:
                continue
            o = sorted({tuple(app(g, z) for z in r) for g in self.G})
            seen.update(o)
            self.orbs.append(o)
            assert len(o) == len(self.G)
        assert len(self.orbs) == 7 and seen == rs
        self.reps = [o[0] for o in self.orbs]
        self.rowinfo = {}
        for r, o in enumerate(self.orbs):
            for g in self.G:
                self.rowinfo.setdefault(tuple(app(g, z) for z in self.reps[r]), (r, g))
        assert len(self.rowinfo) == 28
        self.rowlist = sorted(self.rowinfo)
        order = lambda g: next(k for k in range(1, 9) if _power(g, k) == tuple(range(N)))
        self.gen4 = [g for g in self.G if order(g) == 4]
        assert len(self.G) == 4 and self.gen4, 'Gamma must be Z4'
        g = self.gen4[0]
        self.cyc4 = []
        seen = set()
        for i in range(N):
            if i in seen or g[i] == i:
                continue
            c = [i]
            while g[c[-1]] != i:
                c.append(g[c[-1]])
            assert len(c) == 4
            seen.update(c)
            self.cyc4.append(c)
        self.fix = [i for i in range(N) if g[i] == i]
        # sunflower kernels of the support (pairs a in A, b in B with a matching c outside the rows); Moebius beta
        A = [r[0] for r in self.rows]
        B = [r[1] for r in self.rows]
        C = [r[2] for r in self.rows]
        kern = {}
        for a in A:
            for b in B:
                ab = a & b
                for c in C:
                    if (a & c) == ab and (b & c) == ab and (a, b, c) not in rs:
                        assert ab != 0, 'induced matching violated'
                        kern[ab] = kern.get(ab, 0) + 1
        beta = {}
        for K in sorted(kern, key=lambda K: (bin(K).count('1'), K)):
            s = 1
            for S, bS in beta.items():
                if S != K and (S & K) == S:
                    s += (-1) ** bin(S).count('1') * bS
            beta[K] = (-1) ** (bin(K).count('1') + 1) * s
        self.beta = {K: b for K, b in beta.items() if b}
        self.kernel_profile = [sum(1 for K in kern if bin(K).count('1') == j) for j in range(4)]

    def local_var(self, row, h, i1, i2, q):
        """Local variable in (0..q)^9 of matrix entry (i1, i2) of mode h in packet row ``row``."""
        r, g = self.rowinfo[row]
        a, b, c = self.reps[r]
        first, second = ((bits(b), bits(c)), (bits(c), bits(a)), (bits(b), bits(a)))[h]
        x = [0] * K9
        for idx, pos in ((i1, first), (i2, second)):
            for t in reversed(pos):
                idx, dgt = divmod(idx, q)
                x[g[t]] = dgt + 1
        return tuple(x)

    def index_arrays(self, row, h, q):
        """(s^2 x 9) local variables of all (i1, i2) of mode h in ``row`` (row-major), with i1 and i2."""
        import numpy as np
        r, g = self.rowinfo[row]
        a, b, c = self.reps[r]
        first, second = ((bits(b), bits(c)), (bits(c), bits(a)), (bits(b), bits(a)))[h]
        s = q ** self.m
        i1 = np.repeat(np.arange(s), s)
        i2 = np.tile(np.arange(s), s)
        X = np.zeros((s * s, K9), dtype=np.int64)
        for idx, pos in ((i1.copy(), first), (i2.copy(), second)):
            for t in reversed(pos):
                idx, dgt = np.divmod(idx, q)
                X[:, g[t]] = dgt + 1
        return X, i1, i2


def _power(g, k):
    x = tuple(range(len(g)))
    for _ in range(k):
        x = tuple(g[x[i]] for i in range(len(g)))
    return x


# ------------------------------------------------------------------------------------------ column orbits
def hdig(H):
    d = [0] * K9
    for i in range(K9 - 1, -1, -1):
        d[i] = H % 3
        H //= 3
    return tuple(d)


def hcode(d):
    c = 0
    for v in d:
        c = 3 * c + v
    return c


def hact(g, H):
    d = hdig(H)
    e = [0] * K9
    for t in range(K9):
        e[g[t]] = d[t]
    return hcode(e)


def lact(g, lab):
    e = [0] * K9
    for t in range(K9):
        e[g[t]] = lab[t]
    return tuple(e)


def pair_unrank(idx, M):
    """idx-th pair (x, y), 0 <= x < y < M, lexicographic."""
    lo, hi = 0, M - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid * M - mid * (mid + 1) // 2 <= idx:
            lo = mid
        else:
            hi = mid - 1
    x = lo
    y = x + 1 + (idx - (x * M - x * (x + 1) // 2))
    assert 0 <= x < y < M
    return x, y


def neck_count(A):
    return (A ** 4 + A ** 2 + 2 * A) // 4


def neck_unrank(idx, A):
    """idx-th Z4-orbit of words of length 4 over range(A); returns a representative word."""
    if idx < A:
        return (idx,) * 4
    idx -= A
    P2 = A * (A - 1) // 2
    if idx < P2:
        a, b = pair_unrank(idx, A)
        return (a, b, a, b)
    idx -= P2
    for mn in range(A):
        k = A - 1 - mn
        T = k ** 3 + k ** 2 + k * (k - 1) // 2 + k
        if idx < T:
            break
        idx -= T
    else:
        raise IndexError('necklace index out of range')
    if idx < k ** 3:
        x, r = divmod(idx, k * k)
        y, z = divmod(r, k)
        return (mn, mn + 1 + x, mn + 1 + y, mn + 1 + z)
    idx -= k ** 3
    if idx < k ** 2:
        y, z = divmod(idx, k)
        return (mn, mn, mn + 1 + y, mn + 1 + z)
    idx -= k ** 2
    if idx < k * (k - 1) // 2:
        x, z = pair_unrank(idx, k)
        return (mn, mn + 1 + x, mn, mn + 1 + z)
    idx -= k * (k - 1) // 2
    assert idx < k
    return (mn, mn, mn, mn + 1 + idx)


class OrbitIndex:
    """Z4-orbits of columns (H, labels), H in ``kept`` (a Z4-invariant pattern set); explicit unranking."""

    def __init__(self, sup, kept, q):
        self.sup = sup
        self.q = q
        G = sup.G
        ks = set(kept)
        assert all(hact(g, H) in ks for g in G for H in kept), 'kept set not Z4-invariant'
        self.preps = []
        seen = set()
        for H in sorted(kept):
            if H in seen:
                continue
            orb = {hact(g, H) for g in G}
            seen |= orb
            stab = [g for g in G if hact(g, H) == H]
            ones = [t for t, v in enumerate(hdig(H)) if v == 1]
            self.preps.append(dict(H=H, stab=stab, ones=ones, count=self._count(stab, ones)))
        self.ends = list(accumulate(p['count'] for p in self.preps))
        self.rank = self.ends[-1]

    def _count(self, stab, ones):
        tot = 0
        for g in stab:
            seen = set()
            c = 0
            for i in ones:
                if i in seen:
                    continue
                c += 1
                j = i
                while j not in seen:
                    seen.add(j)
                    j = g[j]
            tot += self.q ** c
        assert tot % len(stab) == 0
        return tot // len(stab)

    @lru_cache(4096)
    def unrank(self, term):
        """term -> representative column (H, labels); labels in 1..q on ordinary coordinates, else 0."""
        if not 0 <= term < self.rank:
            raise IndexError(term)
        j = bisect_right(self.ends, term)
        P = self.preps[j]
        idx = term - (self.ends[j - 1] if j else 0)
        q = self.q
        ones = P['ones']
        lab = [0] * K9
        stab = P['stab']
        if len(stab) == 1:
            for t in reversed(ones):
                idx, dgt = divmod(idx, q)
                lab[t] = dgt + 1
        elif len(stab) == 2:
            s2 = [g for g in stab if g != tuple(range(K9))][0]
            pairs = [(t, s2[t]) for t in ones if t < s2[t]]
            fixed = [t for t in ones if s2[t] == t]
            M = q ** len(pairs)
            Pn = M + M * (M - 1) // 2
            fidx, pidx = divmod(idx, Pn)
            for t in reversed(fixed):
                fidx, dgt = divmod(fidx, q)
                lab[t] = dgt + 1
            assert fidx == 0
            if pidx < M:
                X = Y = pidx
            else:
                X, Y = pair_unrank(pidx - M, M)
            for (t1, t2) in reversed(pairs):
                X, dx = divmod(X, q)
                Y, dy = divmod(Y, q)
                lab[t1] = dx + 1
                lab[t2] = dy + 1
        else:
            cyc = [c for c in self.sup.cyc4 if c[0] in ones]
            for c in cyc:
                assert all(t in ones for t in c)
            fixed = [t for t in ones if t in self.sup.fix]
            A = q ** len(cyc)
            Nn = neck_count(A) if cyc else 1
            fidx, nidx = divmod(idx, Nn)
            for t in reversed(fixed):
                fidx, dgt = divmod(fidx, q)
                lab[t] = dgt + 1
            assert fidx == 0
            if cyc:
                w = neck_unrank(nidx, A)
                for pos in range(4):
                    letter = w[pos]
                    for c in reversed(cyc):
                        letter, dgt = divmod(letter, q)
                        lab[c[pos]] = dgt + 1
        return P['H'], tuple(lab)

    def members(self, H, lab):
        return sorted({(hact(g, H), lact(g, lab)) for g in self.sup.G})


# ------------------------------------------------------------------------------------------ joint packet
class JointPacket(CompactOrbitPacket):
    """The package's ``CompactOrbitPacket`` (joint-orbit-v1) on the mode-permuted rows, driven by the
    q-independent Z4 relation object instead of a per-q proof file.  Original mode h is packet mode
    perm.index(h)."""

    def __init__(self, sup, q, relations):
        d, steps, digest = relations_object(relations)
        self.q = q
        self.m = d['m']
        self.k = 3 * self.m
        self.rows = [list(r) for r in d['rows']]
        self.patched = True
        self.relations_path = relations
        self.perm = perm = list(d['perm'])
        assert sorted(tuple(r[perm[j]] for j in range(3)) for r in sup.rows) == sorted(map(tuple, self.rows)), \
            'relation rows are not the permuted support rows'
        # table-cache key: relation object digest and q (tables depend on q)
        self.proof_bytes = f'{RELATIONS_FORMAT}:{digest}:q={q}'.encode()
        self.proof = d
        self.steps = steps
        self.families = [{r[h] for r in self.rows} for h in range(3)]
        self.initial_kept = allowed(self.rows, self.k, (0, 1))
        self.kept = d['steps'][-1]['kept']
        self.digits = list(product(range(3), repeat=self.k))
        self.ordinary = [sum((v == 1) << i for i, v in enumerate(ds)) for ds in self.digits]
        self.ends = list(accumulate(q ** self.ordinary[H].bit_count() for H in self.kept))
        self.rank = self.ends[-1]
        self.tables = {}
        self.pos = {H: j for j, H in enumerate(self.kept)}

    def index(self, H, lab):
        j = self.pos[H]
        loc = 0
        for t in range(K9):
            if self.ordinary[H] >> t & 1:
                loc = loc * self.q + (lab[t] - 1)
        return (self.ends[j - 1] if j else 0) + loc

    def entry(self, h, H, lab, x):
        return self.factor_entry(self.perm.index(h), self.index(H, lab), tuple(x))

    def table(self, j, zero):
        if (j, zero) not in self.tables:
            self.build_table(j, zero)
        return self.tables[(j, zero)]

    def denominator(self, h):
        """An integer D with D * (every mode-h entry) in Z: lcm of table denominators times q^9."""
        j = self.perm.index(h)
        L = 1
        for zero in self.families[j]:
            for tab in self.table(j, zero).values():
                for v in tab.values():
                    L = L * v.denominator // math.gcd(L, v.denominator)
        return L * self.q ** K9


# ------------------------------------------------------------------------------------------ the unit
class Z4QuotientUnit:
    """Scheme for 7 independent <q^3> products (outer product h <-> row orbit h), rank = number of column orbits.

    ``matrix_variable(h, mode, i1, i2)`` and ``factor_entry(mode, term, variable)`` follow the packet interface
    of ``general_square.MixedSquare``."""
    compression = COMPRESSION

    def __init__(self, q, support, relations):
        self.q = q
        self.m = 3
        self.support_path = support
        self.relations_path = relations
        self.sup = support_object(support)
        self.packet = JointPacket(self.sup, q, relations)
        self.orbits = OrbitIndex(self.sup, self.packet.kept, q)
        self.rank = self.orbits.rank
        self.rows = [list(r) for r in self.sup.reps]
        self.patched = True

    def matrix_variable(self, h, mode, row, col):
        s = self.q ** self.m
        if not 0 <= h < 7 or mode not in (0, 1, 2) or not (0 <= row < s and 0 <= col < s):
            raise IndexError((h, mode, row, col))
        return (h, row, col)

    def factor_entry(self, mode, term, variable):
        h, i1, i2 = variable
        H, lab = self.orbits.unrank(term)
        sup, q, pk = self.sup, self.q, self.packet
        if mode < 2:  # input factor of the representative column on the Z4-invariant embedding
            return sum((pk.entry(mode, H, lab, sup.local_var(rho, mode, i1, i2, q)) for rho in sup.orbs[h]), F(0))
        rep = sup.reps[h]  # output factor: sum over the column orbit at the representative row
        return sum((pk.entry(2, H2, lab2, sup.local_var(rep, 2, i1, i2, q)) for H2, lab2 in self.orbits.members(H, lab)), F(0))

    def spec(self):
        return dict(q=self.q, m=self.m, compression=COMPRESSION, rank=self.rank, support=self.support_path,
                    relations=self.relations_path)


@lru_cache(None)
def _unit(q, support, relations):
    return Z4QuotientUnit(q, support, relations)


def unit_from_spec(spec):
    u = _unit(spec['q'], spec['support'], spec['relations'])
    assert spec.get('m', 3) == 3
    if 'rank' in spec:
        assert u.rank == spec['rank'], 'Z4 unit rank mismatch'
    return u


# ------------------------------------------------------------------------------------------ square format
def load_certificate(data):
    """``cw-z4-quotient-square-v1``: one unit fed to the outer <2;7> (Strassen), optionally padded to n < 2q^3."""
    from .general_square import MixedSquare
    if data['format'] != FORMAT:
        raise ValueError(data['format'])
    unit = data['unit']
    assert data['outer_recipe'] == dict(kind='lille', n=2, rank=7)
    q = unit['q']
    s = MixedSquare(data['outer_recipe'], dict(kind='classical', n=q ** 3), [dict(unit, compression=COMPRESSION)],
                    data['n'], 0)
    assert s.full_n == data['full_n'] == 2 * q ** 3 and s.rank == data['rank'] and s.used == [7]
    return s


def square(q, support, relations, n=None):
    """The <2q^3> scheme (Strassen outer, one unit) at any q, from bundled support/relation objects."""
    return load_certificate(dict(format=FORMAT, n=n or 2 * q ** 3, full_n=2 * q ** 3, outer_recipe=dict(kind='lille', n=2, rank=7),
                                 unit=dict(q=q, m=3, support=support, relations=relations),
                                 rank=_unit(q, support, relations).rank))


def units(engine):
    """All Z4 units reachable from a coefficient engine (MixedSquare packets and nested 'cw' recipes)."""
    from .general_square import MixedSquare, Ordinary
    out = []
    def walk(x):
        if isinstance(x, MixedSquare):
            out.extend(p for p in x.packets if isinstance(p, Z4QuotientUnit))
            walk(x.outer)
            walk(x.leaf)
        elif isinstance(x, Ordinary):
            for name in ('child', 'left', 'right'):
                if hasattr(x, name):
                    walk(getattr(x, name))
    walk(engine)
    return out
