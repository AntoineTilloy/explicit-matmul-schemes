"""Checks for Z4-quotient units and squares (``z4_quotient``).

* ``audit``: support structure, relation-object structure (Z4-stable kept/zero/deleted sets, input modes never
  modified), and three independent rank counts (orbit unranking, Burnside, the closed-form polynomial).
* ``packet_identity``: exact rational evaluation (Fractions) of the packet tensor
  sum_{all columns} u(x) v(y) w(z) at sampled entries, compared with the 28 copies of <q^3>.
* ``orbit_lemma``: the input factor of a term is identical on every member of its column orbit (exact).
* ``tensor_entries``: full square tensor entries T[a,b,c] through the Strassen wiring, modulo primes.
* ``full_exact_check``: every n^6 coefficient of a small square (q = 2), modulo primes whose product exceeds a
  rigorous bound, hence exact over Q.
Packet-tensor evaluation factorises the label sums coordinate-wise: packet mode 0 (never modified by a relation
step, asserted) has the closed product form e0 / L_l = e0 + e_l + s/q / L* = e0 + 3s/(2q); packet modes 1 and 2 use
the generator's own coefficient tables.
"""
from fractions import Fraction as F
import math
import random
import time
from .z4_quotient import K9, hdig, hact, relations_object, Z4QuotientUnit
from .triangular_full_compression import allowed

POLY_GRAM = (1, 17, 122, 392, 493, 131, 2, 2, 0, 0)  # 4 R(q) = q^9 + 17 q^8 + ... + 2 q^2  (docs/construction.md)


def closed_form_rank(q):
    return sum(c * q ** (9 - i) for i, c in enumerate(POLY_GRAM)) // 4


def audit(unit):
    sup, pk, q = unit.sup, unit.packet, unit.q
    d, _, digest = relations_object(unit.relations_path)
    G = sup.G
    inv = lambda S: all(hact(g, H) in S for g in G for H in S)
    rows = [tuple(r) for r in d['rows']]
    initial = set(allowed(rows, 9, (0, 1)))
    assert d['initial_pair'] == [0, 1]
    current = initial
    free_c = d['perm'].index(2)
    for s in d['steps']:
        kept, deleted, zeros = set(s['kept']), set(s['deleted']), set(s['zero_patterns'])
        assert not (kept & deleted or kept & zeros or deleted & zeros) and kept | deleted | zeros == current
        assert inv(kept) and inv(deleted) and inv(zeros), 'relation step is not Z4-stable'
        if deleted:
            assert 3 - sum(s['pair']) == free_c, 'a relation step modifies an input mode'
            for I, b in s['blocks'].items():
                assert set(b['basis']) <= kept and set(b['deleted']) <= deleted
                assert len(b['coefficients']) == len(b['basis']) and all(len(r) == len(b['deleted']) for r in b['coefficients'])
        else:
            assert not s['blocks']
        current = kept
    assert current == set(pk.kept)
    # counts: unranking (per-orbit stabiliser count), Burnside over the whole kept set, closed form
    ordinary = pk.ordinary
    burnside = 0
    for g in G:
        for H in pk.kept:
            if hact(g, H) != H:
                continue
            ones = [t for t in range(K9) if ordinary[H] >> t & 1]
            seen, cycles = set(), 0
            for i in ones:
                if i in seen:
                    continue
                cycles += 1
                j = i
                while j not in seen:
                    seen.add(j)
                    j = g[j]
            burnside += q ** cycles
    assert burnside % len(G) == 0 and burnside // len(G) == unit.rank, 'Burnside count mismatch'
    assert closed_form_rank(q) == unit.rank, 'closed-form rank mismatch'
    unquotiented = sum(q ** ordinary[H].bit_count() for H in pk.kept)
    assert unquotiented == pk.rank
    return dict(q=q, rank=str(unit.rank), unquotiented_rank=str(unquotiented), support_rows=len(sup.rows), row_orbits=len(sup.orbs),
                kernel_profile=sup.kernel_profile, kept_patterns=len(pk.kept), relation_object=digest,
                counts='orbit unranking = Burnside = closed form')


# ------------------------------------------------------------------------------------------ packet tensor
def _zero(v):
    return sum(1 << t for t in range(K9) if v[t] == 0)


def _label_sums(q, v0, v1, v2, one, qi):
    S = []
    for t in range(K9):
        s = [[0, 0], [0, 0]]
        for l in range(1, q + 1):
            u = one if v0[t] == 0 else (one if v0[t] == l else 0) + qi
            d1 = (one if v1[t] == l else 0) - qi
            d2 = (one if v2[t] == l else 0) - qi
            s[0][0] += u
            s[1][0] += u * d1
            s[0][1] += u * d2
            s[1][1] += u * d1 * d2
        S.append(s)
    return S


def tdec_exact(unit, v0, v1, v2, perturb=None):
    """Packet tensor entry sum_{all columns} u(v0) v(v1) w(v2) at packet-mode variables, exactly over Q."""
    pk, q = unit.packet, unit.q
    z1, z2 = _zero(v1), _zero(v2)
    if z1 not in pk.families[1] or z2 not in pk.families[2]:
        return F(0)
    T1, T2 = pk.table(1, z1), pk.table(2, z2)
    qi = F(1, q)
    S = _label_sums(q, v0, v1, v2, F(1), qi)
    total = F(0)
    special = F(3, 2 * q)
    for H, tab1 in T1.items():
        tab2 = T2.get(H)
        if not tab1 or not tab2:
            continue
        d = hdig(H)
        base = F(1)
        for t in range(K9):
            if d[t] == 0 and v0[t] != 0:
                base = F(0)
                break
            if d[t] == 2 and v0[t] != 0:
                base *= special
        if not base:
            continue
        ords = [t for t in range(K9) if d[t] == 1]
        val = F(0)
        for I1, b1 in tab1.items():
            for I2, b2 in tab2.items():
                x = b1 * b2
                for t in ords:
                    x *= S[t][I1 >> t & 1][I2 >> t & 1]
                    if not x:
                        break
                val += x
        if perturb is not None:
            val += perturb
            perturb = None
        total += base * val
    return total


def _target(r1, r2, r3, v0, v1, v2):
    if not r1 == r2 == r3:
        return 0
    for t in range(K9):
        xs = (v0[t], v1[t], v2[t])
        if xs.count(0) != 1:
            return 0
        nz = [x for x in xs if x]
        if nz[0] != nz[1]:
            return 0
    return 1


def _draw_packet_entry(rows, q, kind, rng):
    if kind in ('valid', 'diag'):
        r1 = r2 = r3 = rng.choice(rows)
    elif kind == 'garbage':  # sunflower block triple: the correction / compression must cancel here
        while True:
            r1, r2, r3 = rng.choice(rows), rng.choice(rows), rng.choice(rows)
            a1, b2, c3 = r1[0], r2[1], r3[2]
            K = a1 & b2
            if (a1, b2, c3) not in rows and K and (a1 & c3) == K and (b2 & c3) == K:
                break
    else:
        r1, r2, r3 = rng.choice(rows), rng.choice(rows), rng.choice(rows)
    v0, v1, v2 = [0] * K9, [0] * K9, [0] * K9
    for t in range(K9):
        zs = (r1[0] >> t & 1, r2[1] >> t & 1, r3[2] >> t & 1)
        l = rng.randint(1, q)
        vals = [0 if zs[j] else (l if kind in ('valid', 'garbage') else rng.randint(1, q)) for j in range(3)]
        if kind == 'garbage' and sum(zs) == 0:
            vals = [rng.randint(1, q) for _ in range(3)]
        v0[t], v1[t], v2[t] = vals
    return r1, r2, r3, tuple(v0), tuple(v1), tuple(v2)


def packet_identity(unit, samples=8, seed=1, perturb=None):
    """Exact packet identity at sampled entries (valid, label-perturbed diagonal, sunflower garbage, uniform)."""
    rng = random.Random(seed)
    rows = [tuple(r) for r in unit.packet.rows]
    kinds = ['valid', 'garbage', 'diag', 'uniform']
    t0 = time.monotonic()
    records = []
    for k in range(samples):
        kind = kinds[k % 4]
        r1, r2, r3, v0, v1, v2 = _draw_packet_entry(rows, unit.q, kind, rng)
        target = _target(r1, r2, r3, v0, v1, v2)
        got = tdec_exact(unit, v0, v1, v2, perturb)
        records.append(dict(kind=kind, target=target, value=str(got), ok=got == target))
    return dict(samples=samples, targets_one=sum(r['target'] for r in records), mismatches=sum(not r['ok'] for r in records),
                passed=all(r['ok'] for r in records), exact=True, seconds=time.monotonic() - t0, records=records)


def orbit_lemma(unit, samples=6, seed=1):
    """For random terms, outer products h and entries: the input factor of every member of the term's column orbit
    coincides with the published (representative) unit coefficient, exactly."""
    rng = random.Random(seed)
    sup, q, pk = unit.sup, unit.q, unit.packet
    s = q ** 3
    checks = 0
    for _ in range(samples):
        term = rng.randrange(unit.rank)
        H, lab = unit.orbits.unrank(term)
        for mode in (0, 1):
            h, i1, i2 = rng.randrange(7), rng.randrange(s), rng.randrange(s)
            values = [sum((pk.entry(mode, H2, lab2, sup.local_var(rho, mode, i1, i2, q)) for rho in sup.orbs[h]), F(0))
                      for H2, lab2 in unit.orbits.members(H, lab)]
            assert all(v == values[0] for v in values), ('orbit lemma fails', term, mode)
            assert values[0] == unit.factor_entry(mode, term, unit.matrix_variable(h, mode, i1, i2))
            checks += 1
    return dict(checks=checks, passed=True)


# ------------------------------------------------------------------------------------------ square tensor mod p
def _fm(v, p):
    v = F(v)
    return v.numerator % p * pow(v.denominator % p, -1, p) % p


def tdec_mod(unit, v0, v1, v2, p, cache, perturb=0):
    import numpy as np
    pk, q = unit.packet, unit.q
    z1, z2 = _zero(v1), _zero(v2)
    if z1 not in pk.families[1] or z2 not in pk.families[2]:
        return 0
    def tab(j, z):
        if (j, z) not in cache:
            cache[(j, z)] = {H: (np.array(list(d.keys()), dtype=np.int64), np.array([_fm(v, p) for v in d.values()], dtype=np.int64))
                             for H, d in pk.table(j, z).items() if d}
        return cache[(j, z)]
    T1, T2 = tab(1, z1), tab(2, z2)
    qi = pow(q, -1, p)
    inv2q = pow(2 * q, -1, p)
    S = [np.array([[x % p for x in r] for r in s], dtype=np.int64) for s in _label_sums(q, v0, v1, v2, 1, qi)]
    total = 0
    for H, (I1, B1) in T1.items():
        if H not in T2:
            continue
        I2, B2 = T2[H]
        d = hdig(H)
        base = 1
        for t in range(K9):
            if d[t] == 0:
                base = base * (1 if v0[t] == 0 else 0) % p
            elif d[t] == 2:
                base = base * (1 if v0[t] == 0 else 3 * inv2q % p) % p
        if not base:
            continue
        M = (B1[:, None] * B2[None, :]) % p
        for t in range(K9):
            if d[t] == 1:
                M = M * S[t][((I1 >> t) & 1)[:, None], ((I2 >> t) & 1)[None, :]] % p
        val = int(M.sum() % p)
        if perturb:
            val = (val + perturb) % p
            perturb = 0
        total = (total + base * val) % p
    return total


def _garbage_entry(square, unit, rng):
    """(a, b, c) whose packet images contain a sunflower (off-row) block triple."""
    sup, q, n = unit.sup, unit.q, square.full_n
    s = q ** 3
    while True:
        r1, r2, r3 = rng.choice(sup.rowlist), rng.choice(sup.rowlist), rng.choice(sup.reps)
        a1, b2, c3 = r1[0], r2[1], r3[2]
        K = a1 & b2
        if (a1, b2, c3) != r1 and K and (a1 & c3) == K and (b2 & c3) == K:
            break
    x, y, z = [0] * K9, [0] * K9, [0] * K9
    for t in range(K9):
        zx, zy, zz = a1 >> t & 1, b2 >> t & 1, c3 >> t & 1
        if zx + zy + zz == 3:
            continue
        if zx + zy + zz == 0:
            x[t], y[t], z[t] = rng.randint(1, q), rng.randint(1, q), rng.randint(1, q)
        else:
            l = rng.randint(1, q)
            x[t], y[t], z[t] = (0 if zx else l), (0 if zy else l), (0 if zz else l)
    def inv_local(rho, h, v):
        r, g = sup.rowinfo[rho]
        a, b, c = sup.reps[r]
        first, second = (([i for i in range(9) if b >> i & 1], [i for i in range(9) if c >> i & 1]),
                         ([i for i in range(9) if c >> i & 1], [i for i in range(9) if a >> i & 1]),
                         ([i for i in range(9) if b >> i & 1], [i for i in range(9) if a >> i & 1]))[h]
        i1 = i2 = 0
        for t in first:
            i1 = i1 * q + (v[g[t]] - 1)
        for t in second:
            i2 = i2 * q + (v[g[t]] - 1)
        return i1, i2
    out = []
    for rho, h, v in ((r1, 0, x), (r2, 1, y), (r3, 2, z)):
        i1, i2 = inv_local(rho, h, v)
        orbit = sup.rowinfo[rho][0]
        blk = rng.choice([b for b in range(4) if square.outer.factor_entry(h, orbit, b // 2, b % 2)])
        out.append((blk // 2 * s + i1) * n + blk % 2 * s + i2)
    return tuple(out)


def tensor_entries(square, unit, samples=8, seed=1, primes=(2147483647, 2147483629), perturb=0):
    """Full square tensor entries T[a,b,c] = sum_t U_t(a) V_t(b) W_t(c) via the orbit lemma (sum over all packet
    columns) and coordinate-wise label sums, modulo primes; entries uniform / valid / near-valid / garbage."""
    rng = random.Random(seed)
    sup, q, n = unit.sup, unit.q, square.full_n
    s = q ** 3
    caches = {p: {} for p in primes}
    perm = unit.packet.perm
    t0 = time.monotonic()
    def loc(R, C, h, rows):
        br, i1 = divmod(R, s)
        bc, i2 = divmod(C, s)
        out = []
        for rho in rows:
            c = square.outer.factor_entry(h, sup.rowinfo[rho][0], br, bc)
            if c:
                out.append((int(c), sup.local_var(rho, h, i1, i2, q)))
        return out
    records = []
    for k in range(samples):
        kind = ('uniform', 'valid', 'near', 'garbage')[k % 4]
        if kind == 'uniform':
            A, B, C = rng.randrange(n * n), rng.randrange(n * n), rng.randrange(n * n)
        elif kind == 'garbage':
            A, B, C = _garbage_entry(square, unit, rng)
        else:
            i, kk, j = rng.randrange(n), rng.randrange(n), rng.randrange(n)
            A, B, C = i * n + kk, kk * n + j, i * n + j
            if kind == 'near':
                w, dd = rng.randrange(3), rng.randrange(1, n)
                if w == 0:
                    A = i * n + (kk + dd) % n
                elif w == 1:
                    B = kk * n + (j + dd) % n
                else:
                    C = ((i + dd) % n) * n + j
        (ai, ak), (bk, bj), (ci, cj) = divmod(A, n), divmod(B, n), divmod(C, n)
        target = int(ak == bk and ai == ci and bj == cj)
        X, Y, Z = loc(ai, ak, 0, sup.rowlist), loc(bk, bj, 1, sup.rowlist), loc(ci, cj, 2, sup.reps)
        residues = []
        for p in primes:
            tot = 0
            for ca, x in X:
                for cb, y in Y:
                    for cc, z in Z:
                        orig = (x, y, z)
                        v = tdec_mod(unit, *(orig[perm[j]] for j in range(3)), p, caches[p], perturb)
                        if v:
                            tot = (tot + ca * cb * cc * v) % p
            residues.append(tot)
        ok = all(r == target % p for r, p in zip(residues, primes))
        records.append(dict(a=A, b=B, c=C, kind=kind, target=target, residues=residues, ok=ok))
    return dict(samples=samples, primes=list(primes), mismatches=sum(not r['ok'] for r in records),
                passed=all(r['ok'] for r in records), seconds=time.monotonic() - t0, records=records)


# ------------------------------------------------------------------------------------------ all coefficients
def _primes_below(start, avoid):
    x = start
    while True:
        x -= 1
        if x < 3:
            raise RuntimeError('prime search exhausted')
        if all(x % d for d in range(2, math.isqrt(x) + 1)) and all(a % x for a in avoid):
            yield x


def factor_matrices(square, unit, p=None):
    """Dense factor matrices (rank x n^2, row-major) of a single-unit Strassen square, modulo p or in float64.

    Vectorised per (outer product, row, pattern-orbit representative): P(H, lab; x) = sum_I b_I prod_{i in I}
    ([lab_i = x_i] - 1/q) from the generator's tables; output factors sum over the column orbit (sum over the
    group divided by the column stabiliser)."""
    import numpy as np
    sup, q, pk = unit.sup, unit.q, unit.packet
    n, s, R = square.full_n, q ** 3, unit.rank
    dt = np.int64 if p else np.float64
    conv = (lambda v: _fm(v, p)) if p else float
    qi = pow(q, -1, p) if p else 1.0 / q
    terms = [unit.orbits.unrank(t) for t in range(R)]
    byH = {}
    for t, (H, lab) in enumerate(terms):
        byH.setdefault(H, []).append(t)
    out = [np.zeros((R, n * n), dtype=dt) for _ in range(3)]
    G = sup.G
    def pvals(j, zero, Hs, labs, X):
        """P for columns with patterns Hs (array per column) and labels labs (cols x 9) at variables X (v x 9)."""
        tabs = pk.table(j, zero)
        res = np.zeros((len(labs), len(X)), dtype=dt)
        for H in set(Hs.tolist()):
            sel = np.nonzero(Hs == H)[0]
            for I, val in tabs.get(H, {}).items():
                m = np.full((len(sel), len(X)), conv(val), dtype=dt)
                for i in range(K9):
                    if I >> i & 1:
                        f = (labs[sel, i][:, None] == X[None, :, i]).astype(dt) - qi
                        m = (m * (f % p)) % p if p else m * f
                res[sel] = (res[sel] + m) % p if p else res[sel] + m
        return res
    for mode in range(3):
        j = pk.perm.index(mode)
        for h in range(7):
            rows = sup.orbs[h] if mode < 2 else [sup.reps[h]]
            coef = [square.outer.factor_entry(mode, h, b // 2, b % 2) for b in range(4)]
            for rho in rows:
                X, i1, i2 = sup.index_arrays(rho, mode, q)
                zero = _zero(X[0])
                assert zero == rho[mode]
                if zero not in pk.families[j]:
                    continue
                for H, ts in byH.items():
                    ts = np.array(ts)
                    labs = np.array([terms[t][1] for t in ts], dtype=np.int64)
                    if mode < 2:
                        vals = pvals(j, zero, np.full(len(ts), H), labs, X)
                    else:
                        vals = np.zeros((len(ts), len(X)), dtype=dt)
                        stab = np.zeros(len(ts), dtype=np.int64)
                        for g in G:
                            H2 = hact(g, H)
                            labs2 = np.zeros_like(labs)
                            for t in range(K9):
                                labs2[:, g[t]] = labs[:, t]
                            v = pvals(j, zero, np.full(len(ts), H2), labs2, X)
                            vals = (vals + v) % p if p else vals + v
                            stab += (H2 == H) & np.all(labs2 == labs, axis=1)
                        if p:
                            inv = np.array([pow(int(x), -1, p) for x in stab], dtype=np.int64)
                            vals = vals * inv[:, None] % p
                        else:
                            vals = vals / stab[:, None]
                    for b in range(4):
                        if coef[b]:
                            pos = ((b // 2) * s + i1) * n + (b % 2) * s + i2
                            c = int(coef[b])
                            if p:
                                out[mode][np.ix_(ts, pos)] = (out[mode][np.ix_(ts, pos)] + c * vals) % p
                            else:
                                out[mode][np.ix_(ts, pos)] += c * vals
    return out


def full_exact_check(square, unit, extra_primes=0, perturb=False):
    """All n^6 coefficients of a single-unit Strassen square modulo primes p ~ 2^20 with prod p > 2 Delta (S+1):
    Delta = D_0 D_1 D_2 clears all denominators (D_h = lcm of table denominators times q^9) and
    S = sum_t max|U_t| max|V_t| max|W_t|, so every coefficient of Delta (sum_t U_t V_t W_t - T) is an integer of
    absolute value at most Delta (S + 1): zero modulo the primes means zero.  Feasible for q = 2."""
    import numpy as np
    t0 = time.monotonic()
    n, R, q = square.full_n, unit.rank, unit.q
    n2 = n * n
    D = [unit.packet.denominator(h) for h in range(3)]
    Delta = D[0] * D[1] * D[2]
    Uf = factor_matrices(square, unit)
    S = float(np.sum(np.abs(Uf[0]).max(axis=1) * np.abs(Uf[1]).max(axis=1) * np.abs(Uf[2]).max(axis=1))) * (1 + 1e-6)
    need = 2 * Delta * (S + 1)
    primes, prod = [], 1
    gen = _primes_below(1 << 20, [6 * q] + D)
    while prod <= need or len(primes) < extra_primes:
        primes.append(next(gen))
        prod *= primes[-1]
    mismatches = {}
    for p in primes:
        U, V, W = (M.astype(np.float64) for M in factor_matrices(square, unit, p))
        if perturb:
            U[0, 0] = (U[0, 0] + 1) % p
        bad = 0
        CH = 1 << 12
        for a0 in range(0, n2, 8):
            a1 = min(n2, a0 + 8)
            acc = np.zeros(((a1 - a0) * n2, n2))
            for r0 in range(0, R, CH):
                r1 = min(R, r0 + CH)
                KR = np.mod(U[r0:r1, a0:a1].T[:, None, :] * V[r0:r1].T[None, :, :], p).reshape(-1, r1 - r0)
                acc = np.mod(acc + KR @ W[r0:r1], p)
            acc = acc.reshape(a1 - a0, n2, n2)
            for t, a in enumerate(range(a0, a1)):
                i, j = divmod(a, n)
                E = np.zeros((n2, n2))
                E[j * n + np.arange(n), i * n + np.arange(n)] = 1
                bad += int(np.count_nonzero(np.mod(acc[t] - E, p)))
        mismatches[str(p)] = bad
    exact = all(v == 0 for v in mismatches.values())
    return dict(exact=exact, n=n, rank=R, primes=primes, log2_delta=round(math.log2(Delta), 1), log2_bound=round(math.log2(need), 1),
                log2_prime_product=round(math.log2(prod), 1), mismatches=mismatches, seconds=time.monotonic() - t0)
