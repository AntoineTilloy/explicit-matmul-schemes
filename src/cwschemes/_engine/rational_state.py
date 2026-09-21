"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from itertools import product
import numpy as np
from .rational_arrays import Q, ZERO, qm, from_qm
from .invariant_state import State
from .triangular_full_compression import allowed

def triangular_apply(vector, sizes, q):
    v = vector.reshape(sizes)
    for axis in range(len(sizes)):
        x = np.moveaxis(v, axis, 0)
        y = x.copy()
        y[0] = x[0] - q * x[1] - Q(1, 4) * x[2]
        y[1] = x[1] - Q(3, 4 * q) * x[2]
        v = np.moveaxis(y, 0, axis)
    return v.ravel()

def initial_free(q, k, zero, kept):
    f = np.array([Q(1)], dtype=object)
    codes = np.zeros(1, dtype=np.int64)
    standards = np.zeros(1, dtype=np.int64)
    sizes = []
    for i in range(k):
        if zero >> i & 1:
            raw = [Q(q, 3), Q(1), Q(-4 * q, 3)]
            types = [0, 1, 2]
            standard = [0, 0, 0]
        else:
            raw = [ZERO, Q(2, q), Q(-2), Q(1)]
            types = [0, 1, 2, 1]
            standard = [0, 0, 0, 1 << i]
        sizes.append(len(raw))
        f = np.kron(f, np.array(raw, dtype=object))
        codes = (3 * codes[:, None] + np.array(types)).ravel()
        standards = (standards[:, None] + np.array(standard)).ravel()
    keep = np.zeros(3 ** k, dtype=bool)
    keep[list(kept)] = True
    removed = ~keep[codes]
    term = np.where(removed, f, ZERO)
    x = term.copy()
    for _ in range(2 * k):
        term = np.where(removed, term - triangular_apply(term, tuple(sizes), q), ZERO)
        x += term
        if not np.any(term):
            break
    result = f - triangular_apply(x, tuple(sizes), q)
    assert not np.any(result[removed])
    return (result, codes, standards)

class ExactState(State):

    def __init__(self, q, rows, m):
        self.q = q
        self.rows = rows
        self.m = m
        self.k = 3 * m
        k = self.k
        self.families = [sorted({r[h] for r in rows}) for h in range(3)]
        self.digits = list(product(range(3), repeat=k))
        self.ordinary = [sum(((v == 1) << i for i, v in enumerate(ds))) for ds in self.digits]
        self.kept = sorted(allowed(rows, k, (0, 1)))
        self.reindex()
        self.coefs = [[np.full((len(self.families[mode]), len(ps)), ZERO, dtype=object) for ps in self.patterns] for mode in range(3)]
        pure = [sum(((v == 0) << i for i, v in enumerate(ds))) for ds in self.digits]
        special = [sum(((v == 2) << i for i, v in enumerate(ds))) for ds in self.digits]
        for mode in (0, 1):
            for I, ps in enumerate(self.patterns):
                for z, Z in enumerate(self.families[mode]):
                    if Z & I:
                        continue
                    for j, H in enumerate(ps):
                        if pure[H] & ~Z:
                            continue
                        self.coefs[mode][I][z, j] = Q(2, q) ** (self.ordinary[H] & ~Z & ~I).bit_count() * Q(3, 2 * q) ** (special[H] & ~Z).bit_count()
        for z, Z in enumerate(self.families[2]):
            values, codes, standards = initial_free(q, k, Z, self.kept)
            for address, H in enumerate(codes):
                I = int(standards[address])
                j = self.lookups[I].get(int(H))
                if j is not None:
                    self.coefs[2][I][z, j] = values[address]
            print('exact initial sector', z + 1, 'of', len(self.families[2]), flush=True)

    def grams(self, mode):
        out = []
        q = self.q
        for I, ps in enumerate(self.patterns):
            F = self.coefs[mode][I]
            if not ps:
                out.append(np.empty((0, 0), dtype=object))
                continue
            weights = np.array([Q(q) ** (self.ordinary[H].bit_count() - I.bit_count()) for H in ps], dtype=object)
            G = from_qm(qm(F).transpose() * qm(F))
            out.append(G * Q(q) ** (self.k - self.m - I.bit_count()) * weights[None, :])
        return out

    def pair_gram(self, pair):
        a = self.grams(pair[0])
        b = self.grams(pair[1])
        out = [np.full(G.shape, ZERO, dtype=object) for G in a]
        buckets = {}
        q = self.q
        for H in self.kept:
            for K in self.kept:
                buckets.setdefault(self.ordinary[H] & self.ordinary[K], []).append((H, K))

        def transform(v, ell, inverse):
            x = v.reshape((2,) * ell + (v.shape[-1],))
            for axis in range(ell):
                a = np.moveaxis(x, axis, 0)
                b = a.copy()
                b[0] = a[0] + (q - 1) * a[1]
                b[1] = a[0] - a[1]
                if inverse:
                    b *= Q(1, q)
                x = np.moveaxis(b, 0, axis)
            return x.reshape(v.shape)
        for L, pairs in buckets.items():
            positions = [i for i in range(self.k) if L >> i & 1]
            Is = [sum((bit << i for bit, i in zip(bits, positions))) for bits in product((0, 1), repeat=len(positions))]
            va = []
            vb = []
            indices = []
            for I in Is:
                rr = np.array([self.lookups[I][H] for H, K in pairs])
                cc = np.array([self.lookups[I][K] for H, K in pairs])
                indices.append((rr, cc))
                va.append(a[I][rr, cc])
                vb.append(b[I][rr, cc])
            va = transform(np.array(va, dtype=object), len(positions), True)
            vb = transform(np.array(vb, dtype=object), len(positions), True)
            result = transform(va * vb, len(positions), False)
            scale = np.array([Q(q) ** (-(self.ordinary[K] & ~L).bit_count()) for H, K in pairs], dtype=object)
            result *= scale[None, :]
            for j, I in enumerate(Is):
                out[I][indices[j]] = result[j]
        return out

    def eliminate_exact(self, G, selection, pair):
        free = 3 - sum(pair)
        Kset = set(selection['kept'])
        Jset = set(selection['deleted'])
        new = [[] for _ in range(3)]
        blocks = {}
        checks = 0
        for I, ps in enumerate(self.patterns):
            K = [j for j, H in enumerate(ps) if H in Kset]
            J = [j for j, H in enumerate(ps) if H in Jset]
            for mode in range(3):
                new[mode].append(self.coefs[mode][I][:, K].copy())
            if not J:
                continue
            B = [self.lookups[I][H] for H in selection['bases'][str(I)]]
            C = qm(G[I][np.ix_(B, B)]).solve(qm(G[I][np.ix_(B, J)]))
            assert qm(G[I][:, J]) == qm(G[I][:, B]) * C, ('rational pair relation', I)
            checks += len(ps) * len(J)
            lookup = {old: j for j, old in enumerate(K)}
            dest = [lookup[j] for j in B]
            new[free][I][:, dest] += from_qm(qm(self.coefs[free][I][:, J]) * C.transpose())
            blocks[str(I)] = dict(basis=[ps[j] for j in B], deleted=[ps[j] for j in J], coefficients=[[str(C[i, j]) for j in range(C.ncols())] for i in range(C.nrows())])
        self.kept = selection['kept']
        self.coefs = new
        self.reindex()
        return (blocks, checks)
