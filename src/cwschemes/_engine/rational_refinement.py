"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from itertools import product
import numpy as np
from .rational_state import ExactState
from .rational_arrays import Q, ZERO, qm, from_qm

class ExactMultiSplitState(ExactState):

    @classmethod
    def split(cls, old, c):
        previous = getattr(old, 'coordinates', [])
        assert c not in previous and old.q >= 3
        s = cls.__new__(cls)
        s.q = old.q
        s.qs = getattr(old, 'qs', [old.q] * old.k).copy()
        s.qs[c] -= 1
        s.rows = old.rows
        s.m = old.m
        s.k = old.k
        s.coordinates = previous + [c]
        bit = 1 << c
        s.digits = dict(enumerate(old.digits)) if isinstance(old.digits, list) else old.digits.copy()
        s.ordinary = dict(enumerate(old.ordinary)) if isinstance(old.ordinary, list) else old.ordinary.copy()
        s.kept = []
        origin = {}
        exceptional = set()
        nextid = max(s.digits) + 1
        for H in old.kept:
            s.kept.append(H)
            origin[H] = H
            if old.ordinary[H] & bit:
                new = nextid
                nextid += 1
                ds = list(old.digits[H])
                ds[c] = 3
                s.digits[new] = tuple(ds)
                s.ordinary[new] = old.ordinary[H] ^ bit
                s.kept.append(new)
                origin[new] = H
                exceptional.add(new)
        s.kept.sort()
        s.reindex()
        s.families = []
        q = s.q
        inv = Q(1, q)
        rinv = Q(1, q - 1)
        s.coefs = [[] for _ in range(3)]
        s.output_weights = [[] for _ in range(3)]
        for mode in range(3):
            oldfamily = old.families[mode] if previous else [(Z, 0) for Z in old.families[mode]]
            family = []
            parent = []
            for oz, (Z, E) in enumerate(oldfamily):
                for exceptional_output in [False] if Z & bit else [True, False]:
                    family.append((Z, E | (bit if exceptional_output else 0)))
                    parent.append(oz)
            s.families.append(family)
            for I, ps in enumerate(s.patterns):
                F = np.full((len(family), len(ps)), ZERO, dtype=object)
                weights = []
                J = I & ~bit
                for z, (Z, E) in enumerate(family):
                    valid = not (Z & I or E & I)
                    oz = parent[z]
                    ow = old.output_weights[mode][J][oz] if previous else Q(q) ** (s.k - s.m - J.bit_count())
                    weight = ow if Z & bit else ow * inv * (q - 1 if not (E & bit or I & bit) else 1)
                    weights.append(weight if valid else ZERO)
                    if not valid:
                        continue
                    for j, H in enumerate(ps):
                        OH = origin[H]
                        f0 = old.coefs[mode][J][oz, old.lookups[J][OH]]
                        f1 = old.coefs[mode][J | bit][oz, old.lookups[J | bit][OH]] if OH in old.lookups[J | bit] else ZERO
                        if I & bit:
                            value = f1
                        elif Z & bit or not old.ordinary[OH] & bit:
                            value = f0
                        elif H in exceptional:
                            value = f0 + (1 - inv if E & bit else -inv) * f1
                        else:
                            value = f0 + (-inv if E & bit else rinv - inv) * f1
                        F[z, j] = value
                s.coefs[mode].append(F)
                s.output_weights[mode].append(np.array(weights, dtype=object))
        assert s.rank() == old.rank()
        return s

    def grams(self, mode):
        out = []
        for I, ps in enumerate(self.patterns):
            if not ps:
                out.append(np.empty((0, 0), dtype=object))
                continue
            F = self.coefs[mode][I]
            weights = np.array([Q(self.orbit_size(self.ordinary[H] & ~I)) for H in ps], dtype=object)
            G = from_qm(qm(F).transpose() * qm(F * self.output_weights[mode][I][:, None]))
            out.append(G * weights[None, :])
        return out

    def pair_gram(self, pair):
        a = self.grams(pair[0])
        b = self.grams(pair[1])
        out = [np.full(G.shape, ZERO, dtype=object) for G in a]
        buckets = {}
        for H in self.kept:
            for K in self.kept:
                buckets.setdefault(self.ordinary[H] & self.ordinary[K], []).append((H, K))

        def transform(v, positions, inverse):
            x = v.reshape((2,) * len(positions) + (v.shape[-1],))
            for axis, c in enumerate(positions):
                a = np.moveaxis(x, axis, 0)
                b = a.copy()
                q = self.qs[c]
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
            va = transform(np.array(va, dtype=object), positions, True)
            vb = transform(np.array(vb, dtype=object), positions, True)
            result = transform(va * vb, positions, False)
            result *= np.array([Q(1, self.orbit_size(self.ordinary[K] & ~L)) for H, K in pairs], dtype=object)[None, :]
            for j, I in enumerate(Is):
                out[I][indices[j]] = result[j]
        return out
