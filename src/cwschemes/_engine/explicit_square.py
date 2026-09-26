"""Explicit sparse rational square schemes (format ``explicit-qcsr-square-v1``).

The root object names a bundled ``fmmp.qcsr.v1`` array file: one CSR matrix per factor
(rows = terms, columns = row-major matrix entries, rational numerators/denominators).
U and V read A and B; W writes C = AB directly (same convention as the public API).
"""
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from .._resources import root as resource_root

FORMAT = 'explicit-qcsr-square-v1'


@lru_cache(None)
def _arrays(logical):
    import numpy as np
    with np.load(resource_root() / logical, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


class ExplicitSquare:
    def __init__(self, data):
        if data['format'] != FORMAT:
            raise ValueError(data['format'])
        self.n = int(data['n'])
        self.rank = int(data['rank'])
        self.source = data['source']
        self.arrays = _arrays(self.source)
        meta = json.loads(self.arrays['metadata_json'].item())
        if meta.get('format') != 'fmmp.qcsr.v1' or meta.get('coefficient_field') != 'Q':
            raise ValueError('Unsupported explicit array format')
        if meta['tensor'] != [self.n] * 3 or meta['rank'] != self.rank:
            raise ValueError('Explicit array does not match the descriptor')
        for axis in 'uvw':
            ptr = self.arrays[axis + '_indptr']
            if len(ptr) != self.rank + 1 or int(ptr[0]) != 0:
                raise ValueError('Malformed CSR pointer')
            if (ptr[1:] - ptr[:-1] <= 0).any():
                raise ValueError('Zero factor in an explicit scheme')
            if (self.arrays[axis + '_denominators'] <= 0).any():
                raise ValueError('Nonpositive denominator')

    def factor_entry(self, mode, term, row, col):
        import numpy as np
        if mode not in (0, 1, 2) or not (0 <= term < self.rank and 0 <= row < self.n and 0 <= col < self.n):
            raise IndexError((mode, term, row, col))
        a = self.arrays
        axis = 'uvw'[mode]
        start, end = int(a[axis + '_indptr'][term]), int(a[axis + '_indptr'][term + 1])
        idx = a[axis + '_indices']
        variable = self.n * row + col
        k = start + int(np.searchsorted(idx[start:end], variable))
        if k == end or idx[k] != variable:
            return F(0)
        return F(int(a[axis + '_numerators'][k]), int(a[axis + '_denominators'][k]))

    # ------------------------------------------------------------------ checks
    def _dense_mod(self, axis, p):
        """Factor matrix (rank x n^2) reduced modulo p as int64 in [0, p)."""
        import numpy as np
        a = self.arrays
        ptr = a[axis + '_indptr']
        rows = np.repeat(np.arange(self.rank), np.diff(ptr))
        den = a[axis + '_denominators']
        inverse = {int(d): pow(int(d) % p, -1, p) for d in np.unique(den)}
        inv = np.array([inverse[int(d)] for d in den], dtype=np.int64)
        val = (np.mod(a[axis + '_numerators'], p) * inv) % p
        M = np.zeros((self.rank, self.n * self.n), dtype=np.int64)
        np.add.at(M, (rows, a[axis + '_indices']), val)
        return M % p

    def random_check(self, trials=2, seed=0, p=2147483647):
        """Evaluate the bilinear algorithm on random matrices modulo a large prime and compare with A@B.

        A nonzero error polynomial of degree 2 survives one trial with probability <= 2/p (Schwartz–Zippel)."""
        import numpy as np
        rng = np.random.default_rng(seed)
        n = self.n
        U, V, W = (self._dense_mod(axis, p) for axis in 'uvw')
        for _ in range(trials):
            A = rng.integers(0, p, size=n * n, dtype=np.int64)
            B = rng.integers(0, p, size=n * n, dtype=np.int64)
            a = _matvec_mod(U, A, p)
            b = _matvec_mod(V, B, p)
            m = (a * b) % p
            C = _matvec_mod(W.T, m, p)
            ref = _matmul_mod(A.reshape(n, n), B.reshape(n, n), p).reshape(-1)
            if not np.array_equal(C, ref):
                return False
        return True

    def exact_check(self, block=None):
        """Check every tensor coefficient modulo primes whose product exceeds a rigorous bound; exact over Q.

        Delta = lcm(U dens)*lcm(V dens)*lcm(W dens) clears all denominators, and each coefficient of
        Delta*(sum_t u_t v_t w_t - T) is an integer bounded by Delta*(sum_t |u_t|_inf |v_t|_inf |w_t|_inf + 1)."""
        import numpy as np
        a = self.arrays
        n, n2, R = self.n, self.n * self.n, self.rank
        S = 0.0
        lcm = 1
        for axis in 'uvw':
            ptr = a[axis + '_indptr']
            vals = np.abs(a[axis + '_numerators'].astype(np.float64)) / a[axis + '_denominators']
            mx = np.maximum.reduceat(vals, ptr[:-1])
            S = mx if isinstance(S, float) else S * mx
            L = 1
            for d in np.unique(a[axis + '_denominators']):
                L = L * int(d) // math.gcd(L, int(d))
            lcm *= L
        need = 2 * lcm * (float(np.sum(S)) * (1 + 1e-9) + 1)
        # primes p with p^2 * R < 2^53 so that one float64 GEMM is exact
        limit = int(math.isqrt(int(2 ** 53 // max(R, 1))))
        primes, prod, q = [], 1, limit
        while prod <= need:
            q -= 1
            if q < 3:
                raise RuntimeError('Not enough primes for the CRT bound')
            if all(q % k for k in range(2, math.isqrt(q) + 1)) and lcm % q:
                primes.append(q)
                prod *= q
        block = block or max(1, min(n2, (1 << 25) // max(n2 * R, 1)))
        for p in primes:
            U, V, W = (self._dense_mod(axis, p).astype(np.float64) for axis in 'uvw')
            Wt = W  # (R, n2)
            for a0 in range(0, n2, block):
                a1 = min(n2, a0 + block)
                KR = (U[:, a0:a1].T[:, None, :] * V.T[None, :, :]) % p   # (blk, n2, R)
                P = np.mod(KR.reshape(-1, R) @ Wt, p).reshape(a1 - a0, n2, n2)
                for t, x in enumerate(range(a0, a1)):
                    i, j = divmod(x, n)
                    expected = np.zeros((n2, n2))
                    expected[j * n + np.arange(n), i * n + np.arange(n)] = 1
                    if not np.array_equal(P[t], expected):
                        return dict(exact=False, prime=p)
        return dict(exact=True, primes=primes, delta_bits=round(math.log2(lcm), 1), bound_bits=round(math.log2(need), 1))


def _matvec_mod(M, x, p):
    import numpy as np
    # split x into 16-bit limbs so that int64 products never overflow
    lo = x & 0xFFFF
    hi = x >> 16
    r_lo = np.zeros(M.shape[0], dtype=np.int64)
    r_hi = np.zeros(M.shape[0], dtype=np.int64)
    step = 1 << 14
    for c in range(0, M.shape[1], step):
        Mc = M[:, c:c + step]
        r_lo = (r_lo + (Mc * lo[c:c + step]).sum(axis=1) % p) % p
        r_hi = (r_hi + (Mc * hi[c:c + step]).sum(axis=1) % p) % p
    return (r_lo + (r_hi * 65536) % p) % p


def _matmul_mod(A, B, p):
    import numpy as np
    n = A.shape[0]
    C = np.zeros((n, n), dtype=np.int64)
    for k in range(n):
        C = (C + (A[:, k:k + 1] * (B[k:k + 1, :] & 0xFFFF)) % p + ((A[:, k:k + 1] * (B[k:k + 1, :] >> 16)) % p * 65536) % p) % p
    return C


def load_certificate(data):
    return ExplicitSquare(data)
