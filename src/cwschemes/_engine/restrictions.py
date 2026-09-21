"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from itertools import combinations, product
from fractions import Fraction
from math import prod
SEED = [(3, 12, 48), (33, 10, 20), (6, 24, 33), (34, 5, 24), (12, 17, 34), (40, 20, 3)]

def masks(k, m):
    return [sum((1 << i for i in s)) for s in combinations(range(k), m)]

def audit(rows, m):
    """Check the induced matching and enumerate surviving pure-term subsets."""
    rows = [tuple(r) for r in rows]
    k = 3 * m
    full = (1 << k) - 1
    assert rows and len(set(rows)) == len(rows)
    assert all((all((0 <= a <= full and a.bit_count() == m for a in r)) and r[0] | r[1] | r[2] == full and (not (r[0] & r[1] or r[0] & r[2] or r[1] & r[2])) for r in rows))
    sectors = [set((r[j] for r in rows)) for j in range(3)]
    assert all((len(s) == len(rows) for s in sectors)), 'shared input/output sector'
    edges = {(a, b, full ^ (a | b)) for a in sectors[0] for b in sectors[1] if not a & b and full ^ (a | b) in sectors[2]}
    assert edges == set(rows), f'{len(edges - set(rows))} cross terms'
    downsets = []
    for mode in sectors:
        down = set()
        for a in mode:
            s = a
            while True:
                down.add(s)
                if not s:
                    break
                s = s - 1 & a
        downsets.append(down)
    common = set.intersection(*downsets)
    surviving = [s for d in range(m + 1) for s in masks(k, d) if s in common]
    counts = [sum((s.bit_count() == d for s in surviving)) for d in range(m + 1)]
    return dict(m=m, k=k, products=len(rows), rows=rows, survivor_counts=counts, surviving_pure_subsets=surviving)

def rank_bound(q, witness):
    return sum((c * (q + 1) ** (witness['k'] - d) for d, c in enumerate(witness['survivor_counts'])))

def completion(q):
    """Rank-q+2 exact completion over Q; weights are absorbed in mode Z."""
    forms = [[Fraction(1)] + [Fraction(i == j) + Fraction(1, q) for j in range(q)] for i in range(q)]
    forms += [[Fraction(1)] + [Fraction(3, 2 * q)] * q, [Fraction(1)] + [Fraction(0)] * q]
    weights = [Fraction(1)] * q + [Fraction(-4 * q, 3), Fraction(q, 3)]
    return (forms, weights)

def verify_completion(q):
    f, w = completion(q)
    checked = 0
    for a, b, c in product(range(q + 1), repeat=3):
        if a and b and c:
            continue
        actual = sum((ws * u[a] * u[b] * u[c] for u, ws in zip(f, w)))
        expected = int(a == 0 and b == c != 0 or (b == 0 and a == c != 0) or (c == 0 and a == b != 0))
        assert actual == expected, (q, a, b, c, actual, expected)
        checked += 1
    return checked

def variables(sectors, k, q):
    return [tuple((next(it) if not mask & 1 << i else 0 for i in range(k))) for mask in sorted(sectors) for labels in product(range(1, q + 1), repeat=k - mask.bit_count()) for it in [iter(labels)]]

def expand_mod(q, rows, m, p=101):
    """Materialize a SMALL instance only. No dense real/border approximation."""
    import numpy as np
    w = audit(rows, m)
    f, weights = completion(q)
    mod = lambda x: x.numerator * pow(x.denominator, -1, p) % p
    f = [[mod(x) for x in u] for u in f]
    weights = list(map(mod, weights))
    vs = [variables({r[j] for r in rows}, 3 * m, q) for j in range(3)]
    terms = [t for t in product(range(q + 2), repeat=3 * m) if sum((1 << i for i, s in enumerate(t) if s == q + 1)) in w['surviving_pure_subsets']]
    assert len(terms) == rank_bound(q, w)
    factors = [np.array([[prod((f[s][v[i]] for i, s in enumerate(t))) % p for t in terms] for v in mode], dtype=np.int64) for mode in vs]
    factors[2] = factors[2] * np.array([prod((weights[s] for s in t)) % p for t in terms]) % p
    return (vs, factors)

def verify_expanded(q=2, rows=SEED, m=2, p=101):
    """Exhaustive all-coefficient check, including every zero coefficient."""
    import numpy as np
    vs, (u, v, w) = expand_mod(q, rows, m, p)
    assert max(map(len, vs)) * len(u[0]) * p ** 3 < 2 ** 63
    for i, x in enumerate(vs[0]):
        got = v * u[i] % p @ w.T % p
        target = np.array([[all((a == 0 and b == c != 0 or (b == 0 and a == c != 0) or (c == 0 and a == b != 0) for a, b, c in zip(x, y, z))) for z in vs[2]] for y in vs[1]], dtype=np.int64)
        assert np.array_equal(got, target), f'coefficient slice {i} failed'
    return dict(q=q, prime=p, dimensions=list(map(len, vs)), rank=u.shape[1], coefficients=prod(map(len, vs)))
