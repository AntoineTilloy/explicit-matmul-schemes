"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
def allowed(rows, k, pair):
    full = (1 << k) - 1
    patterns = set()
    for a in {r[pair[0]] for r in rows}:
        for b in {r[pair[1]] for r in rows}:
            pure = a & b
            required = a ^ b
            opt = full ^ (a | b)
            s = opt
            while True:
                ordinary = required | s
                code = 0
                for i in range(k):
                    code = 3 * code + (0 if pure >> i & 1 else 1 if ordinary >> i & 1 else 2)
                patterns.add(code)
                if not s:
                    break
                s = s - 1 & opt
    return patterns
