"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from itertools import combinations, product
from math import prod
from .restrictions import audit, completion, variables, rank_bound, SEED

def plan_group(rows, m, s):
    k = 3 * m
    families = [{r[j] for r in rows if r[j] & s == s} for j in range(3)]
    remaining = set((i for i in range(k) if not s >> i & 1))
    patches = []
    while True:
        found = False
        for i in sorted(remaining):
            bit = 1 << i
            zero = [all((a & bit for a in f)) for f in families]
            positive = [all((not a & bit for a in f)) for f in families]
            if any(zero) and any(positive):
                j = zero.index(True)
                for h in range(3):
                    if h != j:
                        families[h] = {a for a in families[h] if not a & bit}
                patches.append(dict(kind='identity', positions=[i], zero_mode=j))
                remaining.remove(i)
                if not all(families):
                    return dict(pure=s, zero=True, patches=[], remaining=[], families=[[]] * 3)
                found = True
                break
        if not found:
            break
    for i, j in combinations(sorted(remaining), 2):
        if i not in remaining or j not in remaining:
            continue
        bits = 1 << i | 1 << j
        patterns = [{a & bits for a in f} for f in families]
        for h in range(3):
            other = [t for t in range(3) if t != h]
            if patterns[h] == {0} and all((patterns[t] == {1 << i, 1 << j} for t in other)):
                patches.append(dict(kind='two_matvec', positions=[i, j], matrix_mode=h))
                remaining -= {i, j}
                break
    for i, j in combinations(sorted(remaining), 2):
        if i not in remaining or j not in remaining:
            continue
        bits = 1 << i | 1 << j
        patterns = [{a & bits for a in f} for f in families]
        pairs = [(a, b) for a in range(3) for b in range(3) if a != b and patterns[a] == {0} and (0 not in patterns[b])]
        if pairs:
            patches.append(dict(kind='pair_null', positions=[i, j], pair=list(pairs[0])))
            remaining -= {i, j}
    return dict(pure=s, zero=False, patches=patches, remaining=sorted(remaining), families=[sorted(f) for f in families])

def patch_cost(q, patch):
    return q if patch['kind'] == 'identity' else q * q + 2 * q

def group_cost(q, group):
    if group['zero']:
        return 0
    return (q + 1) ** len(group['remaining']) * prod((patch_cost(q, p) for p in group['patches']))

def packet_plan(rows=SEED, m=2, q=16):
    witness = audit(rows, m)
    groups = [plan_group(rows, m, s) for s in witness['surviving_pure_subsets']]
    return dict(witness=witness, q=q, n=q ** m, original_rank=rank_bound(q, witness), rank=sum((group_cost(q, g) for g in groups)), groups=groups)
