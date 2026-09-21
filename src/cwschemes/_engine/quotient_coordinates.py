"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from .compressed_fixed_quotient import CompressedFixedQuotient, fmpq_mat
from .fixed_sector_quotient import FixedSectorQuotient
from .quotient_invariant_gram import representation

def shell(q, rows, m, groups, character):
    spec, bases, orbits, blocks = representation(q, rows, m, groups, character)
    s = CompressedFixedQuotient.__new__(CompressedFixedQuotient)
    s.q = q
    s.m = m
    s.groups = groups
    s.bases = bases
    s.domain_character = character
    permutation = list(range(3 * m))
    flips = []
    for group in groups:
        if len(group) == 1:
            flips.append(group[0])
        else:
            a, b = group
            permutation[a], permutation[b] = (b, a)
    s.base = FixedSectorQuotient(q, rows, m, permutation, (-1, -1, 1), flips)
    s.values = [[{tuple(x): int(v) for x, v in zip(b['labels'], b['vector'])} for b in bs] for bs in bases]
    s.orbit_lookup = {}
    for i, orbit in enumerate(orbits):
        s.orbit_lookup[tuple(orbit['states'])] = s.orbit_lookup[tuple(orbit['image'])] = i
    return (s, blocks, orbits)
