"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
import math
from .quotient_invariant_gram import representation, local_tables, assemble
PRIMES = (1000000007, 1000000009, 998244353, 1004535809)

def reconstruct(metadata, block_ids, threads=1):
    from numba import set_num_threads
    set_num_threads(threads)
    q, rows, m, groups, pair = (metadata[k] for k in ('q', 'rows', 'm', 'groups', 'pair'))
    spec, bases, orbits, blocks = representation(q, rows, m, groups, math.prod(((-1, -1, 1)[i] for i in pair)))
    blocks = [blocks[i] for i in block_ids]
    bound = 16 * len(rows) ** 2 * q ** (10 * m)
    primes = list(PRIMES)
    if math.prod(primes) <= 2 * bound:
        from sympy import nextprime
        while math.prod(primes) <= 2 * bound:
            primes.append(int(nextprime(max(primes))))
    values = None
    modulus = 1
    for prime in primes:
        tables, _ = local_tables(q, spec, bases, groups, tuple(pair), prime)
        matrices = assemble(tables, blocks, tuple(pair), prime)
        if values is None:
            values = [matrix.astype(object) for matrix in matrices]
        else:
            inverse = pow(modulus, -1, prime)
            for value, matrix in zip(values, matrices):
                delta = (matrix - value) % prime * inverse % prime
                value += modulus * delta
        modulus *= prime
        print('CRT modulus bits', modulus.bit_length(), flush=True)
    for value in values:
        value[value > modulus // 2] -= modulus
        assert all((abs(int(x)) <= bound for x in value.flat))
    return (values, bound, modulus, primes)
