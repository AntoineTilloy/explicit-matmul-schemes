"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from fractions import Fraction as F
import json
import math
import numpy as np
from .compressed_fixed_quotient import CompressedFixedQuotient
from .quotient_hadamard_engine import HadamardEngine, update_free, mod_fraction
from .quotient_invariant_gram import standard_problem

def coefficient_bound(proof_path):
    adapter = CompressedFixedQuotient(proof_path, allow_stage=True)
    weights = [weight for records in adapter.transforms.values() for _, _, weight in records]
    denominator = math.lcm(*(w.denominator for w in weights))
    paired = adapter.q - adapter.q % 2
    moment_denominator = math.prod((math.prod(range(adapter.q - 3, adapter.q + 1)) if len(g) == 2 else paired * (paired - 2) for g in adapter.groups))
    group_order = math.prod((math.factorial(adapter.q) if len(g) == 2 else 2 ** (adapter.q // 2) * math.factorial(adapter.q // 2) for g in adapter.groups))
    assert group_order % moment_denominator == 0
    total = 16 * sum(map(abs, weights), F(0))
    bound = (total.numerator + total.denominator - 1) // total.denominator
    return (denominator * moment_denominator, bound)

def residue(q, prime, history, pair, ids, denominator, threads):
    rows, m, groups = standard_problem()
    engine = HadamardEngine(q, rows, m, groups, prime=prime, threads=threads)
    kernels = engine.initial_kernels()
    for path in history:
        proof = json.loads(path.read_text())
        pp = tuple(proof['pair'])
        character = math.prod(((-1, -1, 1)[i] for i in pp))
        metadata = engine.quotients[character]
        relations = {}
        for key, r in proof['relations'].items():
            C = np.array([[mod_fraction(F(value), prime) for value in row] for row in r['coefficients']], dtype=np.int64)
            if proof.get('completion_sign_omitted', True) and 2 in pp:
                records = metadata[int(key)]['records']

                def sign(index):
                    return (-1) ** sum((int(x == 0) for bs, j in zip(engine.bases, records[index]['indices']) for x in bs[j]['labels'][0]))
                C = C * np.array([sign(i) for i in r['basis']])[:, None] * np.array([sign(i) for i in r['deleted']])[None, :] % prime
            relations[int(key)] = dict(basis=r['basis'], deleted=r['deleted'], coefficients=C)
        update_free(engine, kernels, 3 - sum(pp), character, metadata, relations, proof['deleted_orbits'])
    character = math.prod(((-1, -1, 1)[i] for i in pair))
    matrices, metadata = engine.kernel_to_blocks(kernels[pair[0]] * kernels[pair[1]] % prime, character)
    scale = denominator % prime
    return [matrices[i] * scale % prime for i in ids]
