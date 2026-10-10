"""Replay an obstruction to lifting the frozen 4AT binary factors modulo 4.

This script does not import a solver, a saved Jacobian, or a rank report.
It independently rebuilds the tensor and all derivative contractions from
the hash-bound direct-output binary factors and a list of witness rows.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


SOURCE_SHA256 = 'fe420dbd062c3da4ad8d76185ae304dc1abd480744d0692c91c6c35f7264fc6f'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_factors(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    require(digest == SOURCE_SHA256, 'base NPZ does not match the frozen source hash')
    factors = np.zeros((3, 47, 16), dtype=np.int64)
    with np.load(path, allow_pickle=False) as source:
        metadata = json.loads(str(source['metadata_json']))
        require(metadata['tensor'] == [4, 4, 4], 'wrong source dimensions')
        require(metadata['rank'] == 47 and metadata['coefficient_field'] == 'F2',
                'wrong source rank or coefficient field')
        for mode, name in enumerate('uvw'):
            ptr = source[name + '_indptr']
            columns = source[name + '_indices']
            numerators = source[name + '_numerators']
            denominators = source[name + '_denominators']
            require(len(ptr) == 48 and ptr[0] == 0 and ptr[-1] == len(columns),
                    'invalid CSR row pointers')
            require(np.all(np.diff(ptr) >= 0), 'nonmonotone CSR row pointers')
            require(len(columns) == len(numerators) == len(denominators),
                    'inconsistent CSR coefficient lengths')
            require(np.all((columns >= 0) & (columns < 16)), 'invalid factor coordinate')
            require(np.all(denominators % 2 == 1), 'factor has an even denominator')
            for term in range(47):
                for entry in range(int(ptr[term]), int(ptr[term + 1])):
                    factors[mode, term, int(columns[entry])] ^= int(numerators[entry]) % 2
    return factors, digest


def main():
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=directory / 'bases' / '4at.npz')
    parser.add_argument('--witness', type=Path, default=directory / '4at_mod4_obstruction.json')
    args = parser.parse_args()
    factors, source_hash = load_factors(args.base)
    certificate = json.loads(args.witness.read_text())
    require(certificate['source_sha256'] == source_hash, 'witness refers to a different base')
    require(certificate['convention'] == 'direct output row-major', 'wrong witness convention')
    rows = certificate['coefficient_rows']
    require(all(type(row) is int and 0 <= row < 4096 for row in rows), 'invalid witness row')
    require(len(rows) == len(set(rows)) == 520, 'the witness must contain 520 distinct rows')
    require(certificate['weight'] == len(rows), 'incorrect witness weight')
    if 'tensor_coordinates' in certificate:
        coordinates = [[row // 256, (row // 16) % 16, row % 16] for row in rows]
        require(certificate['tensor_coordinates'] == coordinates, 'incorrect witness coordinates')
    witness = np.zeros(4096, dtype=np.int64)
    witness[rows] = 1
    witness = witness.reshape(16, 16, 16)

    # A[i,j] B[j,k] contributes to direct output C[i,k].
    expected = np.zeros((16, 16, 16), dtype=np.int64)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                expected[i * 4 + j, j * 4 + k, i * 4 + k] = 1
    integer_tensor = np.zeros_like(expected)
    for term in range(47):
        integer_tensor += np.multiply.outer(np.multiply.outer(factors[0, term],
                                                              factors[1, term]),
                                           factors[2, term])
    difference = expected - integer_tensor
    require(not np.any(difference % 2), 'base fails an exact F2 tensor coefficient')
    residual = difference // 2

    # These are all 3 * 47 * 16 = 2256 derivative-column contractions.
    a, b, c = factors
    derivatives = (
        np.einsum('ijk,rj,rk->ri', witness, b, c),
        np.einsum('ijk,ri,rk->rj', witness, a, c),
        np.einsum('ijk,ri,rj->rk', witness, a, b),
    )
    for mode, contractions in enumerate(derivatives):
        require(not np.any(contractions % 2), f'witness does not annihilate mode {mode} derivatives')
    pairing = int(np.sum(witness * residual) % 2)
    require(pairing == 1, 'witness does not obstruct the first Hensel equation')
    print(json.dumps(dict(
        verified=True, source_sha256=source_hash, source_rank=47,
        source_tensor_coefficients_verified=4096, witness_weight=len(rows),
        derivative_columns_annihilated=2256, residual_pairing=pairing,
        conclusion='No rank-47 lift modulo 4 with these binary factor residues.',
        scope='Also excludes rational odd-denominator factor lifts reducing to this seed; no general rational, real, or complex rank lower bound.'
    ), indent=2))


if __name__ == '__main__':
    main()
