#!/usr/bin/env python3
"""Source-only Atlas gate: exact grading or fractional-obstruction certificates.

HiGHS proposes solutions; only complete rational/integer replays certify them.
New objective has no size cap. No finite sub-2.65 record is asserted here.
"""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import reduce
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import time
import warnings

EXPECTED = '71d736f4152a58c8769a4b09208454119fd28b48ebba2c4dd66b23a64b48e04d'
START = time.process_time()
INTERNAL_CAP = 8.5


def guard():
    if time.process_time() - START > INTERNAL_CAP:
        raise TimeoutError('internal 8.5 CPU-second cap')


def masks(rows):
    return [[sum(1 << i for i, ch in enumerate(row) if ch == str(h + 1))
             for row in rows] for h in range(3)]


def sunflower_edges(mode_masks):
    size = len(mode_masks[0])
    out = []
    for i, j, k in itertools.product(range(size), repeat=3):
        a, b, c = mode_masks[0][i], mode_masks[1][j], mode_masks[2][k]
        if a & b == a & c == b & c:
            out.append((i, j, k))
        guard()
    return out


def grade(edge, potential):
    return sum(potential[h][edge[h]] for h in range(3))


def potential_certificate(rows, edges, potential):
    size = len(rows)
    if any(grade((i, i, i), potential) != 0 for i in range(size)):
        return None
    bad = [e for e in edges if not e[0] == e[1] == e[2]]
    values = [grade(e, potential) for e in bad]
    if values and min(values) <= 0:
        return None
    assert all(isinstance(v, int) for mode in potential for v in mode)
    common = reduce(math.gcd, (abs(v) for mode in potential for v in mode), 0)
    if common > 1:
        potential = [[v // common for v in mode] for mode in potential]
        values = [v // common for v in values]
    return dict(kind='INTEGER_POSITIVE_GRADING', potentials=potential,
                degree=max(values, default=0), minimum_unwanted_grade=min(values, default=0),
                unwanted_grade_histogram=dict(sorted(Counter(values).items())),
                desired_edges_checked=size, unwanted_edges_checked=len(bad),
                exact_replay='PASS', degree_optimality_claim=False)


def recover_potential(rows, edges, vector):
    size = len(rows)
    if not all(math.isfinite(float(v)) for v in vector):
        return None
    candidates = []
    # Desired grades are eliminated algebraically: c_i=-a_i-b_i.
    # A grid recovery is accepted only after complete exact edge verification.
    for denominator in (1, 2, 4, 8, 16, 64, 256, 4096):
        guard()
        aa = [int(round(float(v) * denominator)) for v in vector[:size]]
        bb = [int(round(float(v) * denominator)) for v in vector[size:2 * size]]
        cc = [-a - b for a, b in zip(aa, bb)]
        certificate = potential_certificate(rows, edges, [aa, bb, cc])
        if certificate:
            certificate['recovery_grid_denominator'] = denominator
            candidates.append(certificate)
    return min(candidates, key=lambda c: c['degree']) if candidates else None


def flow_certificate(size, edges, weights):
    assert len(weights) == len(edges)
    if any(v < 0 for v in weights):
        return None
    sums = [[F(0) for _ in range(size)] for _ in range(3)]
    unwanted = F(0)
    support = []
    for edge, value in zip(edges, weights):
        if not value:
            continue
        for h in range(3):
            sums[h][edge[h]] += value
        if not edge[0] == edge[1] == edge[2]:
            unwanted += value
        support.append(dict(edge=edge, numerator=value.numerator, denominator=value.denominator))
    if not all(v == 1 for mode in sums for v in mode) or unwanted <= 0:
        return None
    return dict(kind='RATIONAL_FRACTIONAL_MATCHING_OBSTRUCTION', support=support,
                unwanted_mass=[unwanted.numerator, unwanted.denominator],
                all_vertex_marginals_one=True, exact_replay='PASS')


def exact_support_solve(size, edges, approximate):
    """Exact rational recovery on the numerical basic support; no LP claims."""
    selected = [i for i, x in enumerate(approximate) if float(x) > 1e-9]
    if len(selected) > 3 * size + 5:
        return None
    count = len(selected)
    matrix = [[F(int(edges[index][h] == v)) for index in selected] + [F(1)]
              for h in range(3) for v in range(size)]
    pivots = []
    at = 0
    for col in range(count):
        guard()
        pivot = next((r for r in range(at, len(matrix)) if matrix[r][col]), None)
        if pivot is None:
            continue
        matrix[at], matrix[pivot] = matrix[pivot], matrix[at]
        divisor = matrix[at][col]
        matrix[at] = [x / divisor for x in matrix[at]]
        for r in range(len(matrix)):
            if r != at and matrix[r][col]:
                coefficient = matrix[r][col]
                matrix[r] = [a - coefficient * b for a, b in zip(matrix[r], matrix[at])]
        pivots.append(col)
        at += 1
    if any(not any(row[:count]) and row[count] for row in matrix):
        return None
    free = set(range(count)) - set(pivots)
    x = [F(0) for _ in range(count)]
    for col in free:
        x[col] = F(float(approximate[selected[col]])).limit_denominator(65536)
    for r, col in enumerate(pivots):
        x[col] = matrix[r][count] - sum(matrix[r][j] * x[j] for j in free)
    full = [F(0) for _ in edges]
    for index, value in zip(selected, x):
        full[index] = value
    return flow_certificate(size, edges, full)


def recover_flow(size, edges, approximate):
    if not all(math.isfinite(float(x)) for x in approximate):
        return None
    for denominator in (16, 256, 4096, 65536):
        guard()
        values = [F(float(x)).limit_denominator(denominator) for x in approximate]
        certificate = flow_certificate(size, edges, values)
        if certificate:
            return certificate
    return exact_support_solve(size, edges, approximate)


def local_G(q):
    forms = [[F(1)] + [F(i == j) + F(1, q) for j in range(q)] for i in range(q)]
    forms += [[F(1)] + [F(3, 2 * q)] * q]
    weights = [F(1)] * q + [F(-4 * q, 3)]
    table = {}
    for x, y, z in itertools.product(range(q + 1), repeat=3):
        table[x, y, z] = sum(c * u[x] * u[y] * u[z] for u, c in zip(forms, weights))
        zero_count = (x == 0) + (y == 0) + (z == 0)
        expected = int((x == 0 and y == z != 0) or (y == 0 and x == z != 0)
                       or (z == 0 and x == y != 0))
        if zero_count in (1, 2):
            assert table[x, y, z] == expected
    assert table[0, 0, 0] == -F(q, 3)
    assert table[1, 1, 1] == 1 + F(3, q) - F(1, 2 * q * q) > 0
    return table


def tiny_word_control(rows, potential, q=2):
    width = len(rows[0])
    mm = masks(rows)
    variables = []
    for h in range(3):
        mode = []
        for owner, mask in enumerate(mm[h]):
            positions = [i for i in range(width) if not mask >> i & 1]
            for labels in itertools.product(range(1, q + 1), repeat=len(positions)):
                word = [0] * width
                for i, value in zip(positions, labels):
                    word[i] = value
                mode.append((owner, word))
        variables.append(mode)
    table = local_G(q)
    checked = 0
    for left, middle, right in itertools.product(*variables):
        owners = (left[0], middle[0], right[0])
        words = (left[1], middle[1], right[1])
        value = math.prod(table[tuple(word[i] for word in words)] for i in range(width))
        g = grade(owners, potential)
        target = int(owners[0] == owners[1] == owners[2]
                     and all(sum(word[i] == 0 for word in words) == 1
                             and len({word[i] for word in words if word[i]}) == 1
                             for i in range(width)))
        if value:
            assert g >= 0
        assert (value if g == 0 else F(0)) == target
        checked += 1
        guard()
    return dict(q=q, dimensions=list(map(len, variables)), coefficients_checked=checked,
                constant_coefficient_identity='PASS')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--puzzles', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    a = parser.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    raw = a.puzzles.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == EXPECTED
    assignments = [n for n in ast.parse(raw.decode()).body if isinstance(n, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == 'P' for t in n.targets)]
    assert len(assignments) == 1
    puzzles = ast.literal_eval(assignments[0].value)
    result = dict(status='RUNNING', new_record=False, size_cap=None,
                  bindings=dict(puzzles_sha256=EXPECTED,
                                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
                  targets=[], controls={})
    def save():
        temporary = a.out / 'result.tmp'
        temporary.write_text(json.dumps(result, indent=2) + '\n')
        temporary.replace(a.out / 'result.json')
    save()
    try:
        # Local identities hold over Q; actual tiny restricted tensor replay
        # checks all zero coefficients as well as wanted physical coefficients.
        for q in (2, 3):
            local_G(q)
        result['controls']['local_G_q2_q3'] = 'PASS'
        for key, potential in [((2, 2), [[0, 2], [0, -1], [0, -1]]),
                               ((3, 3), [[0, 3, 4], [0, -1, -2], [0, -2, -2]])]:
            edges = sunflower_edges(masks(puzzles[key]))
            certificate = potential_certificate(puzzles[key], edges, potential)
            assert certificate
            result['controls'][f'word_{key[0]}_{key[1]}'] = tiny_word_control(puzzles[key], potential)
        # Fractional obstruction control has exact marginals and unwanted mass.
        control_edges = [(0, 0, 0), (1, 1, 1), (0, 1, 0), (1, 0, 1)]
        assert flow_certificate(2, control_edges, [F(0), F(0), F(1), F(1)])
        result['controls']['fractional_obstruction_replay'] = 'PASS'
        save()
        for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            os.environ[name] = '1'
        import numpy as np
        import scipy
        from scipy.optimize import linprog, OptimizeWarning
        from scipy.sparse import coo_matrix
        result['scipy_version'] = scipy.__version__
        result['numpy_version'] = np.__version__
        save()

        def lp(*args, **kwargs):
            guard()
            remaining = INTERNAL_CAP - (time.process_time() - START)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=OptimizeWarning,
                                        message='Unrecognized options detected.*')
                return linprog(*args, method='highs',
                               options=dict(presolve=True, threads=1, time_limit=max(0.01, min(1.5, remaining))),
                               **kwargs)

        for key in ((23, 7), (35, 8), (14, 6), (8, 5), (5, 4)):
            guard()
            rows = puzzles[key]
            size, width = key
            mm = masks(rows)
            item = dict(puzzle=list(key), rows=rows, masks=mm, status='BUILDING_SUPPORT')
            result['targets'].append(item)
            save()
            if not all(len(set(mode)) == size for mode in mm):
                item['status'] = 'DUPLICATE_MODE_MASKS'
                save()
                continue
            edges = sunflower_edges(mm)
            bad = [e for e in edges if not e[0] == e[1] == e[2]]
            assert all((i, i, i) in edges for i in range(size))
            edge_bytes = (json.dumps(edges) + '\n').encode()
            (a.out / f'edges_{size}_{width}.json').write_bytes(edge_bytes)
            item.update(mode_injective=True, physical_structural_edges=len(edges),
                        unwanted_edges=len(bad), edges_sha256=hashlib.sha256(edge_bytes).hexdigest(),
                        status='GRADING_LP_RUNNING')
            save()
            # Variables a_i,b_i,D; c_i=-a_i-b_i gives desired grade exactly0.
            row_ids, columns, entries, rhs = [], [], [], []
            for index, (i, j, k) in enumerate(bad):
                coefficients = Counter()
                coefficients[i] += 1
                coefficients[size + j] += 1
                coefficients[k] -= 1
                coefficients[size + k] -= 1
                for col, value in coefficients.items():
                    if value:
                        row_ids.extend((2 * index, 2 * index + 1))
                        columns.extend((col, col))
                        entries.extend((-value, value))
                row_ids.append(2 * index + 1)
                columns.append(2 * size)
                entries.append(-1)
                rhs.extend((-1, 0))
            matrix = coo_matrix((entries, (row_ids, columns)), shape=(2 * len(bad), 2 * size + 1)).tocsr()
            objective = np.zeros(2 * size + 1)
            objective[-1] = 1
            bounds = [(None, None)] * (2 * size) + [(1, None)]
            bounds[0] = bounds[size] = (0, 0)
            proposal = lp(objective, A_ub=matrix, b_ub=np.array(rhs), bounds=bounds)
            item['grading_lp_status'] = int(proposal.status)
            item['grading_lp_message'] = str(proposal.message)
            if proposal.x is not None:
                item['proposed_degree_float'] = float(proposal.x[-1])
                certificate = recover_potential(rows, edges, proposal.x)
                if certificate:
                    item.update(status='CERTIFIED_POSITIVE_GRADING', certificate=certificate)
                    save()
                    continue
            # A numerical infeasibility flag alone is not a certificate.
            # Seek a rational fractional perfect matching using unwanted edges.
            item['status'] = 'OBSTRUCTION_LP_RUNNING'
            save()
            rr, cc = [], []
            for index, edge in enumerate(edges):
                for h in range(3):
                    rr.append(h * size + edge[h])
                    cc.append(index)
            incidence = coo_matrix((np.ones(len(rr)), (rr, cc)), shape=(3 * size, len(edges))).tocsr()
            objective = np.array([-int(not e[0] == e[1] == e[2]) for e in edges], dtype=float)
            proposal = lp(objective, A_eq=incidence, b_eq=np.ones(3 * size), bounds=(0, None))
            item['obstruction_lp_status'] = int(proposal.status)
            certificate = recover_flow(size, edges, proposal.x) if proposal.x is not None else None
            if certificate:
                item.update(status='CERTIFIED_NO_POSITIVE_GRADING', certificate=certificate)
            else:
                item['status'] = 'UNKNOWN_NO_EXACT_CERTIFICATE'
            save()
        result['status'] = 'COMPLETE'
    except TimeoutError as error:
        result.update(status='UNKNOWN_CPU_CAP', reason=str(error))
    finally:
        result['CPU_seconds'] = time.process_time() - START
        save()


if __name__ == '__main__':
    main()
