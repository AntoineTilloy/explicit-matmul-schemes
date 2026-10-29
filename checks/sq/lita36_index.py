#!/usr/bin/env python3
"""The indexed current-LITA accessor at size 36 agrees with the pinned serializer, term by term.

The outer scheme of the size-13,947,137,604 entry is CurrentLITA(36). The independent
exact CRT replay (checks/sq/independent.json) certifies the serializer lita(36) of the
same pinned generator; this check binds the two. Usage, from the repository root:
    PYTHONPATH=src python checks/sq/lita36_index.py [--out checks/sq/lita36_index.json]
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

from cwschemes._engine.current_lita import CurrentLITA
from cwschemes._engine.sources import lita_current as L


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    n = 36
    digest = hashlib.sha256(Path(L.__file__).read_bytes()).hexdigest()
    compact, expected = CurrentLITA(n, digest), L.lita(n)
    assert compact.rank == expected.rank == L.lita_rank(n) == 19666
    arrays = {}
    for name in 'UVW':
        arrays.update(getattr(expected, name).arrays(name.lower()))
    entries = 0
    for t in range(compact.rank):
        forms = compact.term(t)
        for mode, axis in enumerate('uvw'):
            ptr = arrays[axis + '_indptr']
            want = {}
            for j in range(int(ptr[t]), int(ptr[t + 1])):
                coordinate = divmod(int(arrays[axis + '_indices'][j]), n)
                want[coordinate] = Fraction(int(arrays[axis + '_numerators'][j]), int(arrays[axis + '_denominators'][j]))
            assert forms[mode] == want, (t, axis)
            entries += len(want)
    record = dict(status='PASS', n=n, rank=compact.rank, terms_compared=compact.rank, nonzero_entries=entries,
                  generator_sha256=digest, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  seconds=time.monotonic() - start)
    if args.out:
        args.out.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record))


if __name__ == '__main__':
    main()
