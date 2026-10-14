#!/bin/bash
set -u
release_root=/home/antoine/matmul_oct14/sub265/release_20261014_01
cd "$release_root/repo"
printf '%s\n' "$$" > "$release_root/remaining.pid"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CWSCHEMES_CACHE="$release_root/cache" PYTHONPATH="$release_root/repo/src"
(ulimit -t 420; ulimit -v 3145728; exec taskset -c 7 /usr/bin/time -v -o "$release_root/remaining.resources" /home/antoine/matmul_oct14/sub265/venv/bin/python -m pytest -q \
  tests/test_release.py::test_z4_quotient_catalogue_counts \
  tests/test_release.py::test_z4_quotient_small_exact \
  tests/test_release.py::test_odd_character_survives_mixed_certificate \
  tests/test_release.py::test_uneven_strassen_full_tensor \
  tests/test_release.py::test_split_with_shared_child_full_tensor \
  tests/test_release.py::test_new_rectangular_catalogue_entries \
  tests/test_z4_certified_q.py) > "$release_root/remaining.log" 2> "$release_root/remaining.stderr"
code=$?
printf '{"exit_code":%s,"terminal":true}\n' "$code" > "$release_root/remaining.terminal.json"
exit "$code"
