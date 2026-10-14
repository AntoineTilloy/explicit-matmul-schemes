#!/bin/bash
set -u
release_root=/home/antoine/matmul_oct14/sub265/release_20261014_01
cd "$release_root/repo"
printf '%s\n' "$$" > "$release_root/rectangular_tail.pid"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CWSCHEMES_CACHE="$release_root/cache" PYTHONPATH="$release_root/repo/src"
(ulimit -t 60; ulimit -v 4194304; exec taskset -c 7 /usr/bin/time -v -o "$release_root/rectangular_tail.resources" /home/antoine/matmul_oct14/sub265/venv/bin/python -m pytest -q tests/test_release.py::test_new_rectangular_catalogue_entries) > "$release_root/rectangular_tail.log" 2> "$release_root/rectangular_tail.stderr"
code=$?
printf '{"exit_code":%s,"terminal":true}\n' "$code" > "$release_root/rectangular_tail.terminal.json"
exit "$code"
