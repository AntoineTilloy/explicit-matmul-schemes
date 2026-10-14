#!/bin/bash
set -u
release_root=/home/antoine/matmul_oct14/sub265/release_20261014_01
cd "$release_root/repo"
printf '%s\n' "$$" > "$release_root/tests.pid"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CWSCHEMES_CACHE="$release_root/cache" PYTHONPATH="$release_root/repo/src"
sha256sum src/cwschemes/_engine/graded_packet_core.py src/cwschemes/_engine/graded_square.py tests/test_graded.py catalogue.json > "$release_root/tests.bindings"
(ulimit -t 180; ulimit -v 2097152; exec taskset -c 7 /usr/bin/time -v -o "$release_root/tests.resources" /home/antoine/matmul_oct14/sub265/venv/bin/python -m pytest -q) > "$release_root/tests.log" 2> "$release_root/tests.stderr"
code=$?
printf '{"exit_code":%s,"terminal":true}\n' "$code" > "$release_root/tests.terminal.json"
exit "$code"
