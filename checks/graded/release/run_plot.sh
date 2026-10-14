#!/bin/bash
set -u
release_root=/home/antoine/matmul_oct14/sub265/release_20261014_01
cd "$release_root/repo"
printf '%s\n' "$$" > "$release_root/plot.pid"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sha256sum comparisons/build.py comparisons/graded_plot.py catalogue.json > "$release_root/plot.bindings"
(ulimit -t 120; ulimit -v 1048576; exec taskset -c 8 /usr/bin/time -v -o "$release_root/plot.resources" /home/antoine/matmul_oct14/sub265/venv/bin/python comparisons/build.py) > "$release_root/plot.log" 2> "$release_root/plot.stderr"
code=$?
printf '{"exit_code":%s,"terminal":true}\n' "$code" > "$release_root/plot.terminal.json"
exit "$code"
