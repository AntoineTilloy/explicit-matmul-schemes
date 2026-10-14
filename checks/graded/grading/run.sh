#!/bin/bash
set -u
cd /home/antoine/matmul_oct14/sub265/filter_theory/positive_grading_20261014_01
printf '%s\n' "$$" > controller.pid
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sha256sum source.py al_puzzles.py run.sh > bindings.sha256
(ulimit -t 10; ulimit -v 524288; exec taskset -c 7 /usr/bin/time -v -o resources.txt /home/antoine/miniforge3/bin/python source.py --puzzles al_puzzles.py --out .) > stdout.log 2> stderr.log
rc=$?
printf '{"exit_code":%s,"terminal":true}\n' "$rc" > terminal.json
exit "$rc"
