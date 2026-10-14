#!/bin/bash
set -u
cd /home/antoine/matmul_oct14/sub265/filter_theory/independent_balanced_20261014_02 || exit 97
printf '%s\n' "$$" > controller.pid
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sha256sum source.py grading.json batched.json catalogue.json run.sh scope.txt > bindings.sha256
(ulimit -t 10; ulimit -v 524288; exec taskset -c 7 /usr/bin/time -v -o resources.txt /home/antoine/miniforge3/bin/python source.py --grading grading.json --result batched.json --expected-result-sha 8ccb38e927999662d242cd44b3e4f188fff4ab155806091a4564d283fa2cfa84 --catalogue catalogue.json --out result.json --max-pairs 8000000) >stdout.log 2>stderr.log
rc=$?
printf '{"exit_code":%s,"terminal":true}\n' "$rc" > terminal.json
if test -f result.json; then sha256sum result.json > result.sha256; fi
exit "$rc"
