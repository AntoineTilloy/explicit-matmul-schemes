#!/usr/bin/env bash
set -u
cd /home/antoine/matmul_oct14/sub265/filter_theory/lazy_full_checks_20261014_01
printf '%s\n' "$$" > controller.pid
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sha256sum source.py graded_lazy_scheme.py catalogue.json ../positive_grading_20261014_01/result.json ../full_balanced_20261014_02/result.json > launch_hashes.txt
(ulimit -t 20; ulimit -v 524288; exec taskset -c 7 /usr/bin/time -v -o resources.txt /home/antoine/miniforge3/bin/python source.py --grading ../positive_grading_20261014_01/result.json --grading-sha c2199f0c9b59129ac18eeb71b8da891d433a1a5852a9b3db2f02a7f8dd72fd0c --candidate ../full_balanced_20261014_02/result.json --candidate-sha 0548ffa84ee20c7270f5635abfac7fa3e5f46b1368ac065cec11ffc60f7578be --generator-sha 7c4a122aa94b75c70ca736faa494b45bbf4bf51c2d56cc390632c06f56db8823 --catalogue catalogue.json --catalogue-sha 36922b9293dba32191be65b09f2ba5dd23038f11c778aa33ca9c2c8d0a4df628 --out .) > stdout.log 2> stderr.log
code=$?
printf '{"status":"TERMINAL","exit_code":%s}\n' "$code" > terminal.json
exit "$code"
