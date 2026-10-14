#!/bin/bash
set -u
cd /home/antoine/matmul_oct14/sub265/filter_theory/catalogue_batch_20261014_01 || exit 97
printf '%s\n' "$$" > controller.pid
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sha256sum source.py full.py helper.py cached.json grading.json catalogue.json run.sh > bindings.sha256
/home/antoine/matmul_oct14/sub265/venv/bin/python -c 'import flint,sys; print(sys.executable); print(sys.version); print(flint.__version__)' > runtime.txt
(ulimit -t 20; ulimit -v 524288; exec taskset -c 7 /usr/bin/time -v -o resources.txt /home/antoine/matmul_oct14/sub265/venv/bin/python source.py --full-source full.py --helper helper.py --cached cached.json --grading grading.json --catalogue catalogue.json --out .) >stdout.log 2>stderr.log
rc=$?
printf '{"exit_code":%s,"terminal":true}\n' "$rc" > terminal.json
exit "$rc"
