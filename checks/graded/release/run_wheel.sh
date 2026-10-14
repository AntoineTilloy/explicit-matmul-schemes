#!/bin/bash
set -u
release_root=/home/antoine/matmul_oct14/sub265/release_20261014_01
runtime_python=/home/antoine/matmul_oct14/sub265/venv/bin/python
cd "$release_root/repo"
printf '%s\n' "$$" > "$release_root/wheel.pid"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"$runtime_python" -m pip wheel . --no-deps --no-build-isolation --wheel-dir "$release_root/wheels" > "$release_root/wheel.build.log" 2>&1
code=$?
if test "$code" -eq 0; then
    "$runtime_python" -m pip install --no-deps --target "$release_root/wheel_install" "$release_root"/wheels/*.whl > "$release_root/wheel.install.log" 2>&1
    code=$?
fi
if test "$code" -eq 0; then
    cd /tmp
    export PYTHONPATH="$release_root/wheel_install" CWSCHEMES_CACHE="$release_root/wheel_cache"
    (ulimit -t 30; ulimit -v 1048576; exec taskset -c 8 /usr/bin/time -v -o "$release_root/wheel.resources" "$runtime_python" "$release_root/wheel_smoke.py") > "$release_root/wheel.smoke.json" 2> "$release_root/wheel.stderr"
    code=$?
fi
printf '{"exit_code":%s,"terminal":true}\n' "$code" > "$release_root/wheel.terminal.json"
exit "$code"
