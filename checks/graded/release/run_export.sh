#!/bin/bash
set -u
release_root=/home/antoine/matmul_oct14/sub265/release_20261014_01
input_root=/home/antoine/matmul_oct14/sub265/filter_theory
cd "$release_root"
printf '%s\n' "$$" > export.pid
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CWSCHEMES_CACHE="$release_root/cache"
sha256sum export_graded.py repo/src/cwschemes/_engine/graded_packet_core.py repo/src/cwschemes/_engine/graded_square.py > export.bindings
(ulimit -t 20; ulimit -v 524288; exec taskset -c 7 /usr/bin/time -v -o export.resources /home/antoine/matmul_oct14/sub265/venv/bin/python export_graded.py \
  --repo "$release_root/repo" \
  --candidate "$input_root/catalogue_batch_20261014_01/result.json" --candidate-sha 8ccb38e927999662d242cd44b3e4f188fff4ab155806091a4564d283fa2cfa84 \
  --grading "$input_root/positive_grading_20261014_01/result.json" --grading-sha c2199f0c9b59129ac18eeb71b8da891d433a1a5852a9b3db2f02a7f8dd72fd0c \
  --audit "$input_root/independent_balanced_20261014_02/result.json" --audit-sha 5c95fea9627527e3b0021d16cd2a5a610cdd4926dffa2cfdecaf195c1e87f5d4 \
  --controls "$input_root/lazy_full_checks_20261014_01/result.json" --controls-sha 27a76f59deb10930a2caed31c257bbfc59bf69b2dffc696c2892d64fc7172898 \
  --core-sha 7c4a122aa94b75c70ca736faa494b45bbf4bf51c2d56cc390632c06f56db8823 \
  --adapter-sha 6d0d4c9d1ed6d41e0d9d8dcc461754e9baec97cf98afcf0472a271195adea101 \
  --catalogue-sha 36922b9293dba32191be65b09f2ba5dd23038f11c778aa33ca9c2c8d0a4df628 \
  --index-sha e6a06cbeec23dfcfb112f0a7932ea6c2b21eae0289c7f9471358bf52883157b7 \
  --comparison-sha 125df6dfd47cc1a541d063e314e499366de2a3ff1790e07680a453a6e841ad89) > export.log 2> export.stderr
code=$?
printf '{"exit_code":%s,"terminal":true}\n' "$code" > export.terminal.json
exit "$code"
