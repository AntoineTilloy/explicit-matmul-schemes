"""Bind completed Atlas release checks to the exported finite constructions."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

root=Path(__file__).resolve().parent
repo=root/'repo'
out=repo/'checks/graded/release'
def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

assert read(root/'export.terminal.json')['exit_code']==0
assert read(root/'plot.terminal.json')['exit_code']==0
assert read(root/'wheel.terminal.json')['exit_code']==0
assert read(root/'remaining.terminal.json')['exit_code']==1
assert read(root/'rectangular_tail.terminal.json')['exit_code']==0
assert read(root/'tests.terminal.json')['exit_code']==137
assert (root/'tests.log').read_text().strip()=='.'*40
assert re.search(r'1 failed, 15 passed in ',(root/'remaining.log').read_text())
assert 'MemoryError' in (root/'remaining.log').read_text()
assert re.search(r'1 passed in ',(root/'rectangular_tail.log').read_text())
collected=subprocess.check_output([sys.executable,'-m','pytest','--collect-only','-q'],cwd=repo,text=True)
names=[line for line in collected.splitlines() if '::test_' in line]
assert len(names)==56
assert names[40]=='tests/test_release.py::test_z4_quotient_catalogue_counts'
wheel=read(root/'wheel.smoke.json')
assert wheel['status']=='PASS' and wheel['catalogue_entries']==99
catalogue=read(repo/'catalogue.json')
assert (repo/'catalogue.json').read_bytes()==(repo/'src/cwschemes/data/catalogue.json').read_bytes()
previous=read(repo/'checks/graded/candidate/catalogue.json')
by_n={x['n']:x for x in catalogue}
assert len(catalogue)==99 and len(previous)==97
assert all(by_n[x['n']]==x for x in previous)
graded=[x for x in catalogue if 'graded' in x['provenance']['construction'].lower()]
assert len(graded)==2
rank=4*31681*7**210
assert all(int(x['rank'])==rank and rank**20<int(x['n'])**53 for x in graded)
assert any(int(x['n'])==2**48*6**70 for x in graded)
assert read(repo/'figures/exponents.metadata.json')['catalogue_sha256']==sha(repo/'catalogue.json')
assert read(repo/'checks/graded/audit/result.json')['status']=='PASS'
assert read(repo/'checks/graded/controls/result.json')['status']=='COMPLETE'
out.mkdir(parents=True,exist_ok=True)
receipt_names=['run_export.sh','run_checks.sh','run_remaining.sh','run_rectangular_tail.sh','run_plot.sh','run_wheel.sh',
               'wheel_smoke.py','summarize_release.py']
for prefix in ('export','tests','remaining','rectangular_tail','plot','wheel'):
    receipt_names.extend(p.name for p in root.glob(prefix+'.*') if p.is_file())
for name in sorted(set(receipt_names)):
    shutil.copy2(root/name,out/name)
(out/'collected_tests.txt').write_text(collected)
summary=dict(status='PASS',catalogue_entries=99,prior_catalogue_entries_unchanged=97,
    tests=dict(passed=56,first_run_passed=40,continuation_passed=15,isolated_tail_passed=1,
               resource_retries='First run reached the 180 CPU-second limit during the existing Z4 catalogue test. The continuation passed 15 checks, then the existing rectangular catalogue test hit the 3 GiB combined-process memory cap. That last test passed in a fresh isolated process.',
               scope='All collected regression tests completed across the three runs; only resource-limited tests were repeated.'),
    installed_wheel='PASS from outside the repository with a fresh resource cache',
    exact_threshold='R^20 < n^53 for both new rational entries',
    rank_formula='4 * 31681 * 7^210',natural_size_formula='2^48 * 6^70',
    search_scope='Bounded catalogue-base powers, LITA32/Strassen fillers, up to64 packets; no global minimum claim',
    compute='Atlas only; no new Vast spend',
    limits='The large square tensor and its immense composed coefficient sums are not expanded; exact finite proof and implementation checks have the scope documented in docs/graded_cw.md.',
    catalogue_sha256=sha(repo/'catalogue.json'),
    figures={ext:sha(repo/'figures'/('exponents.'+ext)) for ext in ('svg','png','pdf')},
    receipts={p.name:sha(p) for p in out.iterdir() if p.is_file()})
(repo/'checks/graded/release.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(dict(status='PASS',tests=56,catalogue_entries=99,release_summary_sha256=sha(repo/'checks/graded/release.json'))))
