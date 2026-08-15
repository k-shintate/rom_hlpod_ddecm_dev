#!/usr/bin/env python3
from pathlib import Path
import argparse, json, re, sys

default_thresholds=Path(__file__).resolve().with_name('validated_thresholds.json')
p=argparse.ArgumentParser(description='Verify the validated dG(1), two-slab, split-basis ST-DDROM case.')
p.add_argument('result_dir', type=Path)
p.add_argument('--thresholds', type=Path, default=default_thresholds)
a=p.parse_args()
T=json.loads(a.thresholds.read_text())
E=T['expected']; tol=T['tolerance']
log=a.result_dir/f"fluid_st_ddrom_online_np{E['mpi_ranks']}.log"
if not log.is_file(): raise SystemExit(f'ERROR: missing online log: {log}')
text=log.read_text(errors='replace')
checks=[]
def need(pattern, label):
    ok=re.search(pattern,text,re.M) is not None; checks.append((label,ok)); return ok
need(rf'dG degree\s*:\s*{E["degree"]}\b','dG degree')
need(rf'slabs/window\s*:\s*{E["slabs_per_window"]}\b','slabs/window')
need(rf'online solve windows\s*:\s*{E["online_windows"]}\b','online windows')
need(rf'physical MPI ranks\s*:\s*{E["mpi_ranks"]}\b','MPI ranks')
need(rf'velocity modes/rank={E["velocity_modes"]}\s+pressure modes/rank={E["pressure_modes"]}','split basis modes')

acc_re=re.compile(r'fluid ST-DDROM accuracy window\s+(\d+):\s+projection=([0-9Ee+.-]+)\s+state=([0-9Ee+.-]+)\s+dynamic=([0-9Ee+.-]+)\s+velocity\(state\)=([0-9Ee+.-]+)\s+pressure\(state\)=([0-9Ee+.-]+)\s+coef=([0-9Ee+.-]+)\s+trace=([0-9Ee+.-]+)\s+FOM residual=([0-9Ee+.-]+)\s+FOM projected defect=([0-9Ee+.-]+)')
rows=[]
for m in acc_re.finditer(text):
    rows.append((int(m.group(1)),)+tuple(float(m.group(i)) for i in range(2,11)))
if len(rows)!=E['online_windows']:
    checks.append((f'accuracy rows {len(rows)}/{E["online_windows"]}',False))
else: checks.append(('accuracy row count',True))
keys=['projection','state','dynamic','velocity_state','pressure_state','coefficient','trace','fom_residual','fom_projected_defect']
maxima={k:0.0 for k in keys}
for row in rows:
    for k,v in zip(keys,row[1:]): maxima[k]=max(maxima[k],abs(v))
for k in keys: checks.append((f'{k} max={maxima[k]:.3e} <= {tol[k]:.3e}', maxima[k] <= tol[k]))

complete_re=re.compile(r'fluid ST-DDROM window\s+(\d+) completed:.*?projected rel residual=([0-9Ee+.-]+)\s+full rel residual=([0-9Ee+.-]+)')
completed=[(int(m.group(1)),float(m.group(2)),float(m.group(3))) for m in complete_re.finditer(text)]
checks.append((f'completed windows {len(completed)}/{E["online_windows"]}',len(completed)==E['online_windows']))
if completed:
    mp=max(abs(x[1]) for x in completed); mf=max(abs(x[2]) for x in completed)
    checks.append((f'ROM projected relative max={mp:.3e}',mp<=tol['rom_projected_relative']))
    checks.append((f'ROM full relative max={mf:.3e}',mf<=tol['rom_full_relative']))

snap=a.result_dir/f"fluid_st_windows_s{E['slabs_per_window']}"
for rank in range(E['mpi_ranks']):
    count=0
    while (snap/f"fluid_st_window_{count:04d}_rank_{rank:06d}.bin").is_file(): count+=1
    checks.append((f'snapshots rank {rank}: {count} >= {E["online_windows"]}',count>=E['online_windows']))

basis=a.result_dir/'fluid_st_ddrom_basis'
for rank in range(E['mpi_ranks']):
    checks.append((f'basis rank {rank}',(basis/f"stddrom_basis_rank_{rank:06d}.bin").is_file()))

vtk_groups={}
for f in a.result_dir.glob('st_ddrom_result_*.vtk.*'):
    m=re.match(r'(st_ddrom_result_\d+\.vtk)\.(\d+)$',f.name)
    if m: vtk_groups.setdefault(m.group(1),set()).add(int(m.group(2)))
for stem,ranks in sorted(vtk_groups.items()):
    checks.append((f'VTK rank set {stem}',ranks==set(range(E['mpi_ranks']))))

for label,ok in checks: print(('PASS' if ok else 'FAIL')+': '+label)
failed=[x for x in checks if not x[1]]
print(f'\nvalidated windows parsed: {len(rows)}')
for k in keys: print(f'max {k:24s}: {maxima[k]:.15e}')
if failed: raise SystemExit(1)
print('\nVALIDATION STATUS: PASS')
