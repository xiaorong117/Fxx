#!/usr/bin/env python3
import csv,json,subprocess,sys,tempfile
from pathlib import Path
binary=Path(sys.argv[1]).resolve();base=json.loads((Path(__file__).parent/'data/adaptive_small.json').read_text())
with tempfile.TemporaryDirectory(prefix='adaptive_dt_') as d:
 root=Path(d); fixed=json.loads(json.dumps(base));fixed['io']['output_directory']=str(root/'fixed');fixed.pop('time_step_control',None)
 adaptive=json.loads(json.dumps(base));adaptive['io']['output_directory']=str(root/'adaptive')
 for name,c in [('fixed',fixed),('adaptive',adaptive)]:
  p=root/(name+'.json');p.write_text(json.dumps(c));r=subprocess.run([str(binary),'--config',str(p)],text=True,capture_output=True)
  if r.returncode:raise SystemExit(f'{name} run failed:\n{r.stdout}\n{r.stderr}')
 rows=[]
 with (root/'adaptive'/'diagnostics.tsv').open(newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
 if not rows or not all(float(r['epsilon_N_adjusted'])<=1e-6 and float(r['epsilon_V_adjusted'])<=1e-6 for r in rows):raise SystemExit('adaptive conservation failed')
 if not all(r['dt_change_reason'] for r in rows):raise SystemExit('adaptive dt reasons missing')
 cp=json.loads((root/'adaptive'/'checkpoint.json').read_text())
 if cp['next_dt_s']<1e-6 or cp['next_dt_s']>0.1:raise SystemExit('adaptive checkpoint next dt out of bounds')
 fixed_rows=list(csv.DictReader((root/'fixed'/'diagnostics.tsv').open(),delimiter='\t'))
 print(f"PASS adaptive dt fixed_steps={len(fixed_rows)} adaptive_steps={len(rows)} reasons={sorted(set(r['dt_change_reason'] for r in rows))}")
