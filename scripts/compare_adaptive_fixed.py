import csv,json,sys
from pathlib import Path
fixed,adaptive,out=map(Path,sys.argv[1:])
def rows(p):
 with (p/'diagnostics.tsv').open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
f,a=rows(fixed),rows(adaptive)
if not f or not a:raise SystemExit('missing fixed/adaptive diagnostics')
for name,r in [('fixed',f),('adaptive',a)]:
 for x in r:
  if float(x['epsilon_N_adjusted'])>1e-6 or float(x['epsilon_V_adjusted'])>1e-6:raise SystemExit(name+' conservation failure')
  if float(x['residual_scaled_inf'])>1e-8 or float(x['relative_update_norm'])>1e-9:raise SystemExit(name+' Newton gate failure')
x,y=f[-1],a[-1]
keys=['free_mol_total','dissolved_mol_total','gas_volume_total','internal_pressure_min_Pa','internal_pressure_max_Pa']
result={'fixed_steps':len(f),'adaptive_steps':len(a),'fixed_time_s':float(x['time_s']),'adaptive_time_s':float(y['time_s']),
 'dt_fixed':[float(r['dt_s']) for r in f],'dt_adaptive':[float(r['dt_s']) for r in a],
 'adaptive_reasons':sorted(set(r['dt_change_reason'] for r in a)),
 'adaptive_retries':sum(int(r['timestep_retries']) for r in a),
 'gas_disappearance_events':sum(int(r['gas_disappearance_event']) for r in a),
 'max_epsilon_N':max(float(r['epsilon_N_adjusted']) for r in a),'max_epsilon_V':max(float(r['epsilon_V_adjusted']) for r in a),
 'final_differences':{k:{'absolute':abs(float(x[k])-float(y[k])),'relative':abs(float(x[k])-float(y[k]))/max(abs(float(x[k])),abs(float(y[k])),1e-300)} for k in keys}}
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
