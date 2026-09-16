#!/usr/bin/env python3
"""Finite, model-free staged run. STOP is handled by the solver after an accepted step."""
import argparse
import copy
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT.parent


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def atomic(path, obj):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def source_hashes():
    files = [PROJECT/'CMakeLists.txt', *sorted((PROJECT/'include').rglob('*.hpp')),
             *sorted((PROJECT/'src').glob('*.cpp')),
             *sorted((PROJECT/'scripts').glob('*.py')), *sorted((PROJECT/'tests').glob('*.*'))]
    return {str(p): digest(p) for p in files if p.is_file()}


def source_tree_hash():
    files = sorted([PROJECT/'CMakeLists.txt', *sorted((PROJECT/'include').rglob('*.hpp')),
                    *sorted((PROJECT/'src').glob('*.cpp'))])
    return hashlib.sha256(''.join(f'{p.relative_to(PROJECT)}\t{digest(p)}\n' for p in files).encode()).hexdigest()


class Tail:
    def __init__(self, path):
        self.path = Path(path); self.offset = 0; self.header = None; self.last = {}

    def read(self):
        if not self.path.exists(): return self.last
        with self.path.open() as f:
            f.seek(self.offset)
            while True:
                start = f.tell(); line = f.readline()
                if not line or not line.endswith('\n'):
                    self.offset = start; break
                self.offset = f.tell()
                if self.header is None: self.header = line.rstrip().split('\t')
                else:
                    values=line.rstrip().split('\t')
                    if len(values)==len(self.header): self.last=dict(zip(self.header,values))
        return self.last


def run_process(command, output, log, deadline, progress_path, stage,
                initial_volume=None, run_started=None, elapsed_before=0.0):
    output.mkdir(parents=True, exist_ok=True)
    tail = Tail(output/'diagnostics.tsv')
    with log.open('a') as stream:
        process = subprocess.Popen(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
        stop_reason = None
        while process.poll() is None:
            row = tail.read()
            if time.monotonic() >= deadline: stop_reason = 'wall_clock_limit'
            if sum(p.stat().st_size for p in output.glob('*') if p.is_file()) > 12*1024**3:
                stop_reason = 'output_disk_limit'
            if stop_reason: (output/'STOP').touch(exist_ok=True)
            atomic(progress_path, {'stage':stage, 'status':stop_reason or 'running',
                'pid':process.pid, 'last_accepted_step':int(row.get('step',0)),
                'physical_time_s':float(row.get('time_s',0)), 'dt_s':float(row.get('dt_s',0)),
                'gas_volume_reduction_fraction':(1-float(row['gas_volume_main_flow_component'])/initial_volume)
                    if initial_volume and row else None,
                'free_moles_main':float(row.get('free_mol_main_flow_component',0)),
                'epsilon_N':float(row.get('epsilon_N_adjusted',0)),
                'epsilon_V':float(row.get('epsilon_V_adjusted',0)),
                'elapsed_seconds':elapsed_before + max(0, time.monotonic()-(run_started or time.monotonic())),
                'last_checkpoint':str(output/'checkpoint.json') if (output/'checkpoint.json').exists() else None,
                'failure_reason':stop_reason})
            # No log growth is NOT a failure. No signal/kill is sent to the solver.
            try: process.wait(timeout=60)
            except subprocess.TimeoutExpired: pass
        return process.returncode, stop_reason


def audit(output, expected_end=None):
    import math
    output = Path(output)
    m=json.loads((output/'run_manifest.json').read_text())
    if m['exit_status']!=0: raise RuntimeError('solver manifest reports failure')
    if m['build']['project_source_tree_sha256']!=source_tree_hash():
        raise RuntimeError('manifest/source tree mismatch')
    for mapping in ('inputs','outputs'):
        for p,h in m[mapping].items():
            if digest(p)!=h: raise RuntimeError(f'{mapping} hash mismatch: {p}')
    for p in m.get('input_provenance',{}).values():
        if digest(p['path'])!=p['sha256']: raise RuntimeError('provenance hash mismatch')
    c=m['configuration']; bg=c['pressure_reference']['background_abs_Pa']
    if c['boundary']['species_mode']!='copy_adjacent_concentration': raise RuntimeError('species mode')
    if c['nonlinear']['residual_tolerance']!=1e-8 or c['nonlinear']['update_tolerance']!=1e-9:
        raise RuntimeError('nonlinear tolerances changed')
    rows=list(csv.DictReader((output/'diagnostics.tsv').open(),delimiter='\t'))
    if not rows: raise RuntimeError('no accepted steps')
    for r in rows:
        if None in r or not all(math.isfinite(float(v)) for v in r.values()):
            raise RuntimeError('malformed/nonfinite diagnostics')
        for key,limit in [('residual_scaled_inf',1e-8),('relative_update_norm',1e-9),
                          ('epsilon_N_equation',1e-6),('epsilon_N_adjusted',1e-6),
                          ('epsilon_V_equation',1e-6),('epsilon_V_adjusted',1e-6)]:
            if not 0<=float(r[key])<=limit: raise RuntimeError(f'gate {key}')
        if float(r['boundary_diffusive_component_outflow'])!=0: raise RuntimeError('boundary diffusion')
    last=rows[-1]; cp=json.loads((output/'checkpoint.json').read_text())
    if cp['format'] not in ('bubble_checkpoint_v5','bubble_checkpoint_v6') or cp['background_abs_Pa']!=bg: raise RuntimeError('checkpoint representation')
    if cp['time_s']!=float(last['time_s']): raise RuntimeError('checkpoint/diagnostic time')
    if expected_end and cp['time_s']<expected_end and m['termination_reason']=='configured_end_time':
        raise RuntimeError('stage incomplete')
    params=json.loads((output/'physical_parameter_audit.json').read_text())
    components=json.loads((output/'component_boundary_and_makeup_audit.json').read_text())
    if not params['all_active_bubbles_initially_dissolving']: raise RuntimeError('initial driving force')
    # Cumulative balance is reconstructed independently from diagnostics, including restart segments.
    net=sum(float(r['dt_s'])*float(r['boundary_liquid_outflow_m3_s']) for r in rows)
    component_net=sum(float(r['dt_s'])*float(r['boundary_component_outflow_mol_s']) for r in rows)
    vl=sum(float(r['active_set_liquid_volume_ledger']) for r in rows)
    ml=sum(float(r['active_set_mole_ledger']) for r in rows)
    dv=float(last['cumulative_liquid_volume_change'])
    volume_raw=dv+net; volume_adjusted=volume_raw-vl
    makeup_scale=max(abs(dv),abs(net),1e-18)
    n_initial=params['initial_free_gas_moles']+params['initial_dissolved_moles']
    dn=float(last['free_mol_total'])+float(last['dissolved_mol_total'])-n_initial
    if abs(volume_adjusted)/makeup_scale>1e-6: raise RuntimeError('actual-change makeup gate')
    if abs(dn+component_net-ml)/max(n_initial,1e-21)>1e-6: raise RuntimeError('cumulative moles gate')
    result={'passed':True,'physical_time_s':cp['time_s'],'steps':len(rows),
        'volume_raw_error_m3':volume_raw,'volume_adjusted_error_m3':volume_adjusted,
        'makeup_relative_error':abs(volume_adjusted)/makeup_scale,'makeup_floor_m3':1e-18,
        'volume_ledger_to_actual_change':vl/makeup_scale,
        'mole_raw_error_mol':dn+component_net,'mole_adjusted_error_mol':dn+component_net-ml,
        'mole_ledger_mol':ml,
        'gas_volume_reduction_fraction':1-float(last['gas_volume_main_flow_component'])/components['main_flow_component']['initial_gas_volume_m3'],
        'free_moles_reduction_fraction':1-float(last['free_mol_main_flow_component'])/components['main_flow_component']['initial_free_gas_mol'],
        'last':last,'initial':params,'component_classes':components['classes'],
        'maxima':{key:max(float(r[key]) for r in rows) for key in
            ['newton_iterations','timestep_retries','residual_scaled_inf','relative_update_norm',
             'update_criterion_ratio','epsilon_N_adjusted','epsilon_V_adjusted']},
        'termination_reason':m['termination_reason']}
    atomic(output/'conservation_and_pressure_audit.json',result)
    return result


def accuracy(coarse, fine):
    a=json.loads((Path(coarse)/'checkpoint.json').read_text());b=json.loads((Path(fine)/'checkpoint.json').read_text())
    ra=list(csv.DictReader((Path(coarse)/'diagnostics.tsv').open(),delimiter='\t'))[-1]
    rb=list(csv.DictReader((Path(fine)/'diagnostics.tsv').open(),delimiter='\t'))[-1]
    errors={k:abs(float(ra[k])-float(rb[k]))/max(abs(float(rb[k])),1e-21) for k in
        ['free_mol_main_flow_component','gas_volume_main_flow_component','dissolved_mol_total']}
    selected=[i for i,s in enumerate(a['state']) if s[4]][:16]
    state_error=0.0
    for i in selected:
        for j,scale in [(0,100),(1,10000),(2,1),(3,1)]:
            x=float(a['state_extended_precision'][i][j]);y=float(b['state_extended_precision'][i][j])
            state_error=max(state_error,abs(x-y)/max(abs(y),scale))
    return {'relative_inventory_errors':errors,'representative_state_scaled_error':state_error,
            'representative_node_ids':[i+1 for i in selected],
            'tolerance':0.01,'passed':max([state_error,*errors.values()])<=0.01}


def postprocess(output):
    subprocess.run([sys.executable,str(PROJECT/'scripts/relative_visualize.py'),str(output)],check=True,cwd=ROOT)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--dry-run',action='store_true')
    ap.add_argument('--resume',action='store_true');ap.add_argument('--max-wall-seconds',type=float,default=7200)
    args=ap.parse_args();c=json.loads(args.config.read_text());binary=args.binary.resolve()
    out=(ROOT/c['io']['output_directory']).resolve()
    protected=['full_network_smoke','full_network_physical_um_virtual_boundary','full_network_preview']
    if any(out.name.startswith(n) for n in protected): raise RuntimeError('protected output')
    if not binary.is_file(): raise RuntimeError('missing binary')
    for key in ['pore_file','throat_file','connect_audit_file','input_geometry_audit_file','virtual_boundary_geometry_audit_file']:
        if not (ROOT/c['io'][key]).is_file(): raise RuntimeError(f'missing {key}')
    stages=[1,20,50,100,200]
    plan={'output':str(out),'binary':str(binary),'stages_s':stages,'time_accuracy_tolerance':0.01,
          'dt_policy':'A coarse+half; if >1%, stop before continuation and require an explicit finer-dt fresh run',
          'wall_seconds':args.max_wall_seconds,'stop':'10% main gas volume loss, transfer stop, 200 s, or wall budget',
          'container_path':str(out),'wsl_path':str(out).replace('/workspace/','/home/rong/Mycode/')}
    if args.dry_run: print(json.dumps(plan,indent=2));return
    out.mkdir(parents=True,exist_ok=True)
    lock=(out/'run.lock').open('a+')
    try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: raise SystemExit('run already locked; duplicate launch refused')
    if (out/'run_summary.json').exists() and not args.resume: raise SystemExit('existing run: use --resume')
    if args.resume:
        (out/'STOP').unlink(missing_ok=True)
    original_hashes=source_hashes(); binary_hash=digest(binary)
    metadata={'pid':os.getpid(),'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
              'container':os.environ.get('BUBBLE_CONTAINER_ID',Path('/etc/hostname').read_text().strip()),
              'source_hashes':original_hashes,'binary_sha256':binary_hash,'config_sha256':digest(args.config),**plan}
    prior=json.loads((out/'run_summary.json').read_text()) if args.resume and (out/'run_summary.json').exists() else {}
    if prior and prior.get('source_hashes')!=original_hashes: raise RuntimeError('resume source/script hash mismatch')
    spent=float(prior.get('elapsed_seconds',0));start=time.monotonic();deadline=start+args.max_wall_seconds-spent
    summary={**metadata,'status':'running','model_final_acceptance_completed':False,'stage':'preflight',
             'exit_code':None,'last_checkpoint':None,'elapsed_seconds':spent,'audits':prior.get('audits',[])}
    atomic(out/'run_metadata.json',metadata)
    status='failed';phase='preflight';exit_code=1
    try:
        baseline=json.loads((PROJECT/'relative_pressure_baseline.json').read_text())
        for p,h in baseline.items():
            if Path(p).name.startswith('PNM_') or Path(p).name=='trans-react.cpp':
                if digest(p)!=h: raise RuntimeError(f'protected input changed: {p}')
        cp=out/'checkpoint.json'; current=json.loads(cp.read_text())['time_s'] if args.resume and cp.exists() else 0.0
        for end in stages:
            if current>=end: continue
            if time.monotonic()>=deadline: status='wall_clock_limit';break
            if source_hashes()!=original_hashes or digest(binary)!=binary_hash: raise RuntimeError('code changed during run')
            phase=f'to_{end}s';resolved=copy.deepcopy(c);resolved['time']['t_end_s']=end
            if end==1: resolved['time']['stop_main_gas_volume_fraction']=0.0
            config_path=out/f'config_to_{end}s.json';atomic(config_path,resolved)
            command=[str(binary),'--config',str(config_path)]
            if current>0: command+=['--restart',str(cp)]
            exit_code,limit=run_process(command,out,out/f'solver_to_{end}s.log',deadline,
                out/'progress.json',phase,run_started=start,elapsed_before=spent)
            if exit_code: raise RuntimeError(f'solver exit {exit_code}; see solver_to_{end}s.log')
            check=audit(out,end);summary['audits'].append(check)
            current=check['physical_time_s'];summary['last_checkpoint']=str(cp)
            if limit or check['termination_reason']=='safe_stop_requested': status=limit or 'safely_stopped';break
            if end==1:
                fine=out/'accuracy_half_dt';fc=copy.deepcopy(resolved)
                fc['io']['output_directory']=str(fine);fc['time']['initial_dummy']=0
                del fc['time']['initial_dummy']
                fc['time']['dt_initial_s']*=0.5;fc['time']['dt_max_s']*=0.5
                fp=out/'config_accuracy_half_dt.json';atomic(fp,fc)
                rc,limit=run_process([str(binary),'--config',str(fp)],fine,
                    out/'accuracy_solver.log',deadline,out/'progress.json','accuracy_half_dt',
                    run_started=start,elapsed_before=spent)
                if rc: raise RuntimeError('time accuracy comparison solver failed')
                audit(fine,1);comparison=accuracy(out,fine);atomic(out/'time_accuracy_audit.json',comparison)
                if limit: status=limit;break
                if not comparison['passed']:
                    # Never silently change the continuation signature. A new refined run is required.
                    raise RuntimeError('time accuracy >1%; finer dt_max required; stopped before long run')
            postprocess(out)
            if check['gas_volume_reduction_fraction']>=0.1:
                status='completed';summary['stop_reason']='main_gas_volume_reduced_10_percent';break
            if check['termination_reason']=='no_positive_transfer_for_five_accepted_steps':
                status='completed';summary['stop_reason']='sustained_no_positive_transfer_not_thermodynamic_equilibrium';break
        else: status='completed';summary['stop_reason']='maximum_physical_time_200s'
        exit_code=0 if status=='completed' else 2
    except Exception as error:
        summary['failure_reason']=str(error)
    finally:
        summary.update(status=status,stage=phase,exit_code=exit_code,
            elapsed_seconds=spent+time.monotonic()-start,
            last_checkpoint=str(out/'checkpoint.json') if (out/'checkpoint.json').exists() else None,
            resume_command=f'{sys.executable} {Path(__file__).resolve()} --config {args.config.resolve()} --binary {binary} --resume')
        atomic(out/'run_summary.json',summary)
        (out/'RUN_SUMMARY.md').write_text('# Local staged run\n\n'+json.dumps({k:v for k,v in summary.items() if k not in ('source_hashes','audits')},indent=2)+'\n')
        atomic(out/'progress.json',{'stage':phase,'status':status,'failure_reason':summary.get('failure_reason'),
               'elapsed_seconds':summary['elapsed_seconds'],'last_checkpoint':summary['last_checkpoint']})
    raise SystemExit(exit_code)


if __name__=='__main__': main()
