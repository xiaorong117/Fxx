#!/usr/bin/env python3
"""Finite local AMGX production driver; no model/API dependency."""
import argparse,csv,fcntl,hashlib,json,os
from pathlib import Path
import subprocess,time

PROJECT=Path(__file__).resolve().parents[1];ROOT=PROJECT.parent
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic(path,data):
    p=Path(path);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(data,indent=2)+'\n');q.replace(p)
def latest(path):
    if not path.exists():return {}
    with path.open() as stream:
        rows=list(csv.DictReader(stream,delimiter='\t'))
    return rows[-1] if rows else {}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True)
    ap.add_argument('--resume',action='store_true');ap.add_argument('--dry-run',action='store_true');args=ap.parse_args()
    config=json.loads(args.config.read_text());out=(ROOT/config['io']['output_directory']).resolve();binary=args.binary.resolve()
    if config['linear']['backend']!='amgx':raise SystemExit('production config must select AMGX')
    if config['termination']['maximum_wall_seconds']!=0:raise SystemExit('production wall limit must be disabled')
    if out.name!='full_network_air_water_amgx_sg_threshold_run01':raise SystemExit('unexpected production output')
    plan={'binary':str(binary),'binary_sha256':digest(binary),'config':str(args.config.resolve()),
          'config_sha256':digest(args.config),'output':str(out),'backend':'amgx',
          'target_gas_saturation':config['termination']['target_gas_saturation'],'wall_limit_seconds':0}
    if args.dry_run:print(json.dumps(plan,indent=2));return
    out.mkdir(parents=True,exist_ok=True);lock=(out/'run.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit('duplicate production launch refused by run.lock')
    if (out/'run_summary.json').exists() and not args.resume:raise SystemExit('existing production run requires --resume')
    command=[str(binary),'--config',str(args.config.resolve())]
    if args.resume:
        (out/'STOP').unlink(missing_ok=True);command+=['--restart',str(out/'checkpoint.json')]
    metadata={**plan,'pid':os.getpid(),'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
              'container':Path('/etc/hostname').read_text().strip(),'command':command}
    atomic(out/'run_metadata.json',metadata);started=time.monotonic()
    with (out/'driver.log').open('a') as log:
        process=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        while process.poll() is None:
            row=latest(out/'global_history.tsv')
            atomic(out/'progress.json',{'status':'running','pid':process.pid,'step':int(row.get('step',0)),
              'time_s':float(row.get('time_s',0)),'dt_s':float(row.get('dt_s',0)),
              'main_real_gas_saturation':float(row.get('main_real_gas_saturation',0)),
              'gas_volume_retained_fraction':float(row.get('gas_volume_retained_fraction',1)),
              'free_gas_moles_retained_fraction':float(row.get('free_gas_moles_retained_fraction',1)),
              'component_balance_relative_error':float(latest(out/'balance_history.tsv').get('component_balance_relative_error',0)),
              'liquid_balance_relative_error':float(latest(out/'balance_history.tsv').get('liquid_balance_relative_error',0)),
              'elapsed_seconds':time.monotonic()-started,
              'checkpoint':str(out/'checkpoint.json') if (out/'checkpoint.json').exists() else None})
            try:process.wait(timeout=60)
            except subprocess.TimeoutExpired:pass
    status='completed' if process.returncode==0 else 'failed';manifest={}
    if (out/'run_manifest.json').exists():manifest=json.loads((out/'run_manifest.json').read_text());status=manifest.get('termination_reason',status)
    summary={**metadata,'status':status,'exit_code':process.returncode,'elapsed_seconds':time.monotonic()-started,
             'final':latest(out/'global_history.tsv'),'manifest_termination_reason':manifest.get('termination_reason'),
             'model_final_acceptance_completed':False,
             'resume_command':f'{__file__} --config {args.config.resolve()} --binary {binary} --resume'}
    atomic(out/'run_summary.json',summary);(out/'RUN_SUMMARY.md').write_text('# AMGX production run\n\n```json\n'+json.dumps(summary,indent=2)+'\n```\n')
    atomic(out/'progress.json',{'status':status,'exit_code':process.returncode,'elapsed_seconds':summary['elapsed_seconds'],
                               'checkpoint':str(out/'checkpoint.json') if (out/'checkpoint.json').exists() else None})
    raise SystemExit(process.returncode)
if __name__=='__main__':main()
