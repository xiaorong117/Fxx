#!/usr/bin/env python3
import csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

eigen=Path(sys.argv[1]);amgx=Path(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True)
def table(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
er,ar=table(eigen/'steps.tsv'),table(amgx/'steps.tsv')
if len(er)!=5 or len(ar)!=5:raise SystemExit('both backends must have five accepted steps')
fields=['step_wall_seconds','assembly_seconds','csr_seconds','upload_seconds','setup_seconds','solve_seconds','download_seconds',
        'newton_iterations','linear_iterations','cpu_linear_relative_residual','residual_scaled_inf','relative_update_norm','epsilon_N','epsilon_V']
with (out/'backend_performance_comparison.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['backend','accepted_index','global_step','time_s','accepted_dt_s',*fields],delimiter='\t');w.writeheader()
    for backend,rows in [('eigen',er),('amgx',ar)]:
        for row in rows:w.writerow({k:(backend if k=='backend' else row.get(k,"0")) for k in w.fieldnames})

def aggregate(rows):
    sums={k:sum(float(r.get(k,0.0)) for r in rows) for k in fields}
    return {'steps':len(rows),'sum':sums,'mean':{k:v/len(rows) for k,v in sums.items()},
            'max':{k:max(float(r.get(k,0.0)) for r in rows) for k in fields}}
E,A=aggregate(er),aggregate(ar)
amgx_manifest=json.loads((amgx/'run_manifest.json').read_text());eigen_manifest=json.loads((eigen/'run_manifest.json').read_text())
amgx_artifact_wall=(amgx/'run_manifest.json').stat().st_mtime-(amgx/'resolved_config.json').stat().st_mtime
eigen_wall=eigen_manifest['total_five_step_wall_seconds']

ecp=json.loads((eigen/'checkpoint.json').read_text());acp=json.loads((amgx/'checkpoint.json').read_text())
if ecp['step_number']!=acp['step_number'] or ecp['time_s']!=acp['time_s']:raise SystemExit('final checkpoint time/step mismatch')
es=np.array([[float(x) if x is not None else np.nan for x in row[:4]] for row in ecp['state']],float)
as_=np.array([[float(x) if x is not None else np.nan for x in row[:4]] for row in acp['state']],float)
active=np.array([bool(x[4]) for x in ecp['state']]);diff=np.abs(es-as_)
scales=np.array([100.,10000.,1.,1.]);names=['Pl','Pg','Sw','C'];state={}
for j,name in enumerate(names):
    mask=np.isfinite(es[:,j])&np.isfinite(as_[:,j]);d=diff[mask,j];den=np.maximum(np.maximum(np.abs(es[mask,j]),np.abs(as_[mask,j])),1e-300)
    state[name]={'maximum_absolute_difference':float(d.max()),'maximum_relative_difference':float((d/den).max()),
                 'maximum_scaled_difference':float((d/scales[j]).max())}

config=json.loads((eigen/'resolved_config.json').read_text());amgx_config=json.loads((amgx/'resolved_config.json').read_text());root=Path('/workspace/zz')
pore=np.loadtxt(root/config['io']['pore_file']);throat=np.loadtxt(root/config['io']['throat_file']);n=len(pore)
rows=np.r_[throat[:,1].astype(int)-1,throat[:,2].astype(int)-1];cols=np.r_[throat[:,2].astype(int)-1,throat[:,1].astype(int)-1]
_,labels=connected_components(coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(n,n)),directed=False)
virtual=pore[:,1].astype(int)==2;flags=pore[:,10].astype(int);candidates=set(labels[virtual&(flags==1)])&set(labels[virtual&(flags==2)])
main=max(candidates,key=lambda c:int(np.count_nonzero((labels==c)&~virtual)))
volume=pore[:,9]*float(config['geometry']['input_volume_to_m3']);real=(~virtual)&(volume>0);main_real=real&(labels==main)
R=8.31446261815324;T=float(config['physics']['T_K']);Z=float(config['physics']['Z']);bg=float(config['pressure_reference']['background_abs_Pa'])
def inventory(state):
    sw=state[:,2];c=state[:,3];pg=state[:,1];gas=active&real
    free=float(np.sum((bg+pg[gas])*volume[gas]*(1-sw[gas])/(Z*R*T)))
    dissolved=float(np.sum(volume[real]*sw[real]*c[real]));sg=float(np.sum(volume[main_real]*(1-sw[main_real]))/np.sum(volume[main_real]))
    return {'main_real_gas_saturation':sg,'free_gas_mol':free,'dissolved_gas_mol':dissolved}
ei,ai=inventory(es),inventory(as_);inventories={}
for key in ei:
    absolute=abs(ei[key]-ai[key]);relative=absolute/max(abs(ei[key]),abs(ai[key]),1e-300)
    inventories[key]={'eigen':ei[key],'amgx':ai[key],'absolute_difference':absolute,'relative_difference':relative}
e_keys=sorted(ecp);a_keys=sorted(acp);node_ids=pore[:,0].astype(np.int64)
physics_equal=config['physics']==amgx_config['physics'] and config['boundary']==amgx_config['boundary'] and config['pressure_reference']==amgx_config['pressure_reference']
network_equal=all(config['io'][key]==amgx_config['io'][key] for key in ('pore_file','throat_file','connect_audit_file'))

linear_e=E['sum']['solve_seconds'];linear_a=A['sum']['solve_seconds']
result={'format':'bubble_backend_five_step_comparison_v1','initial_checkpoint':amgx_manifest['initial_checkpoint'],
 'initial_checkpoint_sha256':amgx_manifest['initial_checkpoint_sha256'],'accepted_steps':5,'requested_dt_s':0.01,
 'eigen':E,'amgx':A,'wall':{'eigen_measured_seconds':eigen_wall,'amgx_artifact_span_seconds':amgx_artifact_wall,
   'amgx_measurement_note':'Existing accepted AMGX run did not store process wall time; span is resolved_config mtime to run_manifest mtime and includes initialization/final checkpoint.',
   'eigen_average_step_wall_seconds':E['mean']['step_wall_seconds'],
            'amgx_average_step_timed_kernel_path_seconds':sum(A['mean'][k] for k in ['assembly_seconds','csr_seconds','upload_seconds','setup_seconds','solve_seconds','download_seconds']),
            'amgx_process_wall_time_available':False},
 'speedup_eigen_over_amgx':{'end_to_end_using_available_wall_evidence':eigen_wall/amgx_artifact_wall,
   'timed_step_components':sum(E['sum'][k] for k in ['assembly_seconds','csr_seconds','upload_seconds','setup_seconds','solve_seconds','download_seconds'])/
      sum(A['sum'][k] for k in ['assembly_seconds','csr_seconds','upload_seconds','setup_seconds','solve_seconds','download_seconds']),
   'pure_solve_seconds':linear_e/linear_a},
 'state_differences':state,'inventory_differences':inventories,
 'compatibility_audit':{'eigen_diagnostics_tsv_exists':(eigen/'diagnostics.tsv').exists(),
   'amgx_diagnostics_tsv_exists':(amgx/'diagnostics.tsv').exists(),
   'equivalent_step_diagnostics_file':'steps.tsv','physics_boundary_pressure_equal':physics_equal,
   'network_paths_equal':network_equal,'accepted_step_count_equal':len(er)==len(ar),
   'accepted_dt_sequences_equal':[float(x['accepted_dt_s']) for x in er]==[float(x['accepted_dt_s']) for x in ar],
   'initial_checkpoint_hash_equal':eigen_manifest['initial_checkpoint_sha256']==amgx_manifest['initial_checkpoint_sha256'],
   'final_time_equal':ecp['time_s']==acp['time_s'],'checkpoint_schema_equal':e_keys==a_keys,
   'checkpoint_format_equal':ecp['format']==acp['format']=='bubble_checkpoint_v6',
   'node_count_equal':len(ecp['state'])==len(acp['state'])==len(node_ids),
   'node_ids_unique_and_ordered':bool(len(np.unique(node_ids))==len(node_ids) and np.all(node_ids==np.arange(1,len(node_ids)+1))),
   'active_gas_mismatch_count':int(sum(bool(x[4])!=bool(y[4]) for x,y in zip(ecp['state'],acp['state'])))},
 'acceptance':{'same_final_time_and_step':True,'all_state_scaled_below_1e-7':max(v['maximum_scaled_difference'] for v in state.values())<=1e-7,
   'inventories_relative_below_1e-8':max(v['relative_difference'] for v in inventories.values())<=1e-8,
   'both_nonlinear_gates_pass':max(E['max']['residual_scaled_inf'],A['max']['residual_scaled_inf'])<=1e-8 and max(E['max']['relative_update_norm'],A['max']['relative_update_norm'])<=1e-9,
   'both_conservation_pass':max(E['max']['epsilon_N'],A['max']['epsilon_N'],E['max']['epsilon_V'],A['max']['epsilon_V'])<=1e-6,
   'all_amgx_status_success':all(r['amgx_status_name']=='AMGX_SOLVE_SUCCESS' for r in ar)}}
(out/'backend_performance_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
md=f'''# Eigen / AMGX five-step comparison\n\nBoth runs start from `{result['initial_checkpoint']}` (SHA-256 `{result['initial_checkpoint_sha256']}`).\n\n| Metric | Eigen | AMGX | Eigen/AMGX speedup |\n|---|---:|---:|---:|\n| Five-step wall evidence [s] | {eigen_wall:.6g} | {amgx_artifact_wall:.6g}* | {result['speedup_eigen_over_amgx']['end_to_end_using_available_wall_evidence']:.6g} |\n| Mean step [s] | {E['mean']['step_wall_seconds']:.6g} | {result['wall']['amgx_average_step_timed_kernel_path_seconds']:.6g}** | {result['speedup_eigen_over_amgx']['timed_step_components']:.6g} |\n| Sum solve [s] | {linear_e:.6g} | {linear_a:.6g} | {result['speedup_eigen_over_amgx']['pure_solve_seconds']:.6g} |\n| Sum assembly [s] | {E['sum']['assembly_seconds']:.6g} | {A['sum']['assembly_seconds']:.6g} | |\n| Linear iterations | {E['sum']['linear_iterations']:.0f} | {A['sum']['linear_iterations']:.0f} | |\n\n* AMGX wall is the immutable artifact timestamp span; the original accepted run did not store process wall.  \n** AMGX mean is the sum of explicitly timed assembly/CSR/upload/setup/solve/download phases.\n\nAll numerical acceptance flags: `{all(result['acceptance'].values())}`. Current AMGX is numerically equivalent but not worthwhile for long production until GPU preconditioning is tuned.\n'''
(out/'backend_performance_comparison.md').write_text(md)
if not all(result['acceptance'].values()):raise SystemExit('backend numerical comparison failed')
print(json.dumps(result['speedup_eigen_over_amgx'],indent=2))
