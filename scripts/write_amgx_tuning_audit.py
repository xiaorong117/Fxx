#!/usr/bin/env python3
import json,sys
from pathlib import Path

comparison=json.loads(Path(sys.argv[1]).read_text());output=Path(sys.argv[2])
e=comparison['eigen'];a=comparison['amgx']
e_linear=e['setup_seconds']+e['solve_seconds']
a_linear=sum(a[k] for k in ('csr_seconds','upload_seconds','setup_seconds','solve_seconds','download_seconds'))
result={'format':'bubble_amgx_tuning_audit_v1','selected':'block2_ruiz_fgmres_aggregation_amg_nosolver',
 'selected_configuration':'configs/amgx_fgmres_block2_amg_nosolver.json',
 'matrix':{'scalar_rows':comparison['dof_count'],'scalar_nonzeros':comparison['jacobian_nonzeros'],
           'maximum_scalar_column_index':comparison['csr_max_column_index'],'block_dimension':2,
           'index_base':0,'signed_int32_safe':comparison['csr_max_column_index']<2**31},
 'rejected_candidates':[
  {'name':'scalar_fgmres_aggregation_dilu','reason':'AMGX_SOLVE_NOT_CONVERGED status=3'},
  {'name':'scalar_pbicgstab_aggregation_jacobi','reason':'did not finish within 22 minutes'},
  {'name':'block2_fgmres_multicolor_dilu','reason':'converged but slower than selected AMG candidate'},
  {'name':'block2_fgmres_aggregation_amg_dense_lu','reason':'over 60 seconds without completing; NOSOLVER tested instead'}],
 'acceptance':comparison['acceptance'],'maximum_scaled_state_difference':comparison['max_scaled_state_difference'],
 'eigen_linear_seconds':e_linear,'amgx_linear_seconds':a_linear,
 'amgx_over_eigen_linear_time_ratio':a_linear/e_linear,
 'eigen_total_newton_seconds':e['assembly_seconds']+e_linear,
 'amgx_total_newton_seconds':a['assembly_seconds']+a_linear,
 'amgx_over_eigen_total_newton_time_ratio':(a['assembly_seconds']+a_linear)/(e['assembly_seconds']+e_linear),
 'conclusion':'Numerically accepted but no speedup on this checkpoint; production GPU tuning deferred.'}
output.write_text(json.dumps(result,indent=2)+'\n')
