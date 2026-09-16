#!/usr/bin/env python3
import json,sys
from pathlib import Path
import numpy as np
from scipy.io import mmread

prefix=Path(sys.argv[1]);output=Path(sys.argv[2])
A=mmread(str(prefix)+'_matrix_scaled.mtx').tocsr();J=mmread(str(prefix)+'_jacobian_unscaled.mtx').tocsr()
b=np.asarray(mmread(str(prefix)+'_rhs_scaled.mtx')).reshape(-1)
x0=np.asarray(mmread(str(prefix)+'_initial_solution.mtx')).reshape(-1)
x=np.asarray(mmread(str(prefix)+'_amgx_solution.mtx')).reshape(-1)
row_n=np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).reshape(-1));col_n=np.sqrt(np.asarray(A.multiply(A).sum(axis=0)).reshape(-1))
diag=A.diagonal();result=json.loads(Path(str(prefix)+'_amgx_result.json').read_text())
cpu=lambda y:float(np.linalg.norm(A@y-b)/max(np.linalg.norm(b),1e-300))
audit={'format':'bubble_frozen_linear_system_audit_v1','shape':list(A.shape),'nnz':int(A.nnz),
 'unscaled_jacobian_nnz':int(J.nnz),'finite_values':bool(np.isfinite(A.data).all() and np.isfinite(b).all()),
 'empty_rows':np.flatnonzero(np.diff(A.indptr)==0).tolist(),'zero_diagonal_count':int(np.count_nonzero(diag==0)),
 'tiny_diagonal_below_1e-14_count':int(np.count_nonzero(np.abs(diag)<1e-14)),
 'absolute_diagonal_range':[float(np.min(np.abs(diag))),float(np.max(np.abs(diag)))],
 'nonzero_absolute_value_range':[float(np.min(np.abs(A.data[A.data!=0]))),float(np.max(np.abs(A.data)))],
 'row_l2_norm_range':[float(row_n.min()),float(row_n.max())],
 'column_l2_norm_range':[float(col_n.min()),float(col_n.max())],
 'matrix_market_index_base':1,'runtime_csr_index_base':0,'runtime_integer_width_bits':32,
 'minimum_runtime_column':int(A.indices.min()),'maximum_runtime_column':int(A.indices.max()),
 'row_pointer_monotone':bool(np.all(np.diff(A.indptr)>=0)),'initial_solution_relative_residual':cpu(x0),
 'amgx_solution_relative_residual':cpu(x),'amgx_result':result,
 'classification':'preconditioner_not_converged' if result['status_name']=='AMGX_SOLVE_NOT_CONVERGED' else 'other'}
output.write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
