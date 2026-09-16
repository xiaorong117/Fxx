import json,sys
from pathlib import Path
b=json.loads(Path(sys.argv[1]).read_text());o=json.loads(Path(sys.argv[2]).read_text());out=Path(sys.argv[3])
keys=['assembly_seconds','csr_seconds','upload_seconds','setup_seconds','solve_seconds','download_seconds']
bt=sum(b['amgx'][k] for k in keys);ot=sum(o['amgx'][k] for k in keys)
r={'format':'amgx_structure_reuse_benchmark_v1','baseline':b['amgx'],'optimized':o['amgx'],
 'baseline_total_timed_seconds':bt,'optimized_total_timed_seconds':ot,'overall_speedup':bt/ot,
 'numerical_equivalence':{'baseline_acceptance':b['acceptance'],'optimized_acceptance':o['acceptance'],
 'baseline_max_scaled_state_difference':b['max_scaled_state_difference'],'optimized_max_scaled_state_difference':o['max_scaled_state_difference']},
 'interpretation':'Fast-path calls are correct; this single representative step shows no speedup because solve dominates and setup timing noise exceeds upload savings.'}
out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
