#!/usr/bin/env python3
import hashlib,json,subprocess,sys
from pathlib import Path

output=Path(sys.argv[1]);binary=Path(sys.argv[2]);small_test=Path(sys.argv[3]);config=Path(sys.argv[4])
def run(command):
    p=subprocess.run(command,text=True,capture_output=True);return {'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
nvidia=run(['nvidia-smi','--query-gpu=index,name,memory.total,memory.free,driver_version,compute_cap','--format=csv,noheader'])
cuda=run(['/usr/local/cuda/bin/nvcc','--version']);probe=run([str(small_test),str(config)])
result={'format':'bubble_amgx_environment_audit_v1','gpu_probe':nvidia,'cuda_probe':cuda,
 'amgx_runtime_probe':probe,'amgx':{'header':'/opt/amgx_install/include/amgx_c.h',
 'header_sha256':sha('/opt/amgx_install/include/amgx_c.h'),'library':'/opt/amgx_install/lib/libamgxsh.so',
 'library_sha256':sha('/opt/amgx_install/lib/libamgxsh.so'),'tested_config':str(config.resolve()),
 'tested_config_sha256':sha(config)},'binary':str(binary.resolve()),'binary_sha256':sha(binary),
 'all_probes_passed':nvidia['returncode']==0 and cuda['returncode']==0 and probe['returncode']==0}
output.write_text(json.dumps(result,indent=2)+'\n')
if not result['all_probes_passed']:raise SystemExit('AMGX environment audit failed')
