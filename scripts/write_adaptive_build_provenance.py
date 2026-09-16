import hashlib,json,subprocess,sys
from pathlib import Path
project=Path(sys.argv[1]);binary=Path(sys.argv[2]);output=Path(sys.argv[3])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
files=sorted([project/'CMakeLists.txt',*sorted((project/'include').rglob('*.hpp')),*sorted((project/'src').glob('*.cpp'))])
manifest=''.join(f'{p.relative_to(project)}\t{sha(p)}\n' for p in files)
tree=hashlib.sha256(manifest.encode()).hexdigest();strings=subprocess.run(['strings',str(binary)],text=True,capture_output=True).stdout
tokens=['time_step_control','easy_newton_growth','normal_newton_keep','hard_newton_shrink','gas_timescale_cap','concentration_timescale_cap','gas_disappearance_event','dt_min_failure','linear_solver_retry','nonlinear_retry']
latest=max(p.stat().st_mtime_ns for p in files)
result={'format':'adaptive_dt_build_provenance_v1','source_tree_sha256':tree,'source_file_count':len(files),
 'source_latest_mtime_ns':latest,'binary':str(binary.resolve()),'binary_sha256':sha(binary),
 'binary_mtime_ns':binary.stat().st_mtime_ns,'binary_built_after_latest_source':binary.stat().st_mtime_ns>latest,
 'adaptive_symbols_or_strings':{t:(t in strings) for t in tokens},
 'cmake_build_dir':'/workspace/zz/build-bubble-amgx-release','compiler_definitions':'BUBBLE_AMGX_ENABLED=1;BUBBLE_SANITIZERS_ENABLED=0'}
output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
