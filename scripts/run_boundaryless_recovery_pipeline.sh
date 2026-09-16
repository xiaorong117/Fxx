#!/usr/bin/env bash
set -euo pipefail

project=/workspace/zz/bubble_dissolution_solver_amgx
binary=/workspace/zz/build-bubble-batch-retirement-release/bubble_solver
verify="$project/output/verification/boundaryless_component_recovery_85A15D_p08"
configs="$project/configs/boundaryless_component_recovery_85A15D_p08_v3"
production="$project/output/production/batch35_flow_pe555_pardiso/sands-of-85A15D-mu-0.5-100kpa-e_p08_boundaryless_frozen_recovery"
progress="$verify/progress.json"

write_stage() {
  python3 - "$progress" "$1" "$2" <<'PY'
import datetime,json,sys
path,stage,status=sys.argv[1:]
with open(path+'.tmp','w') as stream:
    json.dump({'stage':stage,'status':status,
               'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},stream,indent=2)
    stream.write('\n')
import os;os.replace(path+'.tmp',path)
PY
}

cd "$project"
write_stage equivalence running
if [[ ! -f "$verify/equivalence_report_t0p019.json" ]]; then
  while tmux has-session -t boundaryless_eq_original 2>/dev/null || \
        tmux has-session -t boundaryless_eq_filtered 2>/dev/null; do
    sleep 15
  done
  python3 scripts/compare_boundaryless_recovery.py \
    --original "$verify/original_equivalence_t0p02" \
    --filtered "$verify/filtered_equivalence_t0p02" \
    --mapping "$verify/input/old_to_new_node_id.tsv" \
    --report "$verify/equivalence_report_t0p019.json" \
    > "$verify/equivalence_compare_t0p019.log"
fi
python3 - "$verify/equivalence_report_t0p019.json" <<'PY'
import json,sys
if json.load(open(sys.argv[1])).get('status') != 'passed':
    raise SystemExit('boundaryless equivalence did not pass')
PY
write_stage cross_original_failure_time running

cross="$verify/filtered_cross_failure_pilot_t0p6"
test ! -e "$cross"
mkdir -p "$cross"
"$binary" --config "$configs/85A15D_p08_filtered_cross_failure_pilot_t0p6.json" \
  >> "$cross/driver.log" 2>&1
python3 - "$cross" <<'PY'
import csv,json,math,sys
from pathlib import Path
root=Path(sys.argv[1])
manifest=json.load(open(root/'run_manifest.json'))
term=json.load(open(root/'termination_events.json'))
with open(root/'diagnostics.tsv',newline='') as stream:
    rows=list(csv.DictReader(stream,delimiter='\t'))
if manifest.get('exit_status') != 0 or not math.isclose(float(rows[-1]['time_s']),0.6,abs_tol=1e-14):
    raise SystemExit('cross-failure pilot did not reach t=0.6')
if max(abs(float(r['epsilon_N_adjusted'])) for r in rows)>1e-6 or \
   max(abs(float(r['epsilon_V_adjusted'])) for r in rows)>1e-6:
    raise SystemExit('cross-failure pilot violated conservation')
if any(r['residual_criterion_passed']!='1' or r['update_criterion_passed']!='1' for r in rows):
    raise SystemExit('cross-failure pilot violated Newton acceptance')
report={'format':'bubble_boundaryless_cross_failure_pilot_v1','status':'passed',
        'accepted_steps':len(rows),'final_time_s':float(rows[-1]['time_s']),
        'termination_reason':term.get('termination_reason'),
        'minimum_liquid_relative_pressure_Pa':min(float(r['internal_pressure_min_Pa']) for r in rows),
        'minimum_liquid_absolute_pressure_Pa':101325+min(float(r['internal_pressure_min_Pa']) for r in rows),
        'maximum_epsilon_N':max(abs(float(r['epsilon_N_adjusted'])) for r in rows),
        'maximum_epsilon_V':max(abs(float(r['epsilon_V_adjusted'])) for r in rows),
        'total_retries':sum(int(float(r['timestep_retries'])) for r in rows)}
(root/'cross_failure_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
PY
python3 tests/test_exclude_boundaryless_component.py "$project" \
  > "$verify/filter_regression.log"
write_stage production_launch running

test ! -e "$production"
mkdir -p "$production"
sha256sum "$binary" \
  "$configs/85A15D_p08_boundaryless_frozen_production.json" \
  "$verify/input/PNM_new_pore_virtual_dynamic.txt" \
  "$verify/input/PNM_new_throat_virtual_dynamic.txt" \
  "$verify/input/PNM_new_connect_virtual_dynamic.txt" \
  > "$production/launch_sha256.txt"
tmux new-session -d -s boundaryless_recovery_production_85A15D_p08 \
  "cd '$project' && exec flock -n '$production/run.lock' '$binary' --config '$configs/85A15D_p08_boundaryless_frozen_production.json' >> '$production/driver.log' 2>&1"
write_stage production_running passed
