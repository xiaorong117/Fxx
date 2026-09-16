#!/usr/bin/env bash
set -euo pipefail

project=/workspace/zz/bubble_dissolution_solver_amgx
input_root="$project/output/batch35_inputs"
production_configs="$project/configs/batch35_flow_pe555_pardiso_v3"
production_root="$project/output/production/batch35_flow_pe555_pardiso"
pilot_configs="$project/configs/batch35_flow_pe555_pardiso_pilots_v3"
pilot_root="$project/output/verification/batch35_flow_pe555_pardiso_pilots_v3"
binary=/workspace/zz/build-bubble-batch-retirement-release/bubble_solver
log="$project/output/batch35_pipeline.log"

cd "$project"
while tmux has-session -t batch35_input_prepare_1200x7p5 2>/dev/null; do
  sleep 15
done
python3 - "$input_root/progress.json" <<'PY'
import json, sys
d=json.load(open(sys.argv[1]))
if d.get('status') != 'passed' or d.get('completed') != 70:
    raise SystemExit(f"input preparation did not pass: {d}")
PY

if [[ ! -f "$production_configs/batch_cases.json" ]]; then
  python3 scripts/generate_batch35_configs.py \
    --inventory "$input_root/batch_input_inventory.json" \
    --template "$project/configs/bubble_solver.validation_batch_microbubble_Q1e-8_pardiso.json" \
    --config-root "$production_configs" \
    --output-root "$production_root" \
    --pardiso-threads 8 \
    --maximum-output-bytes 2147483648 \
    --t-end-s 1000000
fi

python3 scripts/make_batch35_pilots.py \
  --index "$production_configs/batch_cases.json" \
  --config-root "$pilot_configs" \
  --output-root "$pilot_root" \
  --run-id sands-of-100C-mu-0.1-100kpa-e-0_p06 \
  --run-id sands-of-25ABCD-mu-0.5-100kpa-e_p08 \
  --t-end-s 0.002

python3 scripts/run_batch35.py \
  --index "$pilot_configs/pilot_cases.json" \
  --binary "$binary" --jobs 2 --threads-per-job 8

python3 scripts/audit_batch35_pilots.py \
  --index "$pilot_configs/pilot_cases.json" \
  --report "$pilot_root/pilot_audit.json"

mkdir -p "$production_root"
if tmux has-session -t batch35_production_p06_p08 2>/dev/null; then
  echo "refusing to duplicate existing batch35_production_p06_p08" >&2
  exit 1
fi
tmux new-session -d -s batch35_production_p06_p08 \
  "cd '$project' && exec python3 scripts/run_batch35.py --index '$production_configs/batch_cases.json' --binary '$binary' --jobs 8 --threads-per-job 8 >> '$production_root/scheduler.log' 2>&1"
echo "production scheduler launched after two passing t=0 pilots" >> "$log"
