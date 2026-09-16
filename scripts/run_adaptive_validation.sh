#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/workspace/zz
cd "$ROOT"
P=$ROOT/bubble_dissolution_solver_amgx
B=$ROOT/build-bubble-amgx-release
V=$P/output/verification/adaptive_dt
mkdir -p "$V"
progress(){ printf '{"stage":"%s","status":"%s","updated_utc":"%s"}\n' "$1" "$2" "$(date -u +%FT%TZ)" > "$V/progress.json.tmp";mv "$V/progress.json.tmp" "$V/progress.json"; }
trap 'progress "${stage:-unknown}" failed' ERR
stage=build;progress "$stage" running
cmake --build "$B" -j2 > "$V/build.log" 2>&1
START="$V/common_start_checkpoint.json"
"$B/retag_checkpoint_time_control" "$ROOT/bubble_dissolution_solver/output/full_network_air_water_relative_pressure_run01/checkpoint.json" "$P/configs/bubble_solver.fixed_dt_validation.json" "$START" 0.01
stage=fixed_dt;progress "$stage" running
"$B/bubble_solver" --config "$P/configs/bubble_solver.fixed_dt_validation.json" --restart "$START" > "$V/fixed_dt.log" 2>&1
stage=adaptive_dt;progress "$stage" running
"$B/bubble_solver" --config "$P/configs/bubble_solver.adaptive_dt_validation.json" --restart "$START" > "$V/adaptive_dt.log" 2>&1
stage=compare;progress "$stage" running
python3 "$P/scripts/compare_adaptive_fixed.py" "$V/fixed_final4" "$V/adaptive_final4" "$V/compare_adaptive_fixed.json" > "$V/compare_adaptive_fixed.log" 2>&1
stage=complete;progress "$stage" passed
