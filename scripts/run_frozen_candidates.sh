#!/usr/bin/env bash
set -euo pipefail
ROOT=/workspace/zz
PROJECT=$ROOT/bubble_dissolution_solver_amgx
BUILD=$ROOT/build-bubble-amgx-release
DIR=$PROJECT/output/verification/frozen_status3
PREFIX=$DIR/newton0
progress() { printf '{"stage":"%s","status":"%s"}\n' "$1" "$2" > "$DIR/progress.json.tmp"; mv "$DIR/progress.json.tmp" "$DIR/progress.json"; }
progress eigen running
$BUILD/solve_frozen_eigen "$PREFIX"_matrix_scaled.mtx "$PREFIX"_rhs_scaled.mtx "$DIR/eigen_result.json" > "$DIR/eigen_run.log" 2>&1
progress known_failure running
set +e
BUBBLE_AMGX_DUMP_PREFIX=$DIR/candidate_failure $BUILD/solve_frozen_amgx "$PREFIX"_matrix_scaled.mtx "$PREFIX"_rhs_scaled.mtx "$PROJECT/configs/amgx_fgmres_pnm_baseline.json" "$DIR/candidate_failure_wrapper.json" > "$DIR/candidate_failure.log" 2>&1
failure_rc=$?
set -e
if [[ $failure_rc -eq 0 ]]; then progress known_failure unexpected_success; exit 3; fi
progress selected_config running
BUBBLE_AMGX_DUMP_PREFIX=$DIR/candidate_selected $BUILD/solve_frozen_amgx "$PREFIX"_matrix_scaled.mtx "$PREFIX"_rhs_scaled.mtx "$PROJECT/configs/amgx_fgmres_block2_amg_nosolver.json" "$DIR/candidate_selected_wrapper.json" > "$DIR/candidate_selected.log" 2>&1
progress complete passed
