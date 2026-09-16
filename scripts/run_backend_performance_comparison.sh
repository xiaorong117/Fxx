#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/workspace/zz
PROJECT=$ROOT/bubble_dissolution_solver_amgx
BUILD=$ROOT/build-bubble-amgx-release
VERIFY=$PROJECT/output/verification
AMGX_RUN=$VERIFY/amgx_checkpoint_five_steps
EIGEN_RUN=$VERIFY/eigen_checkpoint_five_steps
RESULT=$VERIFY/backend_comparison
CHECKPOINT=$ROOT/bubble_dissolution_solver/output/full_network_air_water_relative_pressure_run01/checkpoint.json
BASE=$PROJECT/configs/bubble_solver.air_water_relative_pressure.json
AMGX_CONFIG=$PROJECT/configs/amgx_fgmres_block2_amg_nosolver.json
PROGRESS=$RESULT/progress.json
mkdir -p "$RESULT"
update(){ printf '{"stage":"%s","status":"%s","updated_utc":"%s"}\n' "$1" "$2" "$(date -u +%FT%TZ)" > "$PROGRESS.tmp"; mv "$PROGRESS.tmp" "$PROGRESS"; }
failed(){ code=$?;update "${stage:-unknown}" failed;exit "$code";};trap failed ERR
stage=build;update "$stage" running
cmake --build "$BUILD" -j2 > "$RESULT/build.log" 2>&1
stage=eigen_five_steps;update "$stage" running
"$BUILD/run_amgx_checkpoint_sequence" "$BASE" "$CHECKPOINT" "$AMGX_CONFIG" "$EIGEN_RUN" 5 0.01 eigen > "$RESULT/eigen_five_steps.log" 2>&1
stage=compare_outputs;update "$stage" running
python3 "$PROJECT/scripts/compare_backend_five_steps.py" "$EIGEN_RUN" "$AMGX_RUN" "$RESULT" > "$RESULT/comparison.log" 2>&1
stage=related_release_tests;update "$stage" running
ctest --test-dir "$BUILD" --output-on-failure -R '^(ad_smoke|bubble_small_cases|amgx_small_nonsymmetric_system|cli_restart_adaptive_equivalence|amgx_scientific_output_integration|amgx_scientific_output_audit)$' --output-log "$RESULT/ctest-release.log" --output-junit "$RESULT/ctest-release.xml" > "$RESULT/ctest-console.log" 2>&1
stage=complete;update "$stage" passed
