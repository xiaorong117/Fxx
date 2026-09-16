#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/workspace/zz
PROJECT=$ROOT/bubble_dissolution_solver_amgx
BUILD=$ROOT/build-bubble-amgx-release
VERIFY=$PROJECT/output/verification
PROGRESS=$VERIFY/amgx_validation_progress.json
CHECKPOINT=$ROOT/bubble_dissolution_solver/output/full_network_air_water_relative_pressure_run01/checkpoint.json
BASE=$PROJECT/configs/bubble_solver.air_water_relative_pressure.json
AMGX=$PROJECT/configs/amgx_fgmres_block2_amg_nosolver.json
SEQUENCE=$VERIFY/amgx_checkpoint_five_steps
update() { printf '{"stage":"%s","status":"%s","updated_utc":"%s"}\n' "$1" "$2" "$(date -u +%FT%TZ)" > "$PROGRESS.tmp"; mv "$PROGRESS.tmp" "$PROGRESS"; }
failed() { code=$?; update "${stage:-unknown}" failed; exit "$code"; }
trap failed ERR
stage=build; update "$stage" running
cmake --build "$BUILD" -j2 > "$VERIFY/post_frozen_build.log" 2>&1
stage=full_newton_backend_comparison; update "$stage" running
"$BUILD/compare_full_checkpoint_backends" "$BASE" "$CHECKPOINT" "$AMGX" "$VERIFY/post_frozen_full_step_comparison.json" 0.01 > "$VERIFY/post_frozen_full_step_comparison.log" 2>&1
stage=amgx_five_accepted_steps; update "$stage" running
"$BUILD/run_amgx_checkpoint_sequence" "$BASE" "$CHECKPOINT" "$AMGX" "$SEQUENCE" 5 0.01 > "$VERIFY/amgx_checkpoint_five_steps.log" 2>&1
stage=related_release_tests; update "$stage" running
ctest --test-dir "$BUILD" --output-on-failure -R '^(ad_smoke|bubble_small_cases|amgx_small_nonsymmetric_system|cli_restart_adaptive_equivalence|amgx_scientific_output_integration|amgx_scientific_output_audit)$' --output-log "$VERIFY/ctest-post-frozen-release.log" --output-junit "$VERIFY/ctest-post-frozen-release.xml"
stage=complete; update "$stage" passed
