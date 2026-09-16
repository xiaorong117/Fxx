cmake_minimum_required(VERSION 3.18)
project(bubble_dissolution_solver VERSION 0.1.0 LANGUAGES CXX)
find_package(Python3 COMPONENTS Interpreter REQUIRED)
find_package(ZLIB REQUIRED)

option(BUBBLE_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(BUBBLE_ENABLE_AMGX "Enable the NVIDIA AMGX double-precision linear backend" OFF)
option(BUBBLE_ENABLE_PARDISO "Enable the oneMKL PARDISO double-precision linear backend" OFF)
set(AMGX_ROOT "/opt/amgx_install" CACHE PATH "AMGX installation prefix")
set(PARDISO_MKL_ROOT "/opt/conda/envs/Mesh" CACHE PATH "oneMKL runtime prefix for PARDISO")
set(FADBAD_ROOT "" CACHE PATH "Directory containing badiff.h from FADBAD++ 2.1")
set(NLOHMANN_JSON_ROOT "" CACHE PATH "Prefix containing nlohmann/json.hpp")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_path(FADBAD_INCLUDE_DIR badiff.h
  HINTS "${FADBAD_ROOT}" ENV FADBAD_ROOT "/opt/FADBAD++-2.1/FADBAD++"
  PATH_SUFFIXES FADBAD++ FADBAD++-2.1/FADBAD++)
find_path(NLOHMANN_JSON_INCLUDE_DIR nlohmann/json.hpp
  HINTS "${NLOHMANN_JSON_ROOT}" ENV NLOHMANN_JSON_ROOT ENV CONDA_PREFIX ENV _CONDA_ROOT
  PATH_SUFFIXES include)

if(NOT FADBAD_INCLUDE_DIR)
  message(FATAL_ERROR "FADBAD++ badiff.h not found; set -DFADBAD_ROOT=/path/to/FADBAD++")
endif()
if(NOT NLOHMANN_JSON_INCLUDE_DIR)
  message(FATAL_ERROR "nlohmann/json.hpp not found")
endif()

function(header_tree_sha256 ROOT RELATIVE_GLOB OUTPUT_VARIABLE)
  file(GLOB_RECURSE TREE_FILES RELATIVE "${ROOT}" "${ROOT}/${RELATIVE_GLOB}")
  list(SORT TREE_FILES)
  set(TREE_MANIFEST "")
  foreach(RELATIVE_PATH IN LISTS TREE_FILES)
    if(NOT IS_DIRECTORY "${ROOT}/${RELATIVE_PATH}")
      file(SHA256 "${ROOT}/${RELATIVE_PATH}" FILE_HASH)
      string(APPEND TREE_MANIFEST "${RELATIVE_PATH}\t${FILE_HASH}\n")
    endif()
  endforeach()
  string(SHA256 TREE_HASH "${TREE_MANIFEST}")
  set(${OUTPUT_VARIABLE} "${TREE_HASH}" PARENT_SCOPE)
endfunction()

get_target_property(EIGEN_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
header_tree_sha256("${EIGEN_INCLUDE_DIR}" "Eigen/*" EIGEN_TREE_SHA256)
header_tree_sha256("${NLOHMANN_JSON_INCLUDE_DIR}" "nlohmann/*" NLOHMANN_TREE_SHA256)
header_tree_sha256("${FADBAD_INCLUDE_DIR}" "*" FADBAD_TREE_SHA256)
file(GLOB_RECURSE BUBBLE_PROJECT_SOURCE_FILES
  RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(SORT BUBBLE_PROJECT_SOURCE_FILES)
# Reconfigure when any hashed source changes so an incremental build cannot
# retain stale provenance in the executable.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  ${BUBBLE_PROJECT_SOURCE_FILES})
set(BUBBLE_PROJECT_SOURCE_MANIFEST "")
foreach(RELATIVE_PATH IN LISTS BUBBLE_PROJECT_SOURCE_FILES)
  file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/${RELATIVE_PATH}" FILE_HASH)
  string(APPEND BUBBLE_PROJECT_SOURCE_MANIFEST "${RELATIVE_PATH}\t${FILE_HASH}\n")
endforeach()
string(SHA256 BUBBLE_PROJECT_SOURCE_TREE_SHA256 "${BUBBLE_PROJECT_SOURCE_MANIFEST}")
message(STATUS "Eigen header-tree SHA-256: ${EIGEN_TREE_SHA256}")
message(STATUS "nlohmann/json header-tree SHA-256: ${NLOHMANN_TREE_SHA256}")
message(STATUS "FADBAD++ distribution-tree SHA-256: ${FADBAD_TREE_SHA256}")
message(STATUS "Project source-tree SHA-256: ${BUBBLE_PROJECT_SOURCE_TREE_SHA256}")

find_package(Git QUIET)
set(BUBBLE_GIT_COMMIT "unknown")
if(Git_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE BUBBLE_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
endif()

add_library(bubble_core
  src/config.cpp
  src/network.cpp
  src/solver.cpp
  src/output.cpp
  src/analysis.cpp
  src/amgx_backend.cpp
  src/pardiso_backend.cpp)
target_include_directories(bubble_core
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
  PUBLIC "${NLOHMANN_JSON_INCLUDE_DIR}")
target_include_directories(bubble_core SYSTEM PUBLIC "${FADBAD_INCLUDE_DIR}")
target_link_libraries(bubble_core PUBLIC Eigen3::Eigen ZLIB::ZLIB)
if(BUBBLE_ENABLE_AMGX)
  find_path(AMGX_INCLUDE_DIR amgx_c.h HINTS "${AMGX_ROOT}/include" REQUIRED)
  find_library(AMGX_LIBRARY NAMES amgxsh HINTS "${AMGX_ROOT}/lib" REQUIRED)
  target_include_directories(bubble_core SYSTEM PRIVATE "${AMGX_INCLUDE_DIR}")
  target_link_libraries(bubble_core PUBLIC "${AMGX_LIBRARY}" dl pthread)
  get_filename_component(AMGX_LIBRARY_REALPATH "${AMGX_LIBRARY}" REALPATH)
  set_target_properties(bubble_core PROPERTIES BUILD_RPATH "${AMGX_ROOT}/lib;/usr/local/cuda/lib64")
else()
  set(AMGX_LIBRARY_REALPATH "unavailable")
endif()

if(BUBBLE_ENABLE_PARDISO)
  find_library(PARDISO_MKL_LIBRARY
    NAMES mkl_rt libmkl_rt.so.2
    HINTS "${PARDISO_MKL_ROOT}/lib"
    REQUIRED)
  target_link_libraries(bubble_core PUBLIC dl pthread)
  get_filename_component(PARDISO_MKL_LIBRARY_REALPATH "${PARDISO_MKL_LIBRARY}" REALPATH)
else()
  set(PARDISO_MKL_LIBRARY_REALPATH "unavailable")
endif()
target_compile_options(bubble_core PRIVATE -Wall -Wextra -Wpedantic)
target_compile_definitions(bubble_core PUBLIC
  BUBBLE_GIT_COMMIT="${BUBBLE_GIT_COMMIT}"
  BUBBLE_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}"
  BUBBLE_COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}"
  BUBBLE_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
  BUBBLE_SANITIZERS_ENABLED=$<BOOL:${BUBBLE_ENABLE_SANITIZERS}>
  BUBBLE_AMGX_ENABLED=$<BOOL:${BUBBLE_ENABLE_AMGX}>
  BUBBLE_AMGX_LIBRARY_PATH="${AMGX_LIBRARY_REALPATH}"
  BUBBLE_PARDISO_ENABLED=$<BOOL:${BUBBLE_ENABLE_PARDISO}>
  BUBBLE_PARDISO_LIBRARY_PATH="${PARDISO_MKL_LIBRARY_REALPATH}"
  BUBBLE_PROJECT_SOURCE_TREE_SHA256="${BUBBLE_PROJECT_SOURCE_TREE_SHA256}"
  BUBBLE_FADBAD_INCLUDE_DIR="${FADBAD_INCLUDE_DIR}"
  BUBBLE_EIGEN_VERSION="${EIGEN3_VERSION_STRING}"
  BUBBLE_EIGEN_TREE_SHA256="${EIGEN_TREE_SHA256}"
  BUBBLE_NLOHMANN_TREE_SHA256="${NLOHMANN_TREE_SHA256}"
  BUBBLE_FADBAD_TREE_SHA256="${FADBAD_TREE_SHA256}")

if(BUBBLE_ENABLE_SANITIZERS)
  target_compile_options(bubble_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(bubble_core PUBLIC -fsanitize=address,undefined)
endif()

add_executable(bubble_solver src/main.cpp)
target_link_libraries(bubble_solver PRIVATE bubble_core)
target_compile_options(bubble_solver PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
# Keep historical numerical outputs immutable. Each build owns regression output.
foreach(CASE_NAME small_cli small_virtual_cli)
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/${CASE_NAME}.json" CASE_JSON)
  string(JSON CASE_JSON SET "${CASE_JSON}" io output_directory
    "\"${CMAKE_CURRENT_BINARY_DIR}/${CASE_NAME}\"")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${CASE_NAME}.json" "${CASE_JSON}")
endforeach()
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/small_relative_visualization.json"
  RELATIVE_VISUALIZATION_JSON)
string(JSON RELATIVE_VISUALIZATION_JSON SET "${RELATIVE_VISUALIZATION_JSON}"
  io output_directory "\"${CMAKE_CURRENT_BINARY_DIR}/relative_visualization_test\"")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/small_relative_visualization.json"
  "${RELATIVE_VISUALIZATION_JSON}")
if(NOT BUBBLE_ENABLE_PARDISO)
  set(PARDISO_DISABLED_JSON "${RELATIVE_VISUALIZATION_JSON}")
  string(JSON PARDISO_DISABLED_JSON SET "${PARDISO_DISABLED_JSON}"
    io output_directory "\"${CMAKE_CURRENT_BINARY_DIR}/pardiso_disabled\"")
  string(JSON PARDISO_DISABLED_JSON SET "${PARDISO_DISABLED_JSON}"
    linear backend "\"pardiso\"")
  string(JSON PARDISO_DISABLED_JSON SET "${PARDISO_DISABLED_JSON}"
    linear pardiso_threads 16)
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pardiso_disabled.json"
    "${PARDISO_DISABLED_JSON}")
  add_test(NAME pardiso_unavailable_rejected
    COMMAND bubble_solver --config
      "${CMAKE_CURRENT_BINARY_DIR}/pardiso_disabled.json")
  set_tests_properties(pardiso_unavailable_rejected PROPERTIES WILL_FAIL TRUE)
endif()
if(BUBBLE_ENABLE_AMGX)
  set(SMALL_AMGX_JSON "${RELATIVE_VISUALIZATION_JSON}")
  string(JSON SMALL_AMGX_JSON SET "${SMALL_AMGX_JSON}" io output_directory
    "\"${CMAKE_CURRENT_BINARY_DIR}/small_amgx_scientific\"")
  string(JSON SMALL_AMGX_JSON SET "${SMALL_AMGX_JSON}" io state_format "\"compressed_vtp\"")
  string(JSON SMALL_AMGX_JSON SET "${SMALL_AMGX_JSON}" linear backend "\"amgx\"")
  string(JSON SMALL_AMGX_JSON SET "${SMALL_AMGX_JSON}" linear amgx_config_file
    "\"${CMAKE_CURRENT_SOURCE_DIR}/configs/amgx_fgmres_block2_amg_nosolver.json\"")
  string(JSON SMALL_AMGX_JSON SET "${SMALL_AMGX_JSON}" termination
    "{\"primary\":\"main_real_gas_saturation\",\"target_gas_saturation\":0.01,\"maximum_physical_time_s\":0.0005,\"maximum_wall_seconds\":0,\"maximum_output_bytes\":100000000}")
  string(JSON SMALL_AMGX_JSON SET "${SMALL_AMGX_JSON}" analysis
    "{\"x_segments\":20,\"vtk_output_pore_volumes\":[]}")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/small_amgx_scientific.json" "${SMALL_AMGX_JSON}")
endif()
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/configs/bubble_solver.smoke.json" SMOKE_JSON)
string(JSON SMOKE_JSON SET "${SMOKE_JSON}" io output_directory
  "\"${CMAKE_CURRENT_BINARY_DIR}/full_network_regression\"")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/smoke.json" "${SMOKE_JSON}")
add_executable(test_ad_smoke tests/test_ad_smoke.cpp)
target_link_libraries(test_ad_smoke PRIVATE bubble_core)
target_compile_options(test_ad_smoke PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME ad_smoke COMMAND test_ad_smoke)

add_executable(test_bubble_solver tests/test_bubble_solver.cpp)
target_link_libraries(test_bubble_solver PRIVATE bubble_core)
target_compile_options(test_bubble_solver PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME bubble_small_cases COMMAND test_bubble_solver)
add_test(NAME backend_comparison_schema_regression COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_backend_comparison_schema.py")
add_test(NAME adaptive_dt_small_regression COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_adaptive_dt.py" "$<TARGET_FILE:bubble_solver>")
set_tests_properties(adaptive_dt_small_regression PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/..")
add_test(NAME adaptive_dt_end_time_roundoff COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_end_time_roundoff.py")
add_test(NAME accelerated_preview_small_regression COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_accelerated_preview.py"
  "$<TARGET_FILE:bubble_solver>" "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME case_dashboard_provenance COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_case_dashboard.py"
  "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME case_dashboard_control COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_case_control.py"
  "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME case_dashboard_live_console COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_dashboard_console.py"
  "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME delivery_documentation COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_delivery_documentation.py"
  "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME boundaryless_component_exclusion COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_exclude_boundaryless_component.py"
  "${CMAKE_CURRENT_SOURCE_DIR}")

if(BUBBLE_ENABLE_AMGX)
  add_executable(test_amgx_backend tests/test_amgx_backend.cpp)
  target_link_libraries(test_amgx_backend PRIVATE bubble_core)
  target_compile_options(test_amgx_backend PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME amgx_small_nonsymmetric_system COMMAND test_amgx_backend
    "${CMAKE_CURRENT_SOURCE_DIR}/configs/amgx_fgmres_aggregation_dilu.json")
  set_tests_properties(amgx_small_nonsymmetric_system PROPERTIES LABELS GPU)
  add_executable(compare_full_checkpoint_backends tests/compare_full_checkpoint_backends.cpp)
  target_link_libraries(compare_full_checkpoint_backends PRIVATE bubble_core)
  target_compile_options(compare_full_checkpoint_backends PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(solve_frozen_eigen tests/solve_frozen_eigen.cpp)
  target_link_libraries(solve_frozen_eigen PRIVATE bubble_core)
  target_compile_options(solve_frozen_eigen PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(solve_frozen_amgx tests/solve_frozen_amgx.cpp)
  target_link_libraries(solve_frozen_amgx PRIVATE bubble_core)
  target_compile_options(solve_frozen_amgx PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(run_amgx_checkpoint_sequence tests/run_amgx_checkpoint_sequence.cpp)
  target_link_libraries(run_amgx_checkpoint_sequence PRIVATE bubble_core)
  target_compile_options(run_amgx_checkpoint_sequence PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(retag_checkpoint_time_control tests/retag_checkpoint_time_control.cpp)
target_link_libraries(retag_checkpoint_time_control PRIVATE bubble_core)
target_compile_options(retag_checkpoint_time_control PRIVATE -Wall -Wextra -Wpedantic)

if(BUBBLE_ENABLE_PARDISO)
  add_executable(test_pardiso_backend tests/test_pardiso_backend.cpp)
  target_link_libraries(test_pardiso_backend PRIVATE bubble_core)
  target_compile_options(test_pardiso_backend PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME pardiso_small_nonsymmetric_reuse COMMAND test_pardiso_backend)
  set_tests_properties(pardiso_small_nonsymmetric_reuse PROPERTIES LABELS CPU_DIRECT)
  add_executable(solve_frozen_pardiso tests/solve_frozen_pardiso.cpp)
  target_link_libraries(solve_frozen_pardiso PRIVATE bubble_core)
  target_compile_options(solve_frozen_pardiso PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(compare_full_checkpoint_pardiso
    tests/compare_full_checkpoint_pardiso.cpp)
  target_link_libraries(compare_full_checkpoint_pardiso PRIVATE bubble_core)
  target_compile_options(compare_full_checkpoint_pardiso PRIVATE
    -Wall -Wextra -Wpedantic)
  add_executable(run_pardiso_checkpoint_sequence
    tests/run_pardiso_checkpoint_sequence.cpp)
  target_link_libraries(run_pardiso_checkpoint_sequence PRIVATE bubble_core)
  target_compile_options(run_pardiso_checkpoint_sequence PRIVATE
    -Wall -Wextra -Wpedantic)
endif()

add_test(NAME research_null_config_rejected
  COMMAND bubble_solver --config
          "${CMAKE_CURRENT_SOURCE_DIR}/../remote_codex_bubble_case/bubble_solver.research-template.json")
set_tests_properties(research_null_config_rejected PROPERTIES WILL_FAIL TRUE)

add_test(NAME full_network_degenerate_geometry_audit
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/scripts/audit_strict_geometry.py"
          "${CMAKE_CURRENT_SOURCE_DIR}/.."
          --pore-file PNM_new_pore_virtual.txt
          --throat-file PNM_new_throat_virtual.txt)

add_test(NAME physical_units_and_virtual_boundary_generation
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_physical_units_virtual_boundaries.py"
          "${CMAKE_CURRENT_SOURCE_DIR}"
          "${CMAKE_CURRENT_SOURCE_DIR}/..")
set_tests_properties(full_network_degenerate_geometry_audit PROPERTIES
  DEPENDS physical_units_and_virtual_boundary_generation)

add_test(NAME batch35_geometry_and_flow
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_batch35_geometry_and_flow.py"
          "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME full_network_five_step_smoke
  COMMAND bubble_solver --config
          "${CMAKE_CURRENT_BINARY_DIR}/smoke.json")
set_tests_properties(full_network_five_step_smoke PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.."
  DEPENDS full_network_degenerate_geometry_audit
  LABELS FULL_NETWORK)
add_test(NAME full_network_output_audit
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/audit_full_network_output.py"
          "${CMAKE_CURRENT_SOURCE_DIR}/.."
          "${CMAKE_CURRENT_BINARY_DIR}/full_network_regression")
set_tests_properties(full_network_output_audit PROPERTIES
  DEPENDS full_network_five_step_smoke
  LABELS FULL_NETWORK)

add_test(NAME cli_output_integration
  COMMAND bubble_solver --config
          "${CMAKE_CURRENT_BINARY_DIR}/small_cli.json")
set_tests_properties(cli_output_integration PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/..")
add_test(NAME cli_output_audit
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/audit_cli_output.py"
          "${CMAKE_CURRENT_BINARY_DIR}/small_cli")
set_tests_properties(cli_output_audit PROPERTIES DEPENDS cli_output_integration)

add_test(NAME virtual_boundary_cli_integration
  COMMAND bubble_solver --config
          "${CMAKE_CURRENT_BINARY_DIR}/small_virtual_cli.json")
set_tests_properties(virtual_boundary_cli_integration PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/..")
add_test(NAME virtual_boundary_cli_audit
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/audit_virtual_cli_output.py"
          "${CMAKE_CURRENT_BINARY_DIR}/small_virtual_cli")
set_tests_properties(virtual_boundary_cli_audit PROPERTIES
  DEPENDS virtual_boundary_cli_integration)

add_test(NAME cli_restart_adaptive_equivalence
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli_restart.py"
          "$<TARGET_FILE:bubble_solver>" "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME unsafe_configs_rejected_before_mesh
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_invalid_configs.py"
          "$<TARGET_FILE:bubble_solver>" "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME failure_manifest_preserves_progress
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_failure_progress_manifest.py"
          "$<TARGET_FILE:bubble_solver>" "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME optional_geometry_correction_template_untouched
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/audit_correction_template.py"
          "${CMAKE_CURRENT_SOURCE_DIR}/.."
          "${CMAKE_CURRENT_SOURCE_DIR}/configs/geometry_corrections.required.json")

add_test(NAME geometry_correction_overlay_end_to_end
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_geometry_overlay.py"
          "${CMAKE_CURRENT_SOURCE_DIR}/scripts/geometry_correction_overlay.py")

add_test(NAME relative_visualization_integration
  COMMAND bubble_solver --config
          "${CMAKE_CURRENT_BINARY_DIR}/small_relative_visualization.json")
set_tests_properties(relative_visualization_integration PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/..")
add_test(NAME relative_visualization_postprocess
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/scripts/relative_visualize.py"
          "${CMAKE_CURRENT_BINARY_DIR}/relative_visualization_test")
set_tests_properties(relative_visualization_postprocess PROPERTIES
  DEPENDS relative_visualization_integration)
add_test(NAME relative_visualization_readback_audit
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/audit_relative_visualization.py"
          "${CMAKE_CURRENT_BINARY_DIR}/relative_visualization_test"
          "${CMAKE_CURRENT_SOURCE_DIR}/scripts/relative_visualize.py")
set_tests_properties(relative_visualization_readback_audit PROPERTIES
  DEPENDS relative_visualization_postprocess)

if(BUBBLE_ENABLE_AMGX)
  add_test(NAME amgx_scientific_output_integration COMMAND bubble_solver --config
    "${CMAKE_CURRENT_BINARY_DIR}/small_amgx_scientific.json")
  set_tests_properties(amgx_scientific_output_integration PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.." LABELS GPU)
  add_test(NAME amgx_scientific_output_audit COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/audit_scientific_outputs.py"
    "${CMAKE_CURRENT_BINARY_DIR}/small_amgx_scientific" 20 5)
  set_tests_properties(amgx_scientific_output_audit PROPERTIES
    DEPENDS amgx_scientific_output_integration LABELS GPU)
endif()
