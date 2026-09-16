# Status

## Implemented

- Independent C++17 project under `bubble_dissolution_solver/`; `../trans-react.cpp` is reference-only and has not been modified.
- Authority audit completed for the Word document and package notes. Equations (29)--(35) agree with `02_控制方程与数值规格.md`: `Rl` is liquid-volume accumulation plus outward liquid flux; `Rd` is dissolved inventory accumulation plus outward advective/diffusive flux minus transfer; `Rg` is free-gas inventory accumulation plus the identical transfer; `Rc=Pg-Pl-Pc`; all are at `n+1`. SI units are respectively m3/s, mol/s, mol/s, Pa. Boundary flag 1 supplies absolute inlet `Pl,C`; flag 2 supplies absolute outlet `Pl`, zero external diffusion, natural outflow and configured backflow concentration.
- The old framework was audited only for its `fadbad::B<T>`/`diff`/`d` reverse-mode pattern, local COO assembly, CSR handoff, Eigen and AMGX interface structure. It solves a different methane/CO2 adsorption/transport problem and cannot be directly modified into this bubble case. None of its hard-coded boundary arrays, reactions, adsorption, GSL, CUDA helpers, or fixed pore counts were copied.
- Strict JSON validation, SI geometry loading, edge-list topology construction, connect-degree cross-check, dynamic `Pl,C[,Pg,Sw]` DOF maps, backward-Euler residuals, dissolution-only frozen transfer branch, frozen upwind branch, local FADBAD++ Jacobian, row/column scaling, fraction-to-boundary, Armijo search, whole-step rollback, SparseLU and BiCGSTAB+ILUT backends.
- Post-convergence irreversible bubble removal with conservative free-to-dissolved transfer and liquid-volume ledger; no reactivation.
- Initial bubbles already below `Sg_off/ng_off` are retired by the same conservative free-to-dissolved transfer, rather than silently discarding their EOS inventory.
- Closed all-liquid connected components receive exactly one standard pressure gauge in place of their redundant volume equation. This retains every supplied control volume while removing the incompressible-pressure nullspace; the selected nodes are exposed as `pressure_gauge` in VTK.
- Legacy VTK node/edge fields, TSV diagnostics, JSON checkpoint/restart, SHA-256, and run-manifest plumbing. Checkpoint v3 retains adaptive `dt`/fast-step state, serializes the nonlinear state as round-trippable extended-precision strings alongside its binary64 compatibility array, and binds restarts to a normalized continuation-config hash plus all three network-input hashes.
- An optional checksum-locked geometry-correction overlay workflow is retained. `configs/geometry_corrections.required.json` enumerates all 516 raw fields but intentionally leaves values/source notes null; none is needed by the degenerate formulation. The applicator refuses incomplete, stale, non-positive, duplicate, or unsupported edits and never mutates raw input.
- Failed configured runs emit `run_manifest.failed.json` with nonzero status, exact command, error, configuration, timestamps, toolchain and hashes of every input that was readable. A later successful run removes a stale failure manifest from the same output directory.
- The overlay workflow is end-to-end tested on temporary synthetic raw files, including exact replacements, immutable source hashes, derived hashes and correction ledger; no test replacement enters the real correction template.
- CMake deterministically hashes the complete FADBAD++ distribution tree, Eigen/nlohmann header trees and the actual `CMakeLists/include/src` project source tree; these hashes, versions, build mode and license/source records are emitted in every success/failure manifest.
- Unsafe line-search, linear-solver, active-set, boundary-flag and blank research-provenance settings fail before mesh loading. Failure manifests retain the actual accepted step/time instead of hard-coding zero.
- Zero half-lengths are eliminated exactly with union-find clusters. Cluster members alias one `Pl/C` DOF pair and contribute to one aggregated `Rl/Rd`; every positive-volume gas member independently retains `Pg/Sw/Rg/Rc`. Zero edges remain in VTK/audit mappings with identically zero numerical flux.
- `V=0,Sw=1` controls are retained as accumulation-free massless junctions. `V<0`, `L<0`, `V=0,Sw<1`, and conflicting-pressure boundary clusters fail fast.
- A stable shared double/FADBAD++ series evaluates `1-exp(-b Sw)` near `Sw_min`, avoiding multi-Pascal cancellation in the Qin closure. The nonlinear state and full-network gas residual accumulation/transfer cancellation use extended precision while the FADBAD++ Jacobian and sparse linear solve remain binary64 and monolithic.
- Newton acceptance requires both `||R_scaled||inf <= residual_tolerance` and `||delta_scaled||inf <= update_tolerance * (1 + ||state_scaled||inf)`. Accepted-step TSV and manifest records include the final total residual, relative update, criterion gates, and maximum scaled `Rl/Rd/Rg/Rc` with source PNM node IDs; a failed gate continues Newton and ultimately triggers timestep retry/cut rather than tolerance relaxation.
- Physical geometry units are explicit and strictly consistent: `input_basis=physical`, `input_length_unit=um`, voxel size `7.5 um`, image dimensions `1200^3`, domain `9000^3 um3`, and SI factors `1e-6/1e-12/1e-18`; area and volume factors must equal the square/cube of the length factor.
- A reproducible generator transfers every real boundary flag to one degree-one `source_type=2` virtual reservoir, copies the old-model radius/geometry rule, appends one positive 7.5 um boundary edge, and writes a machine-readable geometry/overlap/isolation/hash audit. Virtual declared volume is never internal storage.
- Real former boundary pores now have flag zero and ordinary `Rl/Rd[/Rg/Rc]`; only virtual reservoirs hold fixed `Pl/C`. VTK exposes `node_id`, `is_virtual_boundary`, `is_boundary_edge`, SI edge length and edge length in voxels. Per-run `run.log` records Newton convergence, rejected timesteps and elapsed time.

## Verified by actual execution

- Converter regenerated 146,635 storage/control nodes and 224,262 half-edges from 34,504 original pores and 112,131 original throats.
- Supplied validator reported every advertised check PASS; converted storage-volume relative error was `1.377e-16`.
- Release build succeeds with GCC 13.3, Eigen 3.4.0, FADBAD++ 2.1 and nlohmann/json 3.11.2 without project-source warnings.
- FADBAD++ smoke derivative test passes its analytic and central-difference limits.
- Physical Jacobian central-difference audits cover both sides of the dissolution switch and currently report worst relative error `1.30225e-09`, below the `1e-6` requirement.
- Deterministic T1--T8 test executable passes, including SparseLU versus ILUT (`1.08e-19` scaled difference), rollback/retry, conservative initial/post-step disappearance, checkpoint restart, and VTK/checkpoint reread.
- T2 now also advances the production solver on a closed two-control diffusion problem and an isolated liquid control. One pressure gauge per closed gas-free component removes the nullspace while total dissolved moles remain within `1e-12` relative error.
- T4 additionally advances the production `Solver` on a five-node high-Peclet chain in both pressure directions, checking final upwind direction, concentration bounds, equal-and-opposite internal flux and outlet-backflow concentration.
- `dt,dt/2,dt/4` observed orders for total free moles, total dissolved moles, mean Sw, and cumulative outlet dissolved moles were respectively `0.984, 1.053, 2.510, 1.391` in the deterministic chain test.
- The 5-step deterministic CLI/output run has maximum equation-form `epsilon_N=1.158e-16` and `epsilon_V=1.369e-16`.
- Adaptive production-CLI checkpoint/restart is byte-for-byte identical in its diagnostics to an uninterrupted run, including a grown timestep; altered physics/input provenance and forbidden bubble reactivation are rejected. Restart into a new output directory writes a valid TSV header.
- Seventeen unsafe/blank configurations are rejected before mesh loading, including missing unit basis and inconsistent area/volume conversion factors; an induced post-step output failure proves the failure manifest preserves one already accepted step.
- The supplied topology has 711 connected components and 672 without a pressure boundary. The historical base network had 666 isolated controls; after retaining and attaching all former boundary pores one-to-one, the derived network has 631 isolated controls.
- T9--T14 pass: shared two-member and multi-member cluster DOFs, inventory invariance, clustered AD/finite-difference Jacobian, square nonsingular multi-gas Jacobian, analytic two-edge massless junction, degree-one dead end, prohibited degeneracies, and real/virtual boundary DOF-storage-flux-checkpoint behavior.
- Historical pre-unit-correction verification passed 13/13 Release and 11/11 ASan+UBSan; its logs and outputs remain untouched as numerical-regression evidence only.
- Base-network preprocessing classifies all 514 zero half-edges into 219 nontrivial clusters, retains two massless junctions (32123/140815), and produces 146,121 effective clusters plus 223,748 positive edges. The derived network adds 1,509 singleton reservoir clusters and boundary edges without changing those base classifications. No original or converted input was modified.
- The historical nm-interpreted 146,635-control/224,262-edge run accepted 5 steps, but is retained only as a pre-correction numerical-regression record and is not a physical result.
- The corrected derived network contains 148,144 nodes and 225,771 edges: 146,635 real/storage controls, 1,509 virtual reservoirs (724 inlet, 785 outlet), 224,262 original half-edges and 1,509 boundary edges. It retains 514 zero-distance constraints, 219 nontrivial clusters and massless nodes 32123/140815.
- All connected real-virtual clearances are 7.5 um. The old-model independent-offset rule creates 96 virtual/non-connected-real display-envelope overlaps and one virtual/virtual overlap; 75 and one respectively are new relative to the real-boundary layout. No generated edge is zero, no reservoir is shared, and no virtual node enters a zero-distance cluster.
- The 35 isolated former boundary pores (14 inlet, 21 outlet) are retained. Each forms a nonsingular real-pore/virtual-reservoir two-node component and is outside the largest sample component. Thirty-four are all-liquid equilibrium pairs with zero expected flux; inlet pore 18230 has `Sw=0.1874`, may contribute transient boundary flux/dissolution, and remains included in global conservation. None was deleted or connected artificially.
- Strict-warning Release passes 16/16. ASan+UBSan passes 14/14 with only `full_network_five_step_smoke` and its dependent `full_network_output_audit` excluded; full-network unit generation and strict geometry still run under the sanitizer suite.
- The corrected-unit virtual-boundary network accepts five `1e-4 s` steps without retries. Maximum `||R_scaled||inf=2.7244193879434974e-11`, relative update `1.5351302621425878e-10`, equation/adjusted `epsilon_N=1.0033252695258207e-13`, and `epsilon_V=4.135372990154849e-15`. Block maxima are `Rl=2.7244193879434974e-11` at node 544, `Rd=7.784423382026258e-14` at node 3, `Rg=1.3802620372561522e-15` at node 66468, and `Rc=9.004898037355725e-18` at node 52932.

## Unresolved

- The research template intentionally contains null/unconfirmed physical properties and is rejected before mesh loading. The completed full-network run is numerical smoke evidence only, not a scientific gas-dissolution result.
- The corrected run remains a numerical preview: gas/liquid identity, Henry coefficient, mass-transfer coefficient, surface tension/contact angle, Qin effective radius interpretation, pressure boundary values and the `1e-3 mol/m3` background concentration are not experimentally confirmed.
- AMGX remains an optional second-stage backend and is not implemented; CPU/ILUT now has full-network acceptance evidence.

## Completion audit

- [x] clean build -- `output/verification/ctest-physical-release.log`, `output/verification/ctest-physical-asan.log`
- [x] input/unit/virtual geometry audit -- `output/verification/pnm-physical-unit-audit.json`, `output/full_network_physical_um_virtual_boundary/virtual_boundary_geometry_audit.json`
- [x] AD derivative tests -- `test_ad_smoke`, physical Jacobian in `test_bubble_solver`
- [x] T1-T8 -- `test_bubble_solver`, including production closed-component pressure-gauge coverage
- [x] T9-T14 -- cluster/massless/invalid-degenerate/virtual-boundary tests in `test_bubble_solver`
- [x] dt convergence -- `test_bubble_solver`
- [x] corrected full-network >=5 accepted steps -- `output/full_network_physical_um_virtual_boundary/diagnostics.tsv`, `run_manifest.json`
- [x] corrected full-network epsilon_N and epsilon_V thresholds -- maximum `1.0034e-13` / `4.1354e-15`
- [x] VTK/checkpoint re-read -- deterministic test and provenance-bound adaptive CLI restart equivalence
- [x] VTK/TSV/checkpoint-v3 extended-state/manifest/checksum output audits -- `output/cli_virtual_test/`, `output/full_network_physical_um_virtual_boundary/`, `output/verification/ctest-physical-release.log`
- [x] corrected physical-unit and virtual-boundary output audit -- `output/full_network_physical_um_virtual_boundary/`, `output/verification/ctest-physical-release.log`, `ctest-physical-asan.log`
- [ ] scientific parameters confirmed -- explicitly unresolved; smoke parameters are non-scientific
