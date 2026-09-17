# Bubble dissolution pore-network solver

中文文档：

- [甲方交付与操作手册](docs/CLIENT_DELIVERY_GUIDE_ZH.md)
- [全耦合控制方程实现更新版（DOCX）](docs/静止气泡溶解_液相对流扩散孔隙网络模型_全耦合控制方程_实现更新版_20260912.docx)
- [操作手册](docs/USER_GUIDE_ZH.md)
- [物理与数值实现](docs/PHYSICS_IMPLEMENTATION_ZH.md)
- [70个有效算例的气泡耐久性统计分析报告](docs/70样品气泡耐久性差异统计分析报告_20260917.md)
- [70个算例结果报告：宏观输入—孔隙结构—气泡耐久性](docs/70算例结果报告_宏观输入-孔隙结构-气泡耐久性_20260917.md)
- [宏观参数—孔隙结构—气泡耐久性关联分析（级配/密实度/围压）](docs/宏观参数-孔隙结构-气泡耐久性关联分析_20260917.md)
- [69个有效算例的气泡耐久性统计分析报告（上一版）](docs/69样品气泡耐久性差异统计分析报告_20260916.md)
- [GitHub仓库范围说明](docs/GITHUB_REPOSITORY_SCOPE_ZH.md)
- [Codex/自动代理约束](AGENTS.md)

This independent C++17 CPU solver implements the monolithic backward-Euler system for one-way dissolution of immobile trapped gas bubbles. Internal liquid-filled controls retain `Pl,C`; active bubble controls retain `Pl,C,Pg,Sw`. There is no gas-network flux, invasion, migration, snap-off, coalescence, or bubble reactivation.

All calculations use SI internally and absolute pressure in EOS and Henry relations. The supplied pnextract coordinates, radii and lengths are physical micrometres, areas are `um2`, and volumes are `um3`; the explicit conversion is therefore `1e-6/1e-12/1e-18`. The evidence is machine-readable in `output/verification/pnm-physical-unit-audit.json`: all coordinate axes span `3.75--8996.3 um`, the image metadata gives `1200*7.5 um=9000 um`, and the pnextract writer multiplies coordinates/radii/lengths by voxel size and volumes by its cube. The pore `radius` is provisionally used both as hydraulic geometry input and as the Qin effective radius; those meanings are not physically identical and require later calibration. Smoke configuration values are software placeholders and must not be interpreted scientifically.

Disconnected all-liquid components are retained. Because incompressible pressure in such a component is defined only up to a constant, the solver replaces one redundant liquid-volume row per component with a pressure gauge anchored to the previous accepted pressure. The gauge selection is written as the VTK point field `pressure_gauge`; it does not alter volume, concentration, geometry or edge fluxes.

Zero half-lengths are represented as exact zero-distance constraints, not as guessed positive lengths. A union-find quotient gives every cluster one shared `Pl,C` pair; each positive-volume gas-bearing member keeps its own `Pg,Sw`, EOS inventory, interface area and Qin radius. The cluster receives one aggregated `Rl,Rd` pair, while every gas member receives its own `Rg,Rc`. Zero-volume nodes with `Sw=1` remain as algebraic, accumulation-free massless junctions. Negative volume/length and `V=0,Sw<1` remain hard errors.

## Physical units and one-to-one virtual boundaries

The derived network is written without modifying the five original/converted PNM inputs:

- `PNM_new_pore_virtual.txt`: the original 146,635 controls followed by 1,509 `source_type=2` virtual reservoirs. `source_id` points to the corresponding real pore. The real pore flag is cleared and transferred to its virtual pore.
- `PNM_new_throat_virtual.txt`: the original 224,262 half-edges followed by 1,509 boundary edges. Generated edges use `source_old_throat_id=0`, side 1/2 for inlet/outlet, and a positive net length of `extension_voxels*voxel_size_um=7.5 um`.
- `PNM_new_pore_connect_virtual.txt`: derived coordination table; every virtual node has coordination one.
- `PNM_new_virtual_boundary_geometry_audit.json`: counts, clearances, overlap classification, isolated boundary records and all base/derived hashes.

The rule reproduces the old model in `/workspace/Mesh/src/Mesh/Muilti_PNM_load_np_array.py`: copy the real pore radius and place the virtual centre at `2*r_real+L_boundary` outside the sample. In the richer current format, radius/area/shape factor and the declared pore volume are copied from the real pore; virtual volume is explicitly excluded from every internal-storage and conservation sum. Boundary-edge radius/area/shape factor are copied and its source-volume metadata is zero.

The 724 inlet flags occupy the low-x side and the 785 outlet flags occupy the high-x side. The derived virtual centres extend along negative/positive x respectively. Each connected real-virtual clearance is `7.5 um`. Spherical display envelopes produce 96 non-connected virtual-real overlaps and one virtual-virtual overlap; these are reported but not “fixed” by changing edge length because all generated boundary edges remain positive, one-to-one and outside zero-distance clusters.

In ParaView, use point filter `is_virtual_boundary == 1` for reservoirs and cell filter `is_boundary_edge == 1` for boundary connections. Point fields include `node_id`, `source_type`, `boundary_flag`, `active_gas`, `Pl`, `Pg`, `Sw`, and `C`; edge fields include `edge_length`, `edge_length_voxels`, and `source_old_throat_id`.

## Build and tests

From `/workspace/zz`:

```bash
python3 bubble_dissolution_solver/scripts/audit_pnm_physical_units.py \
  --root . --pnextract-source /workspace/pnflow/src/pnm/pnextract/blockNet_write_cnm.cpp \
  --output bubble_dissolution_solver/output/verification/pnm-physical-unit-audit.json
python3 bubble_dissolution_solver/scripts/generate_virtual_boundary_network.py \
  --base-pore PNM_new_pore.txt --base-throat PNM_new_throat.txt \
  --base-connect PNM_new_pore_connect.txt \
  --output-pore PNM_new_pore_virtual.txt \
  --output-throat PNM_new_throat_virtual.txt \
  --output-connect PNM_new_pore_connect_virtual.txt \
  --audit-output PNM_new_virtual_boundary_geometry_audit.json \
  --axis x --extension-voxels 1 --voxel-size-um 7.5 \
  --old-model-pore /workspace/SinglePhase/Project/kong/fina_test/reversed_full_pore.txt \
  --old-model-throat /workspace/SinglePhase/Project/kong/fina_test/reversed_full_throat.txt
cmake -S bubble_dissolution_solver -B build-bubble-physical-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DFADBAD_ROOT=/opt/FADBAD++-2.1/FADBAD++ -DNLOHMANN_JSON_ROOT=/opt/conda
cmake --build build-bubble-physical-release -j
ctest --test-dir build-bubble-physical-release --output-on-failure
```

The optional oneMKL PARDISO backend is enabled independently and never falls
back silently:

```bash
cmake -S bubble_dissolution_solver_amgx -B build-bubble-pardiso-release \
  -DCMAKE_BUILD_TYPE=Release -DBUBBLE_ENABLE_AMGX=OFF \
  -DBUBBLE_ENABLE_PARDISO=ON -DPARDISO_MKL_ROOT=/opt/conda/envs/Mesh
```

Select it with `"linear.backend": "pardiso"` and
`"linear.pardiso_threads": 16`. The implementation retains the same AD
Jacobian and Newton equations, performs symbolic analysis once for a fixed
canonical 2x2 sparsity graph, and repeats only numerical factorization and
solve as coefficients change.

The degenerate-geometry audit classifies all supplied zeros and passes without modifying the raw network. The base topology has 514 zero-distance constraints, 219 nontrivial clusters, 146,121 effective liquid clusters and 223,748 positive production edges. Adding the virtual leaves produces 147,630 clusters and 225,257 positive edges while retaining massless junctions 32123/140815.

An optional, checksum-locked correction overlay remains available at `configs/geometry_corrections.required.json` if authoritative replacement geometry is supplied in the future. Its 516 values and source notes deliberately remain null and are not required by the zero-distance/massless formulation. To apply a future authoritative overlay without touching the originals:

```bash
python3 bubble_dissolution_solver/scripts/geometry_correction_overlay.py \
  --root . --overlay bubble_dissolution_solver/configs/geometry_corrections.required.json \
  --output-dir bubble_dissolution_solver/output/corrected_raw
```

Then run the converter and validators inside that derived directory before changing the smoke configuration paths. The overlay refuses missing, duplicate, stale-checksum, non-positive, or unproven corrections.

The directly reproducible corrected-unit/virtual-boundary preview command is:

```bash
./build-bubble-physical-release/bubble_solver \
  --config bubble_dissolution_solver/configs/bubble_solver.smoke.json
```

The new output is `output/full_network_physical_um_virtual_boundary/` and contains six VTK files, `run.log`, per-step `diagnostics.tsv`, a reread checkpoint, regularization ledger, both geometry audits, and `run_manifest.json`. The full derived network advances five accepted steps with `residual_tolerance=1e-8` and `update_tolerance=1e-9`. Inputs with `Sw<=Sw_min` are moved to a small interior value (`100*Sw_min`, bounded by the admissible interval) and the added liquid/moles are written to `initial_regularization.tsv`. The numerical-preview concentration is `1e-3 mol/m3`, deliberately far below `Hcp*Pg`; it is a numerical boundary-state choice, not a confirmed experimental value.

The old `output/full_network_smoke/` and `output/full_network_preview_20260905/` trees are preserved. Their earlier `nm` interpretation is a historical numerical-regression record and must not be cited as a physical result.

Checkpoint v3 includes the next adaptive timestep, fast-convergence counter, a continuation-config signature and pore/throat/connect hashes. A restart rejects changed physics/numerics or changed input contents; extending `t_end_s` and relocating output are allowed:

```bash
./build-bubble-physical-release/bubble_solver \
  --config bubble_dissolution_solver/configs/case-extended.json \
  --restart bubble_dissolution_solver/output/case/checkpoint.json
```

A research configuration with null, blank or unconfirmed properties fails before reading the network.
# Read-only production case dashboard

The fixed case-information panel is served by `scripts/case_dashboard.py`. It
reads only case directories recursively confined below the resolved production
root (`output/production` by default), and combines the actual config, manifest,
checkpoint, diagnostics, audits, run log and matching solver process. Values in
the parameter table carry `configured`, `effective`, `source` and a
`classification` (`confirmed`, `assumed`, `preview-only`, or `unconfirmed`);
conflicts are reported as warnings and are never silently reconciled.
The lightweight case-list API also reads bounded batch scheduler metadata, so
queued cases are visible before their output directories exist. Large history,
PNM and segment data are loaded only for the selected case. The list can be
filtered by PNM-0.6/PNM-0.8, run status, or case name and summarizes running,
queued, completed and failed counts.

Run a one-shot machine-readable snapshot:

```bash
python3 scripts/case_dashboard.py --snapshot
```

Serve the browser panel locally (read-only, no credentials or environment
variables are read):

```bash
python3 scripts/case_dashboard.py --host 127.0.0.1 --port 8765
```

Each monitored production tmux session may expose a dedicated window named
`live` that runs only `tail -F` on that case's `run.log`. The dashboard's
“Live tmux solver console” captures only this approved window, refreshes every
two seconds, and falls back to a bounded `run.log` tail if the window is not
available. It never captures arbitrary shell panes, command history,
environment variables, control tokens, or files outside `output/production`.

## Conservative batched microbubble retirement

The optional `active_set.batch_microbubble_retirement` fast path prevents a
negligible threshold-crossing bubble from forcing a separate global event
timestep. It is disabled by default. A candidate is batched only when both its
individual free-gas fraction and the aggregate candidate fraction for the
accepted step remain below configured hard limits. Negative gas predictions,
significant bubbles, or an aggregate-limit overflow retain exact event
localization. Every retired bubble's remaining positive moles are transferred
to its cluster dissolved inventory and its gas volume enters the existing
active-set liquid-volume ledger before the unchanged conservation gates run.
Per-bubble evidence is written to `gas_disappearance_events.tsv`.

Restricted control mode is opt-in. It requires a 0600 token file and a fixed
solver binary; it permits safe STOP, checkpoint restart, and versioned changes
to `dt_max_s`/checkpoint/VTK intervals only. See
[`docs/USER_GUIDE_ZH.md`](docs/USER_GUIDE_ZH.md) for setup and safety rules.

Open `http://127.0.0.1:8765/`. The `/api/cases`, `/api/case` and nested console
endpoints are restricted to the same production root and reject traversal or
all other paths; they expose case identity,
network/geometry, physical and numerical parameters, source annotations,
runtime values, checkpoint provenance, process information and read-only
warnings. The dashboard does not stop, restart, or modify a solver.
