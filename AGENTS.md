# Codex / automated-agent instructions

This file governs automated work in `/workspace/zz/bubble_dissolution_solver_amgx`.
Read it completely before inspecting, changing, building, or running the project.
Also read the relevant sections of:

- `docs/CLIENT_DELIVERY_GUIDE_ZH.md`
- `docs/PHYSICS_IMPLEMENTATION_ZH.md`
- `docs/USER_GUIDE_ZH.md`

## 1. Primary rule

Preserve scientific provenance. Do not trade correctness, conservation,
reproducibility, or old evidence for apparent speed or convenience.

## 2. Protected inputs and evidence

Never modify, rename, delete, overwrite, or regenerate in place:

- `/workspace/zz/bubble_dissolution_solver`
- `/workspace/zz/bubble_dissolution_solver/output/full_network_air_water_relative_pressure_run01`
- `/workspace/zz/zz_projects/35级配`
- completed or failed directories under `output/verification/`
- existing directories under `output/production/`
- an active run's config, binary, checkpoint, TSV, VTP, manifest, or audit files

Raw PNM files are immutable. Put every conversion or correction in a new derived
directory and record before/after SHA-256. Do not fabricate a missing PNM variant.

## 3. Before any action

For read-only questions, inspect and report; do not mutate. Before a requested
change or run:

1. inspect `git status`/the actual diff without reverting user changes;
2. list relevant tmux sessions and exact solver/scheduler processes;
3. resolve the executable through `/proc/PID/exe` and hash it;
4. resolve and hash the actual config and PNM inputs;
5. read the latest diagnostics, checkpoint, manifest/failure manifest and batch
   progress;
6. state whether the proposed action affects a running process;
7. use a new build, verification, config and output path.

Do not kill a process merely because its name resembles the target. Resolve PID,
argv, output directory and run lock first. Never stop unrelated user processes.

## 4. Mathematical and physical invariants

There is one AD residual/Jacobian and one Newton workflow. Do not create a second
physics implementation for another linear backend. Preserve:

- monolithic `Rl/Rd/Rg/Rc` equations;
- backward-Euler discretization;
- relative pressure state and absolute pressure in EOS/Henry relations;
- clusters, zero-distance constraints and massless junctions;
- active-set irreversibility unless explicitly redesigning and validating it;
- both Newton acceptance criteria;
- state bounds and nonnegative free-gas moles;
- molar/component and liquid-volume conservation gates;
- identical interfacial transfer in dissolved and free-gas balances.

Never silently fall back from AMGX/PARDISO to Eigen. Never treat an AMGX
`NOT_CONVERGED` status as success. Never loosen Newton, linear, state or
conservation tolerances merely to pass a run.

`mass_transfer_multiplier` must multiply the same conservative interfacial flux
on both sides. Scientific cases use `1.0`; other values are preview-only.

## 5. Unconfirmed physics

Do not describe configured values as experimentally confirmed unless an
authoritative source is added. At minimum keep Henry coefficient, kL,
interfacial-area closure, Qin effective radius/b, contact angle, relative
conductance, diffusion-area/tortuosity, immobile-bubble assumption and selected
Pe/Ca flow range marked assumed or unconfirmed.

Changing a physical property, boundary type/value, inlet concentration,
mass-transfer multiplier, PNM geometry or active-set physics creates a new case
that starts at `t=0`. Do not restart an incompatible checkpoint.

## 6. Builds and tests

Use out-of-source build directories. Do not replace a validated binary while a
production process uses it. Capture the configure command, build log, source
aggregate hash and executable hash.

Run tests proportionate to the change:

- documentation/UI: relevant Python tests plus syntax and API checks;
- solver/config: Release warnings and affected unit/regression tests;
- checkpoint/output: continuous-vs-restart and schema audits;
- CPU memory safety: ASan/UBSan where applicable;
- AMGX: GPU tests only when CUDA driver/device are genuinely available; record a
  precise skip reason otherwise;
- production: representative small network, then short full-network pilot, then
  long run.

Preserve all failed logs and intermediate evidence. Do not delete a failure to
make a later report look clean.

## 7. Time stepping and bubble retirement

The adaptive controller extends the existing retry path; do not add a second,
conflicting dt controller. A step is accepted only after both Newton criteria,
state bounds, nonnegative gas and conservation checks pass.

Concentration change caps must use an explicit physical reference concentration,
not the minimum local `C/|dC/dt|` near zero. A negligible bubble must not hold the
global dt indefinitely.

Conservative batch retirement is allowed only when enabled in config and both
hard fraction limits pass. Remaining positive gas is transferred to dissolved
inventory, retired gas volume enters the active-set liquid ledger, and the
unchanged conservation gates must pass. Batch event time is estimated within the
accepted step; do not call it an exactly re-solved event time.

## 8. Checkpoint/restart

Before restart verify checkpoint JSON, node count, time, step, `next_dt_s`,
continuation hash, PNM hashes, pressure representation, mass-transfer multiplier,
adaptive counters, cumulative ledgers and irreversible active-set metadata.

Do not edit a checkpoint to bypass compatibility. Prevent duplicate TSV times and
preserve cumulative/event continuity. Archive STOP, manifest and transition
evidence instead of erasing them.

## 9. Batch execution

The authoritative current inventory contains 35 PNM-0.6 and 35 PNM-0.8 cases.
Use:

```text
configs/batch35_flow_pe555_pardiso_v3/batch_cases.json
output/production/batch35_flow_pe555_pardiso/batch_progress.json
```

Only one scheduler may own an index/output root. `--jobs 24` means one scheduler
with 24 slots, not three independent schedulers reading the same 70 cases. Every
case uses `run.lock`. If any case fails, preserve it and stop launching more until
the cause is classified.

Do not individually STOP a child while the scheduler is active: it will refill
the slot. A whole-batch concurrency handoff must block queued locks, safely STOP
active cases after accepted steps, verify all checkpoints, archive transition
evidence, release locks, and then start one replacement scheduler. Never signal a
shared tmux process group without proving which descendants will receive it.

## 10. Dashboard and external access

The dashboard may only scan descendants of:

```text
/workspace/zz/bubble_dissolution_solver_amgx/output/production
```

Keep it bound to `127.0.0.1`. Use SSH forwarding or an authenticated HTTPS
gateway. Do not expose port 8765 unauthenticated on `0.0.0.0`. Do not reveal
tokens, passwords, private keys, credential environment variables or unrelated
system information.

Batch children are read-only in the control UI so operations cannot bypass the
scheduler. Dashboard warnings never automatically stop or modify a solver.

## 11. Current validated production identity

The validated batch-retirement production binary is:

```text
/workspace/zz/build-bubble-batch-retirement-release/bubble_solver
SHA-256 6ea4e9c1a3e00f3919e2415edeaa49e527bbd6f7c9042bbb98d7ee44695ebd04
```

It uses PARDISO and conservative batch retirement for the current 70-case run.
Treat PIDs, progress counts and tmux lifetime as volatile: inspect them live; do
not rely on a PID copied from documentation.

## 12. Final report requirements

Report facts backed by final source, logs and machine-readable outputs. Include:

- files changed/added;
- source/config/binary/input hashes;
- backend actually compiled and selected;
- exact commands and output locations;
- test passes, failures and environment-based skips;
- numerical differences and conservation maxima where relevant;
- run/tmux/PID/progress/restart state;
- termination reason and whether the scientific target was reached;
- all still-unconfirmed physical closures.

Never claim completion while a required run is only queued, running, safely
stopped, censored, or failed.
