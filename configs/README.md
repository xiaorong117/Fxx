# Configuration files

This directory contains reusable solver and linear-backend configuration examples.

- `bubble_solver.research-template.json` is the starting point for a new scientific case. Its null physical fields must be replaced with traceable values.
- `amgx_*.json` files are AMGX linear-solver experiments; AMGX is optional and is not the recommended production backend for the current full network.
- Other top-level `bubble_solver.*.json` files document validated numerical or boundary-condition workflows, but their paths and physical closures must be reviewed before reuse.

Generated batch trees, pilot/recovery configurations, original PNM data and production outputs are intentionally excluded from Git. Create those locally from audited external inputs.
