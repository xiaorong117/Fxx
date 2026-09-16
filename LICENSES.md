# Dependency and provenance record

| Component | Version/source | License/provenance | Audited SHA-256 |
|---|---|---|---|
| FADBAD++ | 2.1, user-provided `/opt/FADBAD++-2.1/FADBAD++` | Ole Stauning copyright/dual-license notice in `COPYRIGHT` and `README`; no headers vendored | full distribution tree is hashed by CMake/run manifest; `badiff.h`: `e66d01fecfd5c47db7d2507021b89016b59e6d0c0c1da1a10e58f2b518ca0040`; `fadiff.h`: `ccf4d7b0f98eafdc49a07960d11adbc41e17e7baa5d24d7473cd44f95c97deb0` |
| Eigen | Ubuntu `libeigen3-dev` 3.4.0-4build0.1, http://eigen.tuxfamily.org | MPL-2.0 with component notices in Debian copyright metadata | deterministic `Eigen/` header-tree SHA-256 emitted by CMake/run manifest |
| nlohmann/json | Conda headers 3.11.2, https://github.com/nlohmann/json | MIT/SPDX headers | deterministic `nlohmann/` header-tree SHA-256 emitted by CMake/run manifest |
| Intel oneMKL PARDISO | pip/Conda runtime 2025.2.0, `/opt/conda/envs/Mesh/lib/libmkl_rt.so.2` | Intel Simplified Software License reported by the installed `mkl` package; runtime is dynamically loaded and not vendored | `3cc6ea8cf75c994cf9cc084f17326e4a3610292b64fc2c0d7b23e2e34e5c87c7` |
| pnextract input | Imperial College London `pnflow` project's `pnextract` | scientific data-source attribution supplied with case | input hashes recorded in run manifest |

The new project source has not been assigned a redistribution license by the data owner. FADBAD++ usage and any redistribution must comply with its original notice; this project does not relicense that dependency.
