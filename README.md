# IChem

IChem is a command-line cheminformatics toolkit for analysing protein–ligand and
protein–protein systems. It computes interaction fingerprints (IFP), detects
binding-site cavities (VolSite), fragments molecules, aligns and compares
structures, and related structural analyses.

> **Status:** legacy C++ codebase under incremental cleanup. Expect rough edges
> in the source. Repository hygiene is being addressed first (see below).

## Repository layout

| Path | Contents |
|------|----------|
| `source/`      | C++ source. Modules under `headers/` and `sources/`: `ICMole` (molecular model), `ICCalcs` (interactions, VolSite, fragments, PDB conversion), `ICPars` (parsers/writers), `ICTools` (utilities), `ICRSA`, `ICSVM`. Entry point: `source/main.cpp`. |
| `source/bindings/` | Python bindings (`ichem_ifp_py.cpp`). |
| `external/`    | Vendored third-party headers (e.g. nanoflann). |
| `datas/`       | Reference template data (MOL2) required at runtime. Static; rarely changes. |
| `test/`        | End-to-end test suites driven by `test/test.py` (inputs + reference outputs). |
| `unit_tests/`  | C++ unit tests (CMake). |
| `ci-containers/` | Containerfiles (ubi / ubuntu / debian) used by CI to build release binaries. |
| `.github/workflows/release.yml` | Builds per-distro binaries and publishes GitHub releases on tags. |
| `CHANGELOG.md` | Release history. |
| `User_Guide.pdf` | End-user documentation. |

## Building

Requires CMake ≥ 3.13 and a C++20 compiler (GCC or Clang).

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

This produces the `IChem` executable in `build/`. A `Debug` build is the default
if no build type is given. See `CMakeLists.txt` for options
(`STATIC_ANALYSIS`, `ICHEM_DEBUG`).

Reproducible builds matching the released binaries can be produced inside the
containers in `ci-containers/` (see `release.yml`).

## Running tests

End-to-end tests are orchestrated by a Python harness:

```bash
cd test
python3 test.py
```

Each suite (`ifp`, `volsite`, `frag`, `grim`, `sims`, `realign`, …) runs IChem on
fixed inputs and compares against reference outputs under `*/ref`.

## Releases

Release binaries are **built by CI, not committed**. Pushing a `v*` tag triggers
`.github/workflows/release.yml`, which builds `IChem` for ubi/ubuntu/debian and
attaches the binaries to a GitHub release.
