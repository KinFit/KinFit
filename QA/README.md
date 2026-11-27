# KinFit QA test suite

Automated tests for KinFit, run with CTest. ROOT is the only dependency.

## Running

```bash
cmake -S . -B build                  # KINFIT_BUILD_TESTS is ON by default
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Useful variants:

```bash
ctest --test-dir build -L fitter               # only one label
ctest --test-dir build -R test_kinfitter.pulls # name regex
./build/QA/test_kinfitter --list               # cases of one executable
./build/QA/test_kinfitter pulls_4C             # run a single case directly
```

Each case is its own CTest entry, named `<executable>.<case>`. Labels:
`core`, `fitter`, `vertex`, `analysis`, `io`, `example`.

## Contents

| Executable | Covers |
| --- | --- |
| `test_particle` | `KinFitParticle`: construction, (R, Z) convention and sign of R, covariance, `TClonesArray` I/O round trip |
| `test_vertexfinder` | `KFitVertexFinder` and `KFitVertexFinderConstB`: exact vertices, loopers, ±Bz, B = 0, `fit()` after setters, singular input |
| `test_kinfitter` | `KinFitter`, for each constraint: `exact_*` (truth input), `pulls_*` (smeared input: constraint fulfilled, pulls N(0,1), flat probability), `jacobian_*` (analytic derivatives vs finite differences); plus edge cases |
| `test_decaycandfinder` | `KFitDecayCandFinder`: momentum, direction, R/Z, covariance, use in a 3C fit |
| `test_coordinateconversion` | `CoordinateConversion`: both directions and error propagation |
| `test_analyzer` | `KFitDecayBuilder` combinatorics and the `KFitAnalyzer` file-to-file workflow |
| `example.*` | the self-contained examples as smoke tests (with `-DKINFIT_BUILD_EXAMPLES=ON`) |

Helpers:

- [KinFitTest.h](KinFitTest.h): a minimal test harness (`KF_TEST`,
  `KF_CHECK`, `KF_CHECK_NEAR`, …).
- [ToyTracks.h](ToyTracks.h): a toy generator for straight and helix tracks
  with known truth, Gaussian smearing, Λ → pπ⁻ events, and statistics.

Statistical tests use fixed random seeds and 1000 events, with tolerances
of about 4σ, so they are deterministic and stable.

## Known bugs: disabled tests

A test for a known bug is written with the **correct** expectation and listed
under `DISABLED` in [CMakeLists.txt](CMakeLists.txt), with a one-line reason.
CTest reports it as "Not Run (Disabled)", so the suite stays green and the bug
stays visible. To fix such a bug:

1. Remove the case from `DISABLED`.
2. Run it and see it fail: `./build/QA/<exe> <case>`.
3. Fix the code until it passes, then run the whole suite.

Currently no case is disabled: all known bugs covered by tests have been
fixed (see [TODO.md](../TODO.md), section 3). Use the mechanism for new
findings.

## AddressSanitizer

```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DKINFIT_ENABLE_ASAN=ON
cmake --build build-asan -j
ASAN_OPTIONS=detect_leaks=0 ctest --test-dir build-asan --output-on-failure
```

Notes:

- On macOS, ROOT's TBB library crashes at program exit under ASan. The test
  executables therefore end with `_Exit()` in ASan builds.
