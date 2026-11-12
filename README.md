# dsp-zig

Private OSS library (FFT / Kalman / FIR / windows / Welch PSD). Not a
company, not a primitives product, not a Zig standard-library substitute.

## Proof (re-run this pass)

```sh
zig version                 # 0.16.0
zig build test --summary all
```

The original baseline passed 13 tests. The current target adds allocation-failure
cleanup and invalid-input regressions. Run the commands below for the current
count; the task is a numerical fixture check, not formal verification.

## One task

```sh
zig build test demo -Doptimize=ReleaseSafe --summary all
```

The demo checks 31 complex FFT bins against a direct O(N²) DFT, checks the
inverse values and compares two linear convolutions against direct sums.
Expected result: `DSP_FIXTURE_OK dft_bins=31 inverse_bins=31 convolution_cases=2`.
Absolute tolerance: 1e-9. Nonfinite or mismatching results return nonzero.

`rfft` now frees its intermediate buffer when the output allocation fails,
rejects empty/non-power-of-two input, and convolution refuses empty signals.
Allocation-failure probing is part of the real test target.

Create a keepable pack from committed source:

```sh
python3 tools/selfserve.py selfserve.json --repo . --out /tmp/dsp-selfserve
```

The pack contains matching source, the demo binary, logs, result identity and
checksums. Its workflow provides push/PR/manual/monthly acceptance and an
intentional-tag preview gate. It does not change private visibility or prove
hosted execution simply by existing.

Not measured: FFTW/numpy agreement, large-signal soak, SIMD throughput,
numerical stability on ill-conditioned Kalman, or any hardware DSP path.

## What the source comments do not prove

Several files talk about Lean 4, `proof_manifest.json`, and `proof-gate-zig`.
The `Gate*` types in this tree are empty structs, for example:

```zig
const GateConvolution = struct {};
```

There is no Lean, no `proof_manifest.json`, and no `proof-gate-zig`
dependency in `build.zig.zon`. A refuted formal claim would **not** stop
this package from compiling. Do not call this `audited` or formally
verified.

## Visibility

Keep private. Do not public-flip from this file.
