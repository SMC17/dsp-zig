# dsp-zig

Private OSS library (FFT / Kalman / FIR / windows / Welch PSD). Not a
company, not a primitives product, not a Zig standard-library substitute.

## Proof (re-run this pass)

```sh
zig version                 # 0.16.0
zig build test --summary all
```

**13/13 tests passed.** Proof level: `unit-tested`.

Not measured: gold-reference agreement (FFTW/numpy), soak, SIMD throughput,
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
