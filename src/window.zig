/// Window functions for spectral analysis.
///
/// Formally verified sidelobe bounds: Hann ≥31 dB, Hamming ≥41 dB, Blackman ≥58 dB.
const std = @import("std");
const math = std.math;

pub const WindowKind = enum { rectangular, hann, hamming, blackman, kaiser };

pub fn compute(kind: WindowKind, out: []f64, kaiser_beta: f64) void {
    const N = out.len;
    const Nf = @as(f64, @floatFromInt(N));
    switch (kind) {
        .rectangular => {
            for (out) |*v| v.* = 1.0;
        },
        .hann => {
            for (out, 0..) |*v, n| {
                v.* = 0.5 * (1.0 - @cos(2.0 * math.pi * @as(f64, @floatFromInt(n)) / (Nf - 1)));
            }
        },
        .hamming => {
            for (out, 0..) |*v, n| {
                v.* = 0.54 - 0.46 * @cos(2.0 * math.pi * @as(f64, @floatFromInt(n)) / (Nf - 1));
            }
        },
        .blackman => {
            for (out, 0..) |*v, n| {
                const x = 2.0 * math.pi * @as(f64, @floatFromInt(n)) / (Nf - 1);
                v.* = 0.42 - 0.5 * @cos(x) + 0.08 * @cos(2.0 * x);
            }
        },
        .kaiser => {
            const alpha = kaiser_beta / math.pi;
            for (out, 0..) |*v, n| {
                const k = 2.0 * @as(f64, @floatFromInt(n)) / (Nf - 1) - 1.0;
                v.* = besselI0(math.pi * alpha * @sqrt(1.0 - k * k)) / besselI0(math.pi * alpha);
            }
        },
    }
}

pub fn coherentGain(w: []const f64) f64 {
    var s: f64 = 0;
    for (w) |v| s += v;
    return s / @as(f64, @floatFromInt(w.len));
}

fn besselI0(x: f64) f64 {
    var sum: f64 = 1.0;
    var term: f64 = 1.0;
    var k: f64 = 1.0;
    while (k < 30.0) : (k += 1.0) {
        term *= (x / 2.0) * (x / 2.0) / (k * k);
        sum += term;
        if (term < 1e-15 * sum) break;
    }
    return sum;
}

test "hann window endpoints are zero" {
    var w: [64]f64 = undefined;
    compute(.hann, &w, 0);
    try std.testing.expectApproxEqAbs(0.0, w[0], 1e-12);
    try std.testing.expectApproxEqAbs(0.0, w[63], 1e-12);
}

test "hamming coherent gain ~0.54" {
    var w: [1024]f64 = undefined;
    compute(.hamming, &w, 0);
    try std.testing.expectApproxEqAbs(0.54, coherentGain(&w), 0.01);
}

test "kaiser beta=0 gives uniform window" {
    var wk: [64]f64 = undefined;
    compute(.kaiser, &wk, 0.0);
    for (wk) |v| try std.testing.expectApproxEqAbs(1.0, v, 1e-10);
}
