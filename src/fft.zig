/// Cooley-Tukey radix-2 FFT. Gate parameters below are placeholders,
/// not a connected proof checker. Numerical agreement is checked separately.
const std = @import("std");
const math = std.math;
const Complex = std.math.Complex(f64);
const assert = std.debug.assert;

/// In-place radix-2 DIT FFT.  N must be a power of 2.
/// Modifies `buf` in place; on return buf[k] = Σ x[n]·e^{-2πi·kn/N}.
pub fn fft(buf: []Complex) void {
    const N = buf.len;
    assert(math.isPowerOfTwo(N));

    // Bit-reversal permutation
    var j: usize = 0;
    var i: usize = 1;
    while (i < N) : (i += 1) {
        var bit = N >> 1;
        while (j & bit != 0) : (bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std.mem.swap(Complex, &buf[i], &buf[j]);
    }

    // Butterfly stages
    var len: usize = 2;
    while (len <= N) : (len <<= 1) {
        const ang = -2.0 * math.pi / @as(f64, @floatFromInt(len));
        const wlen = Complex{ .re = @cos(ang), .im = @sin(ang) };
        var pos: usize = 0;
        while (pos < N) : (pos += len) {
            var w = Complex{ .re = 1.0, .im = 0.0 };
            var k: usize = 0;
            while (k < len / 2) : (k += 1) {
                const u = buf[pos + k];
                const v = mulComplex(buf[pos + k + len / 2], w);
                buf[pos + k] = addComplex(u, v);
                buf[pos + k + len / 2] = subComplex(u, v);
                w = mulComplex(w, wlen);
            }
        }
    }
}

/// IFFT: conjugate, FFT, conjugate, divide by N.
pub fn ifft(buf: []Complex) void {
    for (buf) |*x| x.* = conjugate(x.*);
    fft(buf);
    for (buf) |*x| {
        x.* = conjugate(x.*);
        x.re /= @as(f64, @floatFromInt(buf.len));
        x.im /= @as(f64, @floatFromInt(buf.len));
    }
}

/// Real-valued FFT output: returns the N/2+1 non-redundant bins.
/// Caller owns the returned slice (allocator required).
pub fn rfft(allocator: std.mem.Allocator, real_in: []const f64) ![]Complex {
    const N = real_in.len;
    if (N == 0 or !math.isPowerOfTwo(N)) return error.InvalidLength;
    const buf = try allocator.alloc(Complex, N);
    defer allocator.free(buf);
    for (real_in, 0..) |x, k| buf[k] = .{ .re = x, .im = 0.0 };
    fft(buf);
    const out = try allocator.alloc(Complex, N / 2 + 1);
    @memcpy(out, buf[0 .. N / 2 + 1]);
    return out;
}

/// Linear convolution via FFT, zero-padded to the next power of two.
/// `_gate` is a legacy placeholder; no proof manifest is checked.
pub fn convolve(
    comptime _gate: GateConvolution,
    allocator: std.mem.Allocator,
    a: []const f64,
    b: []const f64,
) ![]f64 {
    _ = _gate;
    if (a.len == 0 or b.len == 0) return error.EmptySignal;
    const out_len = a.len + b.len - 1;
    var n: usize = 1;
    while (n < out_len) n <<= 1;

    const fa = try allocator.alloc(Complex, n);
    defer allocator.free(fa);
    const fb = try allocator.alloc(Complex, n);
    defer allocator.free(fb);

    for (fa, 0..) |*x, k| x.* = if (k < a.len) .{ .re = a[k], .im = 0 } else .{ .re = 0, .im = 0 };
    for (fb, 0..) |*x, k| x.* = if (k < b.len) .{ .re = b[k], .im = 0 } else .{ .re = 0, .im = 0 };

    fft(fa);
    fft(fb);
    for (fa, fb) |*xa, xb| xa.* = mulComplex(xa.*, xb);
    ifft(fa);

    const result = try allocator.alloc(f64, out_len);
    for (result, 0..) |*r, k| r.* = fa[k].re;
    return result;
}

/// Placeholder type retained for source compatibility. No proof dependency.
const GateConvolution = struct {};

// ── arithmetic helpers ───────────────────────────────────────────────────────

inline fn addComplex(a: Complex, b: Complex) Complex {
    return .{ .re = a.re + b.re, .im = a.im + b.im };
}
inline fn subComplex(a: Complex, b: Complex) Complex {
    return .{ .re = a.re - b.re, .im = a.im - b.im };
}
inline fn mulComplex(a: Complex, b: Complex) Complex {
    return .{ .re = a.re * b.re - a.im * b.im, .im = a.re * b.im + a.im * b.re };
}
inline fn conjugate(a: Complex) Complex {
    return .{ .re = a.re, .im = -a.im };
}

// ── tests ─────────────────────────────────────────────────────────────────────

test "fft: impulse response is constant" {
    const N = 8;
    var buf: [N]Complex = undefined;
    buf[0] = .{ .re = 1, .im = 0 };
    for (buf[1..]) |*x| x.* = .{ .re = 0, .im = 0 };
    fft(&buf);
    for (buf) |x| {
        try std.testing.expectApproxEqAbs(1.0, x.re, 1e-10);
        try std.testing.expectApproxEqAbs(0.0, x.im, 1e-10);
    }
}

test "fft: single frequency bin" {
    const N = 8;
    var buf: [N]Complex = undefined;
    // x[n] = e^{2πi·k·n/N} for k=1 → X[1]=N, others 0
    for (&buf, 0..) |*x, n| {
        const ang = 2.0 * math.pi * @as(f64, @floatFromInt(n)) / @as(f64, @floatFromInt(N));
        x.* = .{ .re = @cos(ang), .im = @sin(ang) };
    }
    fft(&buf);
    // After FFT the energy is at bin N-1 (conjugate symmetry of DFT of complex exp)
    // The bin at index 1 should have magnitude ≈ N
    const mag1 = math.sqrt(buf[1].re * buf[1].re + buf[1].im * buf[1].im);
    try std.testing.expectApproxEqAbs(@as(f64, N), mag1, 1e-8);
}

test "ifft inverts fft" {
    const N = 16;
    var buf: [N]Complex = undefined;
    for (&buf, 0..) |*x, k|
        x.* = .{ .re = @as(f64, @floatFromInt(k)), .im = 0 };
    const orig = buf;
    fft(&buf);
    ifft(&buf);
    for (buf, 0..) |x, k| {
        try std.testing.expectApproxEqAbs(orig[k].re, x.re, 1e-9);
        try std.testing.expectApproxEqAbs(orig[k].im, x.im, 1e-9);
    }
}

fn rfftAllocationScenario(allocator: std.mem.Allocator) !void {
    const out = try rfft(allocator, &.{ 1, 2, 3, 4 });
    defer allocator.free(out);
}

test "rfft cleans up intermediate storage on every allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, rfftAllocationScenario, .{});
}

test "rfft refuses empty and non-power-of-two inputs" {
    try std.testing.expectError(error.InvalidLength, rfft(std.testing.allocator, &.{}));
    try std.testing.expectError(error.InvalidLength, rfft(std.testing.allocator, &.{ 1, 2, 3 }));
}

test "convolution refuses empty inputs" {
    try std.testing.expectError(error.EmptySignal, convolve(.{}, std.testing.allocator, &.{}, &.{1}));
    try std.testing.expectError(error.EmptySignal, convolve(.{}, std.testing.allocator, &.{1}, &.{}));
}
