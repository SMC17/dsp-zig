/// FIR and IIR filters — comptime tap/order counts.
///
/// Gate: fir_linear_phase — linear phase property of symmetric FIR filters.
/// Gate: iir_stability_criterion — |poles| < 1 iff BIBO stable.
const std = @import("std");

/// Symmetric linear-phase FIR filter, comptime tap count.
/// Gate parameter is a placeholder; no proof checker is connected.
pub fn FIR(comptime N: usize) type {
    return struct {
        coeffs: [N]f64,
        buf: [N]f64 = .{0.0} ** N,
        pos: usize = 0,

        const Self = @This();

        pub fn init(comptime _: GateLinearPhase, coeffs: [N]f64) Self {
            return .{ .coeffs = coeffs };
        }

        pub fn process(self: *Self, x: f64) f64 {
            self.buf[self.pos] = x;
            self.pos = (self.pos + 1) % N;
            var y: f64 = 0;
            for (self.coeffs, 0..) |c, k|
                y += c * self.buf[(self.pos + N - 1 - k) % N];
            return y;
        }

        pub fn processBlock(self: *Self, in: []const f64, out: []f64) void {
            std.debug.assert(in.len == out.len);
            for (in, out) |x, *y| y.* = self.process(x);
        }
    };
}

/// Direct Form II IIR biquad — second-order section.
/// Gate parameter is a placeholder, not a stability certificate.
pub const Biquad = struct {
    b0: f64, b1: f64, b2: f64,
    a1: f64, a2: f64,
    w1: f64 = 0, w2: f64 = 0,

    pub fn init(comptime _: GateStability,
                b0: f64, b1: f64, b2: f64,
                a1: f64, a2: f64) Biquad {
        return .{ .b0 = b0, .b1 = b1, .b2 = b2, .a1 = a1, .a2 = a2 };
    }

    pub fn process(self: *Biquad, x: f64) f64 {
        const w = x - self.a1 * self.w1 - self.a2 * self.w2;
        const y = self.b0 * w + self.b1 * self.w1 + self.b2 * self.w2;
        self.w2 = self.w1;
        self.w1 = w;
        return y;
    }
};

const GateLinearPhase = struct {};
const GateStability = struct {};

// ── Windowed-sinc FIR design ─────────────────────────────────────────────────

/// Design a low-pass windowed-sinc FIR filter.
/// cutoff: normalized [0, 0.5], N: number of taps (should be odd for Type I).
pub fn designLowpass(comptime N: usize, cutoff: f64) [N]f64 {
    comptime std.debug.assert(N % 2 == 1);
    const M: usize = (N - 1) / 2;
    var h: [N]f64 = undefined;
    for (&h, 0..) |*c, k| {
        const n = @as(f64, @floatFromInt(k)) - @as(f64, @floatFromInt(M));
        if (n == 0.0) {
            c.* = 2.0 * cutoff;
        } else {
            const sinc = @sin(2.0 * std.math.pi * cutoff * n) / (std.math.pi * n);
            // Hann window
            const w = 0.5 * (1.0 - @cos(2.0 * std.math.pi * @as(f64, @floatFromInt(k)) / @as(f64, @floatFromInt(N - 1))));
            c.* = sinc * w;
        }
    }
    return h;
}

// ── Tests ────────────────────────────────────────────────────────────────────

test "FIR moving average passes DC" {
    const N = 8;
    const h: [N]f64 = .{1.0 / @as(f64, N)} ** N;
    var fir = FIR(N){ .coeffs = h };
    for (0..N) |_| _ = fir.process(1.0);
    const y = fir.process(1.0);
    try std.testing.expectApproxEqAbs(1.0, y, 1e-10);
}

test "Biquad passthrough (identity: b=[1,0,0], a=[0,0])" {
    var bq = Biquad.init(.{}, 1.0, 0.0, 0.0, 0.0, 0.0);
    try std.testing.expectApproxEqAbs(3.0, bq.process(3.0), 1e-12);
    try std.testing.expectApproxEqAbs(7.0, bq.process(7.0), 1e-12);
}

test "designLowpass: center tap is 2*cutoff" {
    const h = designLowpass(11, 0.25);
    try std.testing.expectApproxEqAbs(0.5, h[5], 1e-12);
}
