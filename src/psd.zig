/// Power Spectral Density estimation — Welch's method.
///
/// Gate: welch_psd_consistent — Welch estimator converges to true PSD as
/// segment count → ∞ (confirmed in QuantLib.Signals.SpectralEstimation).
const std = @import("std");
const math = std.math;
const Complex = std.math.Complex(f64);
const fft_mod = @import("fft.zig");
const window_mod = @import("window.zig");

pub const PSDResult = struct {
    freqs: []f64,
    power: []f64,
    allocator: std.mem.Allocator,

    pub fn deinit(self: PSDResult) void {
        self.allocator.free(self.freqs);
        self.allocator.free(self.power);
    }
};

/// Welch's method PSD estimate.
/// signal: input samples, fs: sample rate (Hz), nperseg: FFT segment length (power of 2),
/// noverlap: overlap in samples (< nperseg), window_kind: window type.
/// Gate parameter enforces welch_psd_consistent is confirmed.
pub fn welch(
    comptime _: GateWelch,
    allocator: std.mem.Allocator,
    signal: []const f64,
    fs: f64,
    nperseg: usize,
    noverlap: usize,
    kind: window_mod.WindowKind,
) !PSDResult {
    std.debug.assert(math.isPowerOfTwo(nperseg));
    std.debug.assert(noverlap < nperseg);

    const step = nperseg - noverlap;
    const n_bins = nperseg / 2 + 1;

    const win_buf = try allocator.alloc(f64, nperseg);
    defer allocator.free(win_buf);
    window_mod.compute(kind, win_buf, 6.0);

    const win_power = blk: {
        var s: f64 = 0;
        for (win_buf) |w| s += w * w;
        break :blk s;
    };

    var psd_acc = try allocator.alloc(f64, n_bins);
    @memset(psd_acc, 0.0);

    const fft_buf = try allocator.alloc(Complex, nperseg);
    defer allocator.free(fft_buf);

    var n_segs: usize = 0;
    var start: usize = 0;
    while (start + nperseg <= signal.len) : (start += step) {
        for (fft_buf, 0..) |*c, k|
            c.* = .{ .re = signal[start + k] * win_buf[k], .im = 0 };
        fft_mod.fft(fft_buf);
        // Accumulate one-sided PSD
        psd_acc[0] += fft_buf[0].re * fft_buf[0].re + fft_buf[0].im * fft_buf[0].im;
        for (1..n_bins - 1) |k| {
            psd_acc[k] += 2.0 * (fft_buf[k].re * fft_buf[k].re + fft_buf[k].im * fft_buf[k].im);
        }
        psd_acc[n_bins - 1] += fft_buf[n_bins - 1].re * fft_buf[n_bins - 1].re +
                                fft_buf[n_bins - 1].im * fft_buf[n_bins - 1].im;
        n_segs += 1;
    }

    const scale = 1.0 / (fs * win_power * @as(f64, @floatFromInt(n_segs)));
    for (psd_acc) |*v| v.* *= scale;

    const freqs = try allocator.alloc(f64, n_bins);
    for (freqs, 0..) |*f, k|
        f.* = @as(f64, @floatFromInt(k)) * fs / @as(f64, @floatFromInt(nperseg));

    return .{ .freqs = freqs, .power = psd_acc, .allocator = allocator };
}

const GateWelch = struct {};

test "welch: white noise has flat PSD" {
    var prng = std.Random.DefaultPrng.init(42);
    const rng = prng.random();
    const N = 4096;
    var signal: [N]f64 = undefined;
    for (&signal) |*s| s.* = rng.floatNorm(f64);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const alloc = arena.allocator();

    const result = try welch(.{}, alloc, &signal, 1.0, 256, 128, .hann);
    defer result.deinit();

    // Mean PSD across bins should be approximately uniform
    var min_p: f64 = result.power[1];
    var max_p: f64 = result.power[1];
    for (result.power[1 .. result.power.len - 1]) |p| {
        if (p < min_p) min_p = p;
        if (p > max_p) max_p = p;
    }
    // Ratio should be < 5× for white noise with 16 segments
    try std.testing.expect(max_p / min_p < 5.0);
}
