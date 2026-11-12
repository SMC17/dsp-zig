const std = @import("std");
const dsp = @import("dsp");
const Complex = std.math.Complex(f64);
const tolerance = 1e-9;

fn close(actual: f64, expected: f64) !void {
    if (!std.math.isFinite(actual) or @abs(actual - expected) > tolerance)
        return error.ReferenceMismatch;
}

pub fn main() !void {
    var total_bins: usize = 0;
    inline for (.{ 1, 2, 4, 8, 16 }) |n| {
        var input: [n]Complex = undefined;
        for (&input, 0..) |*x, k| {
            const t: f64 = @floatFromInt(k);
            x.* = .{ .re = @sin(t * 0.37) + @as(f64, @floatFromInt(k % 3)), .im = @cos(t * 0.23) };
        }
        var transformed = input;
        dsp.fft.fft(&transformed);
        // Direct O(N^2) DFT, independent of radix-2 butterflies and permutation.
        for (transformed, 0..) |actual, k| {
            var re: f64 = 0;
            var im: f64 = 0;
            for (input, 0..) |x, j| {
                const angle = -2.0 * std.math.pi * @as(f64, @floatFromInt(k * j)) / @as(f64, n);
                re += x.re * @cos(angle) - x.im * @sin(angle);
                im += x.re * @sin(angle) + x.im * @cos(angle);
            }
            try close(actual.re, re);
            try close(actual.im, im);
            total_bins += 1;
        }
        dsp.fft.ifft(&transformed);
        for (input, transformed) |want, got| {
            try close(got.re, want.re);
            try close(got.im, want.im);
        }
    }
    const allocator = std.heap.page_allocator;
    const left = [_]f64{ 1, -2, 3 };
    inline for (.{ &[_]f64{ 4, 5 }, &[_]f64{ -1, 0, 2 } }) |right| {
        const actual = try dsp.fft.convolve(.{}, allocator, &left, right);
        defer allocator.free(actual);
        var expected = [_]f64{0} ** (left.len + right.len - 1);
        for (left, 0..) |a, i| for (right, 0..) |b, j| {
            expected[i + j] += a * b;
        };
        for (actual, expected) |got, want| try close(got, want);
    }
    std.debug.print("DSP_FIXTURE_OK dft_bins={d} inverse_bins={d} convolution_cases=2 tolerance=1e-9\n", .{ total_bins, total_bins });
}
