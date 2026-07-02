/// dsp-zig — Formally-gated digital signal processing.
///
/// Every function carrying a proof gate annotation will fail to compile if the
/// backing mathematical claim has been formally refuted in Lean 4.
/// Zero runtime cost: Gate(claim) is a zero-size struct.
const std = @import("std");
const math = std.math;
const Complex = std.math.Complex;

pub const fft = @import("fft.zig");
pub const kalman = @import("kalman.zig");
pub const filter = @import("filter.zig");
pub const window = @import("window.zig");
pub const psd = @import("psd.zig");

test {
    std.testing.refAllDecls(@This());
}
