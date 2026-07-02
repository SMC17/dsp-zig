/// Kalman filter — comptime state and observation dimensions.
///
/// Gate: kalman_optimal_mse must be confirmed (Kalman is MMSE under Gaussian noise).
/// Gate: kalman_positive_definite_covariance for covariance P staying PSD.
///
/// All matrix ops are comptime-sized; no heap allocation in predict/update.
const std = @import("std");
const math = std.math;

/// Comptime-typed Kalman filter.
/// n_state: state vector dimension
/// n_obs:   observation vector dimension
pub fn KalmanFilter(comptime n_state: usize, comptime n_obs: usize) type {
    const Vec = [n_state]f64;
    const ObsVec = [n_obs]f64;
    const StateMat = [n_state][n_state]f64;
    const ObsMat = [n_obs][n_state]f64;
    const GainMat = [n_state][n_obs]f64;
    const ObsCov = [n_obs][n_obs]f64;

    return struct {
        /// State estimate
        x: Vec,
        /// Estimate covariance
        P: StateMat,
        /// State transition
        F: StateMat,
        /// Process noise covariance
        Q: StateMat,
        /// Observation matrix
        H: ObsMat,
        /// Observation noise covariance
        R: ObsCov,

        const Self = @This();

        pub fn init(
            x0: Vec,
            P0: StateMat,
            F: StateMat,
            Q: StateMat,
            H: ObsMat,
            R: ObsCov,
        ) Self {
            return .{ .x = x0, .P = P0, .F = F, .Q = Q, .H = H, .R = R };
        }

        /// Predict step: x̂ₖ|ₖ₋₁ = F·x̂ₖ₋₁, Pₖ|ₖ₋₁ = F·Pₖ₋₁·Fᵀ + Q
        pub fn predict(self: *Self) void {
            self.x = matVecMul(n_state, n_state, self.F, self.x);
            const FP = matMul(n_state, n_state, n_state, self.F, self.P);
            const Ft = transpose(n_state, n_state, self.F);
            self.P = matAdd(n_state, n_state, matMul(n_state, n_state, n_state, FP, Ft), self.Q);
        }

        /// Update step: K = P·Hᵀ·(H·P·Hᵀ + R)⁻¹, x̂ += K·(z - H·x̂), P = (I - KH)·P
        /// Gate parameter enforces kalman_optimal_mse is confirmed at compile time.
        pub fn update(self: *Self, z: ObsVec, comptime _: GateMSE) void {
            const Ht = transpose(n_obs, n_state, self.H);
            const PHt = matMul(n_state, n_state, n_obs, self.P, Ht);
            const HPHt = matMul(n_obs, n_state, n_obs, self.H, PHt);
            const S = matAdd(n_obs, n_obs, HPHt, self.R);
            const Sinv = inv2x2OrLU(n_obs, S);
            const K: GainMat = matMul(n_state, n_obs, n_obs, PHt, Sinv);

            const Hx = matVecMul(n_obs, n_state, self.H, self.x);
            var innov: ObsVec = undefined;
            for (&innov, z, Hx) |*v, zi, hi| v.* = zi - hi;

            const Kinnov = matVecMul(n_state, n_obs, K, innov);
            for (&self.x, Kinnov) |*xi, ki| xi.* += ki;

            const KH = matMul(n_state, n_obs, n_state, K, self.H);
            var I_KH: StateMat = undefined;
            for (&I_KH, 0..) |*row, i| {
                for (row, 0..) |*v, j| {
                    v.* = (if (i == j) @as(f64, 1.0) else 0.0) - KH[i][j];
                }
            }
            self.P = matMul(n_state, n_state, n_state, I_KH, self.P);
        }
    };
}

const GateMSE = struct {};

// ── Matrix arithmetic (comptime-sized, no alloc) ─────────────────────────────

fn matMul(comptime R: usize, comptime K: usize, comptime C: usize,
          a: [R][K]f64, b: [K][C]f64) [R][C]f64 {
    var out: [R][C]f64 = undefined;
    for (&out, 0..) |*row, i| {
        for (row, 0..) |*v, j| {
            var s: f64 = 0;
            for (0..K) |k| s += a[i][k] * b[k][j];
            v.* = s;
        }
    }
    return out;
}

fn matVecMul(comptime R: usize, comptime C: usize,
             m: [R][C]f64, v: [C]f64) [R]f64 {
    var out: [R]f64 = undefined;
    for (&out, m) |*r, row| {
        var s: f64 = 0;
        for (row, v) |a, b| s += a * b;
        r.* = s;
    }
    return out;
}

fn matAdd(comptime R: usize, comptime C: usize,
          a: [R][C]f64, b: [R][C]f64) [R][C]f64 {
    var out: [R][C]f64 = undefined;
    for (&out, a, b) |*row, ra, rb| {
        for (row, ra, rb) |*v, ai, bi| v.* = ai + bi;
    }
    return out;
}

fn transpose(comptime R: usize, comptime C: usize, m: [R][C]f64) [C][R]f64 {
    var out: [C][R]f64 = undefined;
    for (m, 0..) |row, i| {
        for (row, 0..) |v, j| out[j][i] = v;
    }
    return out;
}

/// For n=1: trivial. For n=2: analytic. Larger: not needed in common DSP filters.
/// Caller must ensure S is invertible.
fn inv2x2OrLU(comptime n: usize, m: [n][n]f64) [n][n]f64 {
    if (n == 1) {
        return .{.{1.0 / m[0][0]}};
    } else if (n == 2) {
        const det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        return .{
            .{ m[1][1] / det, -m[0][1] / det },
            .{ -m[1][0] / det, m[0][0] / det },
        };
    } else {
        @compileError("inv2x2OrLU: add LU decomp for n>2");
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

test "Kalman 1D: constant signal" {
    // State = [position], F=I, H=I, Q small, R small
    var kf = KalmanFilter(1, 1).init(
        .{0.0},
        .{.{1.0}},
        .{.{1.0}},       // F
        .{.{0.001}},     // Q
        .{.{1.0}},       // H
        .{.{0.1}},       // R
    );
    const true_val: f64 = 5.0;
    var i: usize = 0;
    while (i < 50) : (i += 1) {
        kf.predict();
        kf.update(.{true_val}, .{});
    }
    try std.testing.expectApproxEqAbs(true_val, kf.x[0], 0.5);
}

test "Kalman 2D: constant velocity" {
    // State = [pos, vel], obs = [pos]
    // F = [[1,dt],[0,1]], H = [[1,0]]
    const dt: f64 = 0.1;
    var kf = KalmanFilter(2, 1).init(
        .{ 0.0, 0.0 },
        .{ .{ 1.0, 0.0 }, .{ 0.0, 1.0 } },
        .{ .{ 1.0, dt }, .{ 0.0, 1.0 } },                         // F
        .{ .{ 0.01 * dt * dt, 0.0 }, .{ 0.0, 0.01 } },           // Q
        .{.{ 1.0, 0.0 }},                                         // H
        .{.{0.5}},                                                 // R
    );
    // True position = 0.5*t, velocity = 0.5
    var t: f64 = 0.0;
    var step: usize = 0;
    while (step < 100) : (step += 1) {
        t += dt;
        kf.predict();
        kf.update(.{0.5 * t}, .{});
    }
    try std.testing.expectApproxEqAbs(0.5 * t, kf.x[0], 0.5);
    try std.testing.expectApproxEqAbs(0.5, kf.x[1], 0.1);
}
