const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    _ = b.addModule("dsp_zig", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib = b.addLibrary(.{
        .name = "dsp",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/root.zig"),
            .target = target,
            .optimize = optimize,
        }),
        .linkage = .static,
    });
    b.installArtifact(lib);

    const tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/root.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const run_tests = b.addRunArtifact(tests);
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_tests.step);

    const demo_mod = b.createModule(.{
        .root_source_file = b.path("examples/reference_demo.zig"),
        .target = target,
        .optimize = optimize,
    });
    demo_mod.addImport("dsp", lib.root_module);
    const demo = b.addExecutable(.{ .name = "dsp-reference", .root_module = demo_mod });
    b.installArtifact(demo);
    const run_demo = b.addRunArtifact(demo);
    const demo_step = b.step("demo", "Check FFT and convolution against direct references");
    demo_step.dependOn(&run_demo.step);
}
