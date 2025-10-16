#!/usr/bin/env julia
# PAM SystemImage Builder
# Minimal script based on FUSE install_pam.jl logic

println("=" ^ 70)
println("PAM SystemImage Builder")
println("=" ^ 70)

# Parse command-line arguments
cpu_target = get(ENV, "JULIA_CPU_TARGET", "native")
if length(ARGS) > 0 && startswith(ARGS[1], "--cpu-target=")
    cpu_target = split(ARGS[1], "=")[2]
end

println("\n### Configuration")
println("    CPU Target: $cpu_target")
println("    Project: PAM")

# Setup environment
import Pkg

println("\n### Activating PAM environment")
project_dir = dirname(dirname(@__FILE__))
Pkg.activate(project_dir)

println("\n### Installing PackageCompiler")
if !haskey(Pkg.project().dependencies, "PackageCompiler")
    Pkg.add("PackageCompiler")
end
using PackageCompiler

println("\n### Creating precompile script")
sysimage_dir = joinpath(project_dir, "sysimage")
mkpath(sysimage_dir)

precompile_execution_file = joinpath(sysimage_dir, "precompile_script.jl")
precompile_cmds = """
# Precompile script for PAM
using PAM, IMAS, DifferentialEquations, Plots

println("Loading example data...")
dd = IMAS.json2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.json"))

println("Running PAM simulation...")
pellet = PAM.run_PAM(dd;
    t_start=0.0,
    t_finish=0.0045,
    time_step=0.0001,
    drift_model=:none,
    Bt_dependance=true,
    update_plasma=false
)

println("Running tests...")
include(joinpath(pkgdir(PAM), "test", "runtests.jl"))

println("Precompilation complete!")
"""
write(precompile_execution_file, precompile_cmds)
println("    Precompile script: $precompile_execution_file")

println("\n### Building PAM sysimage (this may take 5-10 minutes)...")
# Platform-specific extension: .dylib (macOS), .so (Linux), .dll (Windows)
sysimage_ext = Sys.isapple() ? "dylib" : (Sys.iswindows() ? "dll" : "so")
sysimage_path = joinpath(sysimage_dir, "sys_pam.$sysimage_ext")
println("    Platform: $(Sys.KERNEL)")
println("    Output: sys_pam.$sysimage_ext")
println("\n    Build process:")
println("      1. Precompiling packages → Julia IR (1-2 min)")
println("      2. Executing precompile_script.jl → function tracing (2-3 min)")
println("      3. Compiling traced functions → native code (2-5 min)")
println("      4. Bundling into sysimage file")
create_sysimage(
    ["PAM", "IMAS", "DifferentialEquations", "Plots"];
    sysimage_path,
    precompile_execution_file,
    cpu_target
)

println("\n### Creating launcher script")
launcher_script = joinpath(sysimage_dir, "pam")
launcher_content = """
#!/bin/bash
# PAM launcher with precompiled sysimage
SCRIPT_DIR="\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="\$(dirname "\$SCRIPT_DIR")"

echo "Starting Julia with PAM sysimage..."
julia --project=\$PROJECT_DIR --sysimage=\$SCRIPT_DIR/sys_pam.$(Sys.isapple() ? "dylib" : "so") "\$@"
"""
write(launcher_script, launcher_content)
chmod(launcher_script, 0o755)

println("\n" * "=" ^ 70)
println("✅ PAM SystemImage build complete!")
println("=" ^ 70)
println("\nUsage:")
println("  1. Direct usage:")
println("     julia --project=. --sysimage=sysimage/sys_pam.$(Sys.isapple() ? "dylib" : "so")")
println("\n  2. Using launcher script:")
println("     ./sysimage/pam")
println("\n  3. Interactive REPL:")
println("     ./sysimage/pam -i")
println("\nExpected speedup: ~10x faster package loading")
println("=" ^ 70)
