#!/usr/bin/env julia
# PAM SystemImage Builder
# Minimal script based on FUSE install_pam.jl logic

println("PAM SystemImage Builder")

# Parse command-line arguments
cpu_target = get(ENV, "JULIA_CPU_TARGET", "native")
output_dir = nothing

for arg in ARGS
    if startswith(arg, "--cpu-target=")
        global cpu_target = split(arg, "=")[2]
    elseif startswith(arg, "--output=") || startswith(arg, "-o=")
        global output_dir = split(arg, "=")[2]
    elseif arg in ["--help", "-h"]
        println("""
        Usage: julia deploy/build_sysimage.jl [OPTIONS]

        Options:
          --cpu-target=TARGET    Set CPU target (default: native)
          --output=DIR, -o=DIR   Set output directory (default: sysimage/)
          --help, -h             Show this help message

        Examples:
          julia deploy/build_sysimage.jl
          julia deploy/build_sysimage.jl --output=/path/to/custom/dir
          julia deploy/build_sysimage.jl --cpu-target=generic --output=~/pam_build
        """)
        exit(0)
    end
end

println("CPU Target: $cpu_target")
!isnothing(output_dir) && println("Output: $output_dir")

# Setup environment
import Pkg

project_dir = dirname(dirname(@__FILE__))
Pkg.activate(project_dir)

if !haskey(Pkg.project().dependencies, "PackageCompiler")
    println("Installing PackageCompiler...")
    Pkg.add("PackageCompiler")
end
using PackageCompiler

println("Creating precompile script...")
# Determine output directory: use custom path if provided, otherwise default to project/sysimage
if isnothing(output_dir)
    sysimage_dir = joinpath(project_dir, "sysimage")
else
    # Expand ~ and resolve relative paths
    sysimage_dir = abspath(expanduser(output_dir))
end
mkpath(sysimage_dir)
println("Output: $sysimage_dir")

precompile_execution_file = joinpath(sysimage_dir, "precompile_script.jl")
precompile_cmds = """
# Precompile script for PAM
using PAM, IMAS, DifferentialEquations, Plots

println("Loading example data...")
dd_json = IMAS.json2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.json"))
dd_hdf5 = IMAS.hdf2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.h5"))

println("Running PAM simulation...")
# Match actual usage in PAM_test.jl (t_finish=0.0125)
pellet = PAM.run_PAM(dd_json;
    t_start=0.0,
    t_finish=0.0125,  # ← Changed to match actual usage
    time_step=0.0001,
    drift_model=:none,
    Bt_dependance=true,
    update_plasma=false
)

println("Precompiling CLI driver patterns...")
# Simulate typical CLI usage patterns
maximum(pellet.density_source)  # Common output operation
println("Maximum density source: \$(maximum(pellet.density_source))")

println("Running tests...")
include(joinpath(pkgdir(PAM), "test", "runtests.jl"))

println("Precompilation complete!")
"""
write(precompile_execution_file, precompile_cmds)

println("Building sysimage (5-10 minutes)...")
sysimage_ext = Sys.isapple() ? "dylib" : (Sys.iswindows() ? "dll" : "so")
sysimage_path = joinpath(sysimage_dir, "sys_pam.$sysimage_ext")
create_sysimage(
    ["PAM", "IMAS", "DifferentialEquations", "Plots"];
    sysimage_path,
    precompile_execution_file,
    cpu_target,    
    incremental=false,
)

println("Creating launcher script...")
launcher_script = joinpath(sysimage_dir, "pam")
launcher_content = """
#!/bin/bash
SYSIMAGE_DIR="$sysimage_dir"
PROJECT_DIR="$project_dir"
julia --project=\$PROJECT_DIR --sysimage=\$SYSIMAGE_DIR/sys_pam.$sysimage_ext "\$@"
"""
write(launcher_script, launcher_content)
chmod(launcher_script, 0o755)

println("\n✅ Build complete!")
println("Usage: $launcher_script [script.jl] or julia --project=$project_dir --sysimage=$sysimage_path")
