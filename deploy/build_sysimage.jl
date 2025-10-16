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
    elseif startswith(arg, "--outdir=")
        global output_dir = split(arg, "=")[2]
    elseif arg in ["--help", "-h"]
        println("""
        Usage: julia deploy/build_sysimage.jl [OPTIONS]

        Options:
          --cpu-target=TARGET    Set CPU target (default: native)
          --outdir=DIR           Set output directory (default: sysimage/)
          --help, -h             Show this help message

        Examples:
          julia deploy/build_sysimage.jl
          julia deploy/build_sysimage.jl --outdir=/path/to/custom/dir
          julia deploy/build_sysimage.jl --cpu-target=generic --outdir=~/pam_build
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

# Determine output directory: use custom path if provided, otherwise default to project/sysimage
if isnothing(output_dir)
    sysimage_dir = joinpath(project_dir, "sysimage")
else
    # Expand ~ and resolve relative paths
    sysimage_dir = abspath(expanduser(output_dir))
end
mkpath(sysimage_dir)
println("Output: $sysimage_dir")

# Use precompile script from deploy directory
precompile_execution_file = joinpath(dirname(@__FILE__), "precompile_script.jl")
if !isfile(precompile_execution_file)
    error("Precompile script not found: $precompile_execution_file")
end
println("Using precompile script: $precompile_execution_file")

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
