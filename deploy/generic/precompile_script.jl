#!/usr/bin/env julia
# Precompile script for PAM system image
# This script runs during sysimage creation to trace and precompile frequently-used code paths

println("Precompiling PAM packages...")

# Load packages to trigger precompilation
using PAM
using IMAS
using DifferentialEquations
using Plots

println("✓ Packages loaded")

# Try to load example data (graceful fallback if not available)
println("Loading example data...")
dd = nothing
try
    dd = IMAS.json2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.json"))
    println("✓ Example data loaded from JSON")
catch e
    println("⚠️  Could not load example data: $e")
    println("   Skipping precompile simulation (package imports still precompiled)")
end

# Only run simulation if data loaded successfully
if !isnothing(dd)
    println("Running PAM simulation...")
    try
        pellet = PAM.run_PAM(dd;
            t_start=0.0,
            t_finish=0.0125,
            time_step=0.0001,
            drift_model=:none,
            Bt_dependance=true,
            update_plasma=false
        )
        println("✓ Simulation complete")

        # Precompile common operations
        maximum(pellet.density_source)
        println("✓ Common operations precompiled")
    catch e
        println("⚠️  Simulation failed: $e")
        println("   (Package code still precompiled)")
    end
end

println("Precompilation complete!")
