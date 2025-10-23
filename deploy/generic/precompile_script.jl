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
dd_json = IMAS.json2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.json"))
dd_hdf5 = IMAS.hdf2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.h5"))

# Run simulation if JSON data loaded
println("Running basic PAM simulation...")
try
    pellet = PAM.run_PAM(dd_json;
        t_start=0.0,
        t_finish=0.0125,
        time_step=0.0001,
        drift_model=:none,
        Bt_dependance=true,
        update_plasma=false
    )
    println("✓ Basic simulation complete")

    # Precompile common operations
    maximum(pellet.density_source)
    println("✓ Common operations precompiled")
catch e
    println("⚠️  Basic simulation failed: $e")
end

println("Precompiling plot recipes...")
try
    # Run the simulation again to get pellet for plotting
    pellet = PAM.run_PAM(dd_json;
        t_start=0.0,
        t_finish=0.0125,
        time_step=0.0001,
        drift_model=:none,
        Bt_dependance=true,
        update_plasma=false
    )

    # Precompile plot code paths
    p = plot(dd_hdf5.equilibrium; cx=true, dpi=300)
    plot!(pellet; plot_cloud=false)
    savefig(p, "test_D3D_pam.png")
    sleep(0.5)
    rm("test_D3D_pam.png"; force=true)
    println("✓ Plot recipes precompiled")
catch e
    println("⚠️  Plot precompilation failed: $e")
end

# Run full test suite for comprehensive precompilation
println("Running full test suite...")
try
    include(joinpath(pkgdir(PAM), "test", "runtests.jl"))
    println("✓ Tests complete")
catch e
    println("⚠️  Tests failed: $e")
end

println("Precompilation complete!")
