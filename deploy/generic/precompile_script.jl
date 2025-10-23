#!/usr/bin/env julia
# Precompile script for PAM system image
# This script runs during sysimage creation to trace and precompile frequently-used code paths

# Load packages to trigger precompilation
using PAM
using IMAS
using DifferentialEquations
using Plots

function precompile_pam()
    println("Precompiling PAM packages...")
    println("✓ Packages loaded")

    # Set GR to headless mode BEFORE loading any plots
    ENV["GKSwstype"] = "nul"

    # Try to load example data (graceful fallback if not available)
    println("Loading example data...")
    dd_json = nothing
    dd_hdf5 = nothing

    try
        dd_json = IMAS.json2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.json"))
        println("✓ JSON data loaded")
    catch e
        println("⚠️  Could not load JSON: $e")
    end

    try
        dd_hdf5 = IMAS.hdf2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.h5"))
        println("✓ HDF5 data loaded")
    catch e
        println("⚠️  Could not load HDF5: $e")
    end

    # Run simulation if JSON data loaded
    if !isnothing(dd_json)
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
    else
        println("⚠️  Skipping basic simulation (no JSON data)")
    end

    # Precompile plot recipes if both data sources available
    if !isnothing(dd_json) && !isnothing(dd_hdf5)
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

            # Clean up
            rm("test_D3D_pam.png"; force=true)

            # Close all plots to prevent GR backend from hanging
            closeall()

            println("✓ Plot recipes precompiled")
        catch e
            println("⚠️  Plot precompilation failed: $e")
        end
    else
        println("⚠️  Skipping plot precompilation (missing JSON or HDF5 data)")
    end

    # Run full test suite for comprehensive precompilation
    run_tests = get(ENV, "PAM_PRECOMPILE_TESTS", "false") == "true"
    if run_tests
        println("Running full test suite (PAM_PRECOMPILE_TESTS=true)...")
        try
            include(joinpath(pkgdir(PAM), "test", "runtests.jl"))
            println("✓ Tests complete")
        catch e
            println("⚠️  Tests failed: $e")
        end
    else
        println("⚠️  Skipping full test suite (set PAM_PRECOMPILE_TESTS=true to enable)")
    end

    println("Precompilation complete!")
end

# Execute precompilation
precompile_pam()
