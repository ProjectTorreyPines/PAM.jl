#!/usr/bin/env julia
# Precompile script for PAM system image
# This script runs during sysimage creation to trace and precompile frequently-used code paths

using PAM, IMAS, DifferentialEquations, Plots

println("Loading example data...")
dd_json = IMAS.json2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.json"))
dd_hdf5 = IMAS.hdf2imas(joinpath(pkgdir(PAM), "examples", "template_D3D_1layer_2species.h5"))

println("Running PAM simulation...")
# Match actual usage in PAM_test.jl (t_finish=0.0125)
pellet = PAM.run_PAM(dd_json;
    t_start=0.0,
    t_finish=0.0125,
    time_step=0.0001,
    drift_model=:none,
    Bt_dependance=true,
    update_plasma=false
)

# Plot recipes
p=plot(dd_hdf5.equilibrium; cx=true, dpi=300)
plot!(pellet; plot_cloud=false)
savefig(p,"test_D3D_pam.png")
sleep(5)
rm("test_D3D_pam.png")


println("Precompiling CLI driver patterns...")
# Simulate typical CLI usage patterns
maximum(pellet.density_source)  # Common output operation
println("Maximum density source: $(maximum(pellet.density_source))")

println("Running tests...")
include(joinpath(pkgdir(PAM), "test", "runtests.jl"))

println("Precompilation complete!")
