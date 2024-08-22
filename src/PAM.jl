module PAM

import IMAS
using Plots



function pellet_position(starting_position::Vector{Float64}, velocity_vector::Vector{Float64}, time::AbstractVector, tinj::Float64)
    x = starting_position[1] .+ velocity_vector[1] .* (time .- tinj)
    y = starting_position[2] .+ velocity_vector[2] .* (time .- tinj)
    z = starting_position[3] .+ velocity_vector[3] .* (time .- tinj)
    return x, y, z
end

function projected_coordinate(x::AbstractVector{Float64}, y::AbstractVector{Float64}, z::AbstractVector{Float64})
    r = sqrt.(x .^ 2.0 .+ y .^ 2.0)
    return r, z
end

function run_pam(dd::IMAS.dd, inputs)

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    pelt = dd.pellets.time_slice[].pellet

    _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)

    # coordinates of the pellet

    R1 = pelt[1].path_geometry.first_point.r
    Z1 = pelt[1].path_geometry.first_point.z
    ϕ1 = pelt[1].path_geometry.first_point.phi
    X1 = R1 * cos(ϕ1)
    Y1 = R1 * sin(ϕ1)

    R2 = pelt[1].path_geometry.second_point.r
    Z2 = pelt[1].path_geometry.second_point.z
    ϕ2 = pelt[1].path_geometry.second_point.phi
    X2 = R2 * cos(ϕ2)
    Y2 = R2 * sin(ϕ2)

    starting_position = [X1, Y1, Z1]
    velocity_vector = [X2 - X1, Y2 - Y1, Z2 - Z1] .* pelt[1].velocity_initial

    # --- check the pellet position
    rho_start = RHO_interpolant.(R1, Z1)
    if (rho_start < 1)
        error("Pellet starting inside plasma at rho = $rho_start")
    end

    #---- position of the pellet on each time step (times will be writtent from dd)----------
    time = range(inputs.t0, inputs.tf, step=inputs.dt)
    tinj = time[1] #dd.pulse_schedule.pellets.???   # will be taken from pulse_schedule

    x, y, z = pellet_position(starting_position, velocity_vector, time, tinj)

    r, z = PAM.projected_coordinate(x, y, z)

    # plot the pellet direction in time

    if inputs.debug
        plot(dd.equilibrium; cx=true)
        plot!(r, z)
        display(scatter!([r[1]], [z[1]]))

        rmin, rmax = extrema(dd.equilibrium.time_slice[].boundary.outline.r)
        plot(x, y)
        t = 0:0.01:2π
        plot!(rmin .* cos.(t), rmin .* sin.(t))
        display(plot!(rmax .* cos.(t), rmax .* sin.(t); aspect_ratio=:equal))
    end

    # get plasma profiles in rho coordinates and applay mask outside of plasma

    ρ = RHO_interpolant.(r, z)
    mask = zero(ρ)
    mask[ρ.<1.0] .= 1.0

    ne = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density).(ρ).*mask
    Te = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(ρ).*mask

    if inputs.debug
        display(scatter(ρ, ne; alpha=0.1))
        display(scatter(ρ, Te; alpha=0.1))
    end

    
end

end
