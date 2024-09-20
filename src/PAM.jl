module PAM

import IMAS
using Plots
using DifferentialEquations


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

function pellet_position_t(starting_position::Vector{Float64}, velocity_vector::Vector{Float64}, t::Float64, tinj::Float64)
    if t > tinj
        x = starting_position[1] + velocity_vector[1] * (t - tinj)
        y = starting_position[2] + velocity_vector[2] * (t - tinj)
        z = starting_position[3] + velocity_vector[3] * (t - tinj)

    else
        x = starting_position[1]
        y = starting_position[2]
        z = starting_position[3]

    end
    r = sqrt(x^2.0 + y^2.0)
    return r, z
end



#function rpdot(y::Vector{Float64}, t::Float64, ipel::Vector{Float64}, it::Float64, update::Bool)
function rpdot(t::Float64, tinj::Float64)
    if t > tinj
        x, y, z = pellet_position(starting_position, velocity_vector, time, tinj)

        r, z = PAM.projected_coordinate(x, y, z)
    end
    return r, z
end

"""
get_ilayer: This function returns the layer index based on current pellet radius and radii of the pellet players.
"""

function get_ilayer(pellet_radius::Float64, layers_radii::Vector{Any}, nComp::Int) #, pellet_index::Int, t::Float64)
    pos = pellet_radius .- layers_radii

    if maximum(pos) >= 0.0

        layer_index = maximum([i for i in 1:length(pos) if pos[i] >= 0.0])
    else
        layer_index = 1  # change it from 0 as was in Python
    end
    if layer_index > length(layers_radii)   # Do we actually need this? remove -1 as Julia starts counting from 1, while python from 0, change nComp to number of layers??
        layer_index = length(layers_radii)
    end
    return layer_index
end





mutable struct Pellet1{A,T}
    properties::IMAS.pellets__time_slice___pellet
    time::A
    t::T
    Bt::T
    r::A
    z::A
    x::A
    y::A
    ρ::A
    Te::A
    ne::A
    radius::A
    ablation_rate::A
end

function Pellet1(pelt::IMAS.pellets__time_slice___pellet, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, time::Vector{Float64})
    # coordinates of the pellet
    # first point
    R1 = pelt.path_geometry.first_point.r
    Z1 = pelt.path_geometry.first_point.z
    ϕ1 = pelt.path_geometry.first_point.phi
    X1 = R1 * cos(ϕ1)
    Y1 = R1 * sin(ϕ1)
    # second point
    R2 = pelt.path_geometry.second_point.r
    Z2 = pelt.path_geometry.second_point.z
    ϕ2 = pelt.path_geometry.second_point.phi
    X2 = R2 * cos(ϕ2)
    Y2 = R2 * sin(ϕ2)

    starting_position = [X1, Y1, Z1]
    velocity_vector = [X2 - X1, Y2 - Y1, Z2 - Z1] .* pelt.velocity_initial

    # trajectory for all time steps
    x, y, z = pellet_position(starting_position, velocity_vector, time, time[1])
    r, z = PAM.projected_coordinate(x, y, z)

    # rho interpolant
    _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)

    # check the pellet position
    rho_start = RHO_interpolant.(R1, Z1)
    if (rho_start < 1)
        error("Pellet1 starting inside plasma at rho = $rho_start")
    end

    Bt = eqt.global_quantities.magnetic_axis.b_field_tor
    # get plasma profiles in rho coordinates and set to zero outside of the plasma
    ρ = RHO_interpolant.(r, z)
    ne = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density).(ρ)
    ne[ρ.>1.0] .= 0.0
    Te = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(ρ)
    Te[ρ.>1.0] .= 0.0

    radii = cumsum([layer.thickness for layer in pelt.layer])
    @assert pelt.shape.size[1] == radii[end] "The layer's thickness don't add up to the total radius"
    radius = fill(radii[end], size(time))
    ablation_rate = fill(0.0, size(time))

    return Pellet1(pelt, time, time[1], Bt, r, z, x, y, ρ, Te, ne, radius, ablation_rate)
end

function drift(pelt::Pellet1, k::Int)
    # just a simple assumtion now
    dr_drift = 0.5*pelt.radius[k]

    return dr_drift
end

function dr_dt(pelt::Pellet1, k::Int)
    # model for each layer, now works for DT
    layer = pelt.properties.layer[1]

    A = Float64[]
    Z = Float64[]
    fractions = Float64[]
    fractD = 0
    fractT = 0
    AD = 0
    AT = 0
    #  following numbers (rho) are from OMFIT script, are they really correct? 
    ρ_D = 0.2
    ρ_T = 0.318


    for species in layer.species

        tmp = IMAS.ion_properties(Symbol(species.label))
        push!(A, tmp.a)
        push!(Z, tmp.z_n)
        push!(fractions, species.fraction)


        if species.label == "D"

            fractD = species.fraction
            AD = tmp.a
        elseif species.label == "T"

            fractT = species.fraction
            AT = tmp.a
        else
            println("Pelet model works only for DT pellet")
        end
    end
    @assert sum(fractions) == 1.0 "Species fractions dont sum up to 1.0"



    Wratio = (1 - fractD) * AT / AD + fractD
    # following equation is from the paper, not OMFIT script
    ρ_mean = (1 - fractD + fractD * AD / AT) * ((1 - fractD) / ρ_T + (fractD * AD / AT) / ρ_D)^(-1)
    # probably Bt should be something
    Bt_exp = 1.0
    c0 = 8.358 * Wratio^0.6667 * (abs(pelt.Bt) / 2.0)^Bt_exp
    #Wavg = sum(A.*fractions)
    cont=1e-18 # just for testing and plotting
    dr_dt = -c0 / ρ_mean * pelt.Te[k]^1.6667 * pelt.ne[k]^0.3333 / abs(pelt.radius[k])^0.6667*cont
    dr_dt *= 1e-3  # to make everything in meters?
    G = -dr_dt * (4 * pi * ρ_mean * pelt.radius[k] * pelt.radius[k])
    # drpdt *=1e-3   #???? do we need this? cm/s to cm/ms?
    # Gd = -fracD*G*N_A*1e3/Wavg*ZD
    # Gt = -fracT*G*N_A*1e3/Wavg*ZT
    
    return (dr_dt, G)   
end

function pellet_density(pelt::Pellet1, eqt::IMAS.equilibrium__time_slice, k::Int)
    
    
    R=dd.equilibrium.time_slice[].global_quantities.magnetic_axis.r
    Z = dd.equilibrium.time_slice[].global_quantities.magnetic_axis.z
    shiftR = drift(pelt[k])
    rcloudR=pelt.radius[k]*cloudFactorR  # from original pellet dd or should it be in Pellet1?
    rcloudZ=pelt.radius[k]*cloudFactorZ
    nsource = exp(-0.5 * ((pelt.r[k] - R + 0.5 * shiftR) ^ 2 / (rcloudR + 0.25 * shiftR) ^ 2 + (pelt.z[k] - Z) ^ 2 / rcloudZ^2))

    # need to normilize on the surface area under the gaussian shape
     nsource /= (2 * π) ^ 2 * (rcloudR + 0.25 * shiftR) * (pelt.r[k] + 0.5 * shiftR) * rcloudZ
   
    nsource *= pelt.ablation_rate[k] 


    # ----------itegrate over flux surface

    #nsource *=dt # from input ????

    # the total source should go to dd source?

 return nsource
end

# pellet radius function calculates the radius of the pellet based on eq(1) and eq(3) from the paper ... by solving the ODE equation

function pellet_radius!(pelt::Pellet1)

    constA=0.2  # take all values which are constants from eq (1) and eq(3)
    f(u, p, t)=u^2*1.01
    tspan=(pelt.time[1], pelt.time[end])
    dt=pelt.time[2]-pelt.time[1]
    #tspan=(0.0, 1.0)
    u0=pelt.radius[1]
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8, saveat=dt)
    
    pelt.radius = sol[:]
    return pelt
end

function ablate!(pelt::Pellet1)
    for k in 2:length(pelt.time)
        pelt.t = pelt.time[k]
        pelt.radius[k] = pelt.radius[k-1]
        if pelt.radius[k] > 0.0 && pelt.ρ[k] < 1.0

            pelt.radius[k] += dr_dt(pelt, k)[1] 
            
            pelt.ablation_rate[k] = dr_dt(pelt, k)[2]
           
            if pelt.radius[k] < 0.0
                pelt.radius[k] = 0.0
            end
        end
    end
    return pelt
end




@recipe function plot_pellet(pelt::Pellet1)
    deposition = abs.(IMAS.gradient(pelt.radius))
    @series begin
        label := "pellet $(IMAS.index(pelt.properties))"
        linewidth := deposition ./ maximum(deposition) * 5 .+ 0.1
        pelt.r, pelt.z
    end
    @series begin
        primary := false
        seriestype := :scatter
        [pelt.r[1]], [pelt.z[1]]
    end
end

function run_pam(dd::IMAS.dd, inputs)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    # Bt = eqt.global_quantities.magnetic_axis.b_field_tor
    # ### Bt_exp?????????
    # Bt_exp = 1.0


    time = collect(range(inputs.t0, inputs.tf, step=inputs.dt))

    pellet = Pellet1(dd.pellets.time_slice[].pellet[1], eqt, cp1d, time)

    #ablate!(pellet)
    pellet_radius!(pellet)

    return pellet






end

end
