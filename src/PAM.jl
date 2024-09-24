module PAM

import IMAS
using Plots
using DifferentialEquations
using Statistics

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

function get_ilayer(pellet_radius::Float64, layers_radii::Vector{Any}, nComp::Int) 
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





mutable struct Pellet1{A,T, N}
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
    density_source::N
end

# the number of g/m^-3 # the data of mass density and atomic weight is sourced from PAM model within OMFIT
function pellet_mass_density(species::String)
    material_density = Dict("DT" => 0.257, "D" => 0.2, "T" => 0.318, "C" => 3.3, "Ne" => 1.44)
    return material_density[species] * 1E-12
end


function Pellet1(pelt::IMAS.pellets__time_slice___pellet, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, time::Vector{Float64}, surfaces::Vector{IMAS.FluxSurface} )
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
    density_source = fill(0.0, (length(time), length(surfaces)))
    return Pellet1(pelt, time, time[1], Bt, r, z, x, y, ρ, Te, ne, radius, ablation_rate, density_source)
end

function drift(pelt::Pellet1, k::Int)
    # just a simple assumtion now
    dr_drift = 0.0

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
    #ρ_pellet=Float64[]

    for species in layer.species
        tmp = IMAS.ion_properties(Symbol(species.label))
        push!(A, tmp.a)
        push!(Z, tmp.z_n)
        push!(fractions, species.fraction)
       # push!(ρ_pellet, pellet_mass_density(Symbol(species.label)) )

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
    @assert ρ_mean > 0.0
   
   
    A_mean=AT*fractT+AD*fractD
    
    # following equation is from the paper, not OMFIT script
    ρ_zero = (1 - fractD + fractD * AD / AT) * ((1 - fractD) / ρ_T + (fractD * AD / AT) / ρ_D)^(-1)
   
    G=39*(1e+2*pelt.radius[k]/0.2)^(4/3)*(A_mean/AD)^(2/3)*(pelt.Te[k].*1e-3).^(5/3)* (pelt.ne[k]*1e-26)^(1/3)
    
    dr_dt = - G / (4 * pi * ρ_mean * (pelt.radius[k]*1e+2)^2)
    
    return (radius_variation=dr_dt, ablation_rate=G)   
end

function pellet_density(pelt::Pellet1, surface::IMAS.FluxSurface, k::Int)
    cloudFactor = 0.1
    cloudFactorR=cloudFactor
    cloudFactorZ=cloudFactor
    
    shiftR = drift(pelt,k)
    rcloudR=pelt.radius[k]*cloudFactorR  # from original pellet dd or should it be in Pellet1?
    rcloudZ=pelt.radius[k]*cloudFactorZ
    nsource = exp.(-0.5 .* ((pelt.r[k] .- surface.r .+ 0.5 .* shiftR) .^ 2 ./ (rcloudR .+ 0.25 .* shiftR) .^ 2 .+ (pelt.z[k] .- surface.z) .^ 2 ./ rcloudZ.^2))

    # need to normilize on the surface area under the gaussian shape
    nsource ./= ((2 * π) ^ 2 * (rcloudR + 0.25 * shiftR) * (pelt.r[k] + 0.5 * shiftR) * rcloudZ)
    nsource .*= pelt.ablation_rate[k] 

    return IMAS.flux_surface_avg(nsource, surface)
end


function ablate!(eqt::IMAS.equilibrium__time_slice, pelt::Pellet1, surfaces::Vector{IMAS.FluxSurface})

    pellet_source = zeros(length(pelt.time),length(surfaces))
    rho_source = IMAS.interp1d(eqt.profiles_1d.psi, eqt.profiles_1d.rho_tor_norm).([surface.psi for surface in surfaces])
    @show length(pellet_source)


    for k in 2:length(pelt.time)
        dt = pelt.time[k] - pelt.time[k-1]
        pelt.t = pelt.time[k]
        

        if pelt.radius[k-1] == 0.0 || pelt.ρ[k] > 1.0
            pelt.radius[k] = pelt.radius[k-1]

        else
            pelt.radius[k] = pelt.radius[k-1] + dr_dt(pelt, k).radius_variation * dt
            
            pelt.ablation_rate[k] = dr_dt(pelt, k).ablation_rate * dt

            
            for (ks, surface) in enumerate(surfaces)
                tmp = pellet_density(pelt, surface, k)
                pellet_source[k,ks] += tmp
            end
          
       
         
        end
    pelt.density_source = pellet_source
    end

    return 
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
   

    time = collect(range(inputs.t0, inputs.tf, step=inputs.dt))

    surfaces = IMAS.trace_surfaces(eqt, IMAS.first_wall(dd.wall)...)

    pellet = Pellet1(dd.pellets.time_slice[].pellet[1], eqt, cp1d, time, surfaces)

    pellet_source = ablate!(eqt, pellet, surfaces)

  

   
    return pellet
end

end
