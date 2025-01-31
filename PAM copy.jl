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

"""
get_ilayer: This function returns the layer index based on current pellet radius and thikness radii of the pellet radii.
"""

function get_ilayer(pelt::Pellet1,k::Int) 
   
    layers=cumsum([layer.thickness for layer in pelt.properties.layer])
   
    pos = pelt.radius[k].-layers
 
    if maximum(pos) >= 0.0

        layer_index = maximum([i for i in 1:length(pos) if pos[i] >= 0.0])
    else
        layer_index = 1  
    end
       
    return layer_index
end



# the number of g/m^-3 # the data of mass density and atomic weight is sourced from PAM model within OMFIT
function pellet_mass_density(species::String)
    material_density = Dict("DT" => 0.257, "D" => 0.2, "T" => 0.318, "C" => 3.3, "Ne" => 1.44)
    return material_density[species]
end



function drift(pelt::Pellet1, k::Int)
    # just a simple assumtion now
    dr_drift = 0.0

    return dr_drift
end

function dr_dt(pelt::Pellet1, k::Int)
    # model for each layer, now works for DT only
    ilayer= get_ilayer(pelt, k)
    layer = pelt.properties.layer[ilayer]

    A = Float64[]
    Z = Float64[]
    fractions = Float64[]
    species_list = String[]
      

    for species in layer.species
        tmp = IMAS.ion_properties(Symbol(species.label))
        push!(A, tmp.a)
        push!(Z, tmp.z_n)
        push!(fractions, species.fraction)
        push!(species_list,species.label)      
   
    end
    A_mean=sum(fractions.*A)
    @assert sum(fractions) == 1.0 "Species fractions dont sum up to 1.0"
    Bexp=39*(2/pelt.Bt)^0.35  #scaling of the magnetic field 
    
    # ablation model for the DT pellet
    if  ("D" in species_list) & ("T" in species_list)       
               
       if species_list[1]=="D"
         AD=A[1]
         AT=A[2]
         FD=fractions[1]
         
       else
         AD=A[2]
         AT=A[1]
         FD=fractions[2]
         
       end
       # according equation #1 in J.McClenaghan at al Nucl.Fusion 63 (2023) 036015 (two coefficients are missing) and normalization to FUSE units
       G=Bexp*(A_mean/AD)^(2/3)*(1e+2*pelt.radius[k]/0.2)^(4/3)*(pelt.Te[k]*1e-3/2)^(5/3)*(pelt.ne[k]*1e-20)^(1/3)  #[g/s]
       ρ_zero = (1 - FD + FD * AD / AT) * ((1 - FD) / pellet_mass_density("T") + (FD * AD / AT) /pellet_mass_density("D"))^(-1) #[g cm^-3]
       dr_dt =  - G / (4 * pi * ρ_zero * (pelt.radius[k]*1e+2)^2) #[cm/s]
       
    else
      println("No ablation model is implemented for such combination of species")
    end
    
    # normilize return variables to make in standart FUSE usnits
   return (radius_variation=dr_dt*1e-2, ablation_rate=G*1e-3)   
end

function pellet_density(pelt::Pellet1, surface::IMAS.FluxSurface, k::Int)
    # calculations of the pellet cloud size
    cloudFactor = 0.1
    cloudFactorR=cloudFactor
    cloudFactorZ=cloudFactor
    
    rcloudR=pelt.radius[k]*cloudFactorR  
    rcloudZ=pelt.radius[k]*cloudFactorZ
    #  shift of the pellet and change of the cloud due to ExB drift 
    shiftR = drift(pelt,k)
    
    # Assume density source as a 2D Gaussian function
    nsource = exp.(-0.5 .* ((pelt.r[k] .- surface.r .+ 0.5 .* shiftR) .^ 2 ./ (rcloudR .+ 0.25 .* shiftR) .^ 2 .+ (pelt.z[k] .- surface.z) .^ 2 ./ rcloudZ.^2))

    # need to normilize on the surface area under the gaussian shape
    nsource ./= ((2 * π) ^ 2 * (rcloudR + 0.25 * shiftR) * (pelt.r[k] + 0.5 * shiftR) * rcloudZ)
    nsource .*= pelt.ablation_rate[k] 

    return IMAS.flux_surface_avg(nsource, surface)
end


function ablate!(eqt::IMAS.equilibrium__time_slice, pelt::Pellet1, surfaces::Vector{IMAS.FluxSurface})

    pellet_source = zeros(length(pelt.time),length(surfaces))
    rho_source = IMAS.interp1d(eqt.profiles_1d.psi, eqt.profiles_1d.rho_tor_norm).([surface.psi for surface in surfaces])
    


    for k in 2:length(pelt.time)
        dt = pelt.time[k] - pelt.time[k-1]
             

        if pelt.radius[k-1] < 0.0 || pelt.ρ[k] > 1.0
            pelt.radius[k] = pelt.radius[k-1]

        else
            pelt.radius[k] = pelt.radius[k-1] + dr_dt(pelt, k).radius_variation * dt
       
            pelt.ablation_rate[k] = dr_dt(pelt, k).ablation_rate 

            
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
    # generate time array for the simulations: t0 start of the pellet modeling, tf - end of the simulations)
    time = collect(range(inputs.t0, inputs.tf, step=inputs.dt))

    # define flux surfaces, will be needed for the pellet source calculations
    surfaces = IMAS.trace_surfaces(eqt, IMAS.first_wall(dd.wall)...)

    # initialize the pellet structure 
    pellet = Pellet1(dd.pellets.time_slice[].pellet[1], eqt, cp1d, time, surfaces)

    pellet_source = ablate!(eqt, pellet, surfaces)

  

   
    return pellet
end

end
