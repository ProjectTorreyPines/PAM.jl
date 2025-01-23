module PAM

import IMAS
using Plots
using DifferentialEquations
#using Statistics
#using Sundials


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


mutable struct Pellet1{A,T, N, S, B, X}
    properties::IMAS.pellets__time_slice___pellet
    Btdep::B
    drift_model::S
    time::A
    t::T
    Bt::A    
    velocity_vector::X
    r::A
    R_drift::A
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
dt

function Pellet1(pelt::IMAS.pellets__time_slice___pellet, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, time::Vector{Float64}, surfaces::Vector{IMAS.FluxSurface}, drift_model::Symbol, BtDependance::Bool )
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
    dist = sqrt((X2 - X1)^2 + (Y2 - Y1)^2 + (Z2 - Z1)^2)
    velocity_vector = [X2 - X1, Y2 - Y1, Z2 - Z1] .* pelt.velocity_initial ./(dist)
   
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
   
    Btdep = BtDependance
    #Bt = abs(eqt.eqt.global_quantities.magnetic_axis.r .b_field_tor) .* eqt.global_quantities.magnetic_axis.r ./r
    Bt = abs(eqt.global_quantities.vacuum_toroidal_field.b0) .* eqt.global_quantities.vacuum_toroidal_field.r0 ./r
 
    # get plasma profiles in rho coordinates and set to zero outside of the plasma
    ρ = RHO_interpolant.(r, z)
    ne =  IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density).(ρ)
    ne[ρ.>1.0] .= 0.0
    Te = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(ρ)
    Te[ρ.>1.0] .= 0.0

    radii = cumsum([layer.thickness for layer in pelt.layer])
    @assert pelt.shape.size[1] == radii[end] "The layer's thickness don't add up to the total radius"
    radius = fill(radii[end], size(time))
    ablation_rate = fill(0.0, size(time))
    R_drift = fill(0.0, size(time))
    density_source = fill(0.0, (length(time), length(surfaces)))
    
    return Pellet1(pelt, Btdep, drift_model, time, time[1], Bt , velocity_vector, r, R_drift, z, x, y, ρ, Te, ne, radius, ablation_rate, density_source)
end

"""
get_ilayer: This function returns the layer index based on current pellet radius and thikness of the pellet radii.
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



function drift!(pelt::Pellet1, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d,  k::Int)
    

    Raxis=eqt.global_quantities.magnetic_axis.r
    Zaxis=eqt.global_quantities.magnetic_axis.z
    _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)
    if pelt.drift_model == :Parks
       
       # ---- constants ----------------
        m_p = 1.67262192369e-27
        e = 1.602176634e-19
        mu_0 = 1.25663706212e-6

        mi=m_p*2
        T0=2.0    # eV, temperature at the pellet surface (channel entrance)
        M0 = 1  # Mach number at the channel entrance, from 2007 Roman's paper, in 2000 paper M0=0.8
        c_perp=0.02 # pellet cloud width factor taken as a constant from OMFIT
       #---------------------------------------------------------
        Ca_inf=pelt.Bt[k]/sqrt(mu_0*pelt.ne[k]*mi)
        C_s=sqrt(e*T0/mi)
        
        if pelt.Btdep
            attenuation_factor = pelt.Bt[k]^0.7
        else
            attenuation_factor = 2.0^0.7
        end
       
        rperp_calc=c_perp*sqrt(pelt.ablation_rate[k]/(pelt.ne[k]*pelt.Te[k]^1.5*attenuation_factor))
        
        rperp = max(pelt.radius[1], rperp_calc)
        
        Lc0 = sqrt(Raxis*rperp)
      
        

        pfact_exp = 1.69423774
        Sigma0_exp = 0.85


        loglam=23.5-log((pelt.ne[k]*1e-6)^0.5/pelt.Te[k]^(5/6))
        tau_inf = 2.248882e16*pelt.Te[k]^2/loglam
        a1=0.01*M0
        
        n0_0= pelt.ablation_rate[k]/(2*π*rperp^2*C_s)/(a1*(T0/pelt.ne[k]/pelt.Te[k])^(pfact_exp)*(Lc0/tau_inf)^(Sigma0_exp))
               
        n0=n0_0^(1/(1+pfact_exp+Sigma0_exp))
                    
        Sigma0=n0*Lc0/tau_inf
        pe=pelt.Te[k]*pelt.ne[k]*e
     
        pe_rho = e.*cp1d.electrons.density.*cp1d.electrons.temperature
        

        function vRdot!(du,  u, p, t)
          
          pfact_exp = 1.69423774
          Sigma0_exp = 0.85  
         
          R0 = u[1]
    
          vR0 = u[2]  
                     
          Z0 = pelt.z[k]
         
          tbar = t/Lc0 * C_s
          
          rho_cloud = RHO_interpolant.(R0, Z0)[1]
         
          if rho_cloud > 2.0
            
            du[1] = 0.0
            du[2] = 0.0
          end
           
          
         
            
            peinf =  IMAS.interp1d(cp1d.grid.rho_tor_norm, pe_rho).(rho_cloud)
            
             
            peinf = min(e*n0*T0, max(pe, peinf)) 
             
            Pfact = e*n0*T0/peinf
             


            a1 = 0.52296257 / Lc0/(Sigma0^Sigma0_exp*Pfact^pfact_exp)
            a2 = 17.19090852
            Lc = 1.0+(Sigma0^Sigma0_exp*Pfact^pfact_exp)*tbar
             
           
            PSI = a1*(exp(-a2*(Lc-Lc0)/Lc0/Pfact^(2*pfact_exp)/Sigma0^(2*Sigma0_exp))-1/Pfact)*Lc/Lc0+(1-a1)*(1-1/Pfact)
            PSI = max(0.0, PSI)
           
            if PSI <= 0.0
            
             du[1] = 0.0
             du[2] = 0.0
           
            else 
            
          
            du[1]=vR0
                 
            du[2] = -2 * pelt.Bt[k]^2/Ca_inf/mu_0  * vR0 /(mi*n0*Lc0) +2/R0*PSI*(C_s^2/Lc0)
            
            end
        
        
        end
          
              
        t_eval = pelt.time[k:end] .- pelt.time[k]      
       
        t_span = [0, (t_eval[end]-t_eval[1])]
      
        
      
      
        u0 = [pelt.r[k],pelt.velocity_vector[1]]
        
        prob = ODEProblem(vRdot!, u0, t_span)
       
        sol = solve(prob, RadauIIA5(autodiff=false),verbose=false)
       
        dr_drift= sol[1,end]-pelt.r[k]
        
    elseif pelt.drift_model== :HPI2
        c1 = 0.116
        c2 = 0.120
        c3 = 0.368
        c4 = 0.041
        c5 = 0.015
        c6 = 1.665
        c7 = 0.439
        c8 = 0.217
        c9 = -0.038
        c10 = 0.493
        c11 = 0.193
        c12 = -0.346
        c13 = -0.204
       
       vp=sqrt(pelt.velocity_vector[1]^2+pelt.velocity_vector[2]^2+pelt.velocity_vector[3]^2)
       
       rp=pelt.radius[1]*1e3  # from m to mm
       ne0=pelt.ne[1]*1e-1 #to 1e19 m^-3
       te0=pelt.Te[1]*1e-3
      
       R0=eqt.global_quantities.magnetic_axis.r
       B0=eqt.global_quantities.magnetic_axis.b_field_tor
       a0 = eqt.boundary.minor_radius
       kappa= eqt.boundary.elongation
       alpha=atan(pelt.z[k]-Zaxis, pelt.r[k]-Raxis)


       dr_drift=c1*(vp/100)^c2 *rp^c3*ne0^c4*te0^c5 *(abs(abs(alpha)-c6)+c8)^c7*(1-pelt.ρ[k])^c9*a0^c10*R0^c11*B0^c12*kappa^c13
   
    
    elseif pelt.drift_model== :none
        dr_drift = 0.0

    end
    
    pelt.R_drift[k] = dr_drift
    return 
end

function dr_dt!(pelt::Pellet1, k::Int)
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
    if pelt.Btdep
      Bt_exp=-0.35  #scaling of the magnetic field 
    else 
      Bt_exp=0
    end
    
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
      
       
       ρ_zero = (1 - FD + FD * AD / AT) * ((1 - FD) / pellet_mass_density("T") + (FD * AD / AT) /pellet_mass_density("D"))^(-1) #[g cm^-3]
       
      
       Bt=pelt.Bt[k]
      
       Wratio=(1-FD)*AT/AD+FD
       c0 = 8.358 * Wratio^0.6667 * (abs(Bt) / 2.0) ^ Bt_exp
       
       dr_dt=-c0/ρ_zero*(pelt.Te[k]*1e-3)^(1.6667)*(pelt.ne[k]*1e-20)^(0.3333)/(pelt.radius[k-1]*1e2)^0.6667
       
       pelt.radius[k] = max(0.0, pelt.radius[k-1] + dr_dt * 1e-2 * (pelt.time[k] - pelt.time[k-1]))
       G = -dr_dt * (4 * π * ρ_zero * (pelt.radius[k-1]*1e2)^2)
     
       G *= 6.022e23*FD/A_mean*2
       pelt.ablation_rate[k] = G
      
                  
    else
      println("No ablation model is implemented for such combination of species")
    end
   
   
   return 
end

function pellet_density(pelt::Pellet1, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, surface::IMAS.FluxSurface, k::Int)
    # calculations of the pellet cloud size
    cloudFactor = 5
    cloudFactorR=cloudFactor
    cloudFactorZ=cloudFactor
    
    rcloudR=pelt.radius[k]*cloudFactorR  
    rcloudZ=pelt.radius[k]*cloudFactorZ
   
   
    
   
    # Assume density source as a 2D Gaussian function
    nsource = exp.(-0.5 .* ((pelt.r[k] .- surface.r .+ 0.5 .* pelt.R_drift[k]) .^ 2 ./ (rcloudR .+ 0.25 .* pelt.R_drift[k]) .^ 2 .+ (pelt.z[k] .- surface.z) .^ 2 ./ rcloudZ.^2))

    # need to normilize on the surface area under the gaussian shape
    nsource ./= (2 * π) ^ 2 * (rcloudR + 0.25 * pelt.R_drift[k]) * (pelt.r[k] + 0.5 * pelt.R_drift[k]) * rcloudZ
    
  
    nsource .*= pelt.ablation_rate[k]
   
    return  IMAS.flux_surface_avg(nsource, surface)
end


function ablate!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, pelt::Pellet1, surfaces::Vector{IMAS.FluxSurface})

    pellet_source = zeros(length(pelt.time),length(surfaces))
    #rho_source = IMAS.interp1d(eqt.profiles_1d.psi, eqt.profiles_1d.rho_tor_norm).([surface.psi for surface in surfaces])

    for k in 2:length(pelt.time)
       dt= pelt.time[k] - pelt.time[k-1]
         
      
       if  pelt.ρ[k] > 1.0
           pelt.radius[k] = pelt.radius[k-1]
      
       else
             
            dr_dt!(pelt,k) 
           
            
                 
           
            drift!(pelt, eqt, cp1d, k)

            for (ks, surface) in enumerate(surfaces)
                tmp = pellet_density(pelt, eqt, cp1d, surface, k)
               
                pellet_source[k,ks] += tmp*dt

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

function run_PAM(dd::IMAS.dd, inputs)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    # generate time array for the simulations: t0 start of the pellet modeling, tf - end of the simulations)
    time = collect(range(inputs.t0, inputs.tf, step=inputs.dt))
    
    # define flux surfaces, will be needed for the pellet source calculations
    surfaces = IMAS.trace_surfaces(eqt, IMAS.first_wall(dd.wall)...)
    
    drift_model=inputs.drift_model
    BtDependance=inputs.BtDependance

    # initialize the pellet structure 
    pellet = Pellet1(dd.pellets.time_slice[].pellet[1], eqt, cp1d, time, surfaces, drift_model, BtDependance)
    
    ablate!(eqt,cp1d, pellet, surfaces)
   
    return pellet
end

end
