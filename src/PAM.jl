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

function ablation_model(pellet_radius::Float64, rho_pellet::Float64, Bt::Float64, Bt_exp::Float64, layer_index::Int, ne::Vector{Float64}, Te::Vector{Float64}, model::String)
    # Bt_exp?????
    #------------ initialize everything with zeros------------------
    G = Gd = Gt = GC = GNe = GLi = 0.0  # should rename to ablation rate?
    drpdt = 0.0
    N_A = 6.022140857e23




    if (rho_pellet <= 1.0) & (pellet_radius > 0.0)

        println("Ablation model is applyed for the layer ", layer_index)
        """
        if model=="DT"
         Wratio = (1-fracD)*WT/WD + fracD
         Wavg = fracT*WT+fracD*WD
         ρ_D=0.2   #!! g/cm^3
         ρ_T = 0.318  # g/cm^3
         ρ_mean= (1-fracD+fracD*WD/WT)*((1-fracD)/ρ_T+(fracD*WD/WT)/ρ_D)^(-1) 
         c0 = 8.358*Wratio^0.6667*(abs(Bt)/2.0)^Bt_exp

         if pellet_radius<0
            drpdt=0
         else
            drpdt = -c0/ρ_mean*Te_interp^1.6667*ne_interp^0.3333/pellet_radius^0.6667
         # ne, Te  - interpolation ...... to the pellet trajectory; just need Te, ne (\rho) ?
            G = -drpdt*(4*pi*ρ_mean*pellet_radius*pellet_radius)
            drpdt *=1e-3   #???? do we need this? cm/s to cm/ms?
            Gd = -fracD*G*N_A*1e3/Wavg*ZD
            Gt = -fracT*G*N_A*1e3/Wavg*ZT
         end
        end 

        if model=="NeD"
          Wratio = (1 - fractD)*WNe/WD + fracD
          Wavg = fracNe*WNe+fracD*WD
          ρ_D=0.2
          ρ_Ne = 1.44
          ρ_mean= (1-fracD+fracD*WD/WNe)*((1-fracD)/ρ_Ne+(fracD*WD/WNe)/ρ_D)^(-1) 
          X = fracD/(2-fracD)
          AoX = 27.0 + tan(1.48*X)
          c0 = AoX/(4*pi)*(abs(Bt)/2.0)^Bt_exp
          
          if pellet_radius<0.0
            drpdt=0.0
          else
            drpdt = -c0/ρ_mean*Te_interp^1.6667*ne_interp^0.3333/pellet_radius^0.6667
         # ne, Te  - interpolation ...... to the pellet trajectory; just need Te, ne (\rho) ?
            G = -drpdt*(4*pi*ρ_mean*pellet_radius*pellet_radius)
            drpdt *=1e-3
            Gd = -fracD*G*N_A*1e3/Wavg*ZD
            GNe = -fracNe*G*N_A*1e3/Wavg*ZNe
          end
          
        end

        if model=="C"
            # need to normilize density and temperature?
            Te_interp *=1.0e3
            ne_interp *=1.0e14
            C0= 8.146777e-9
            gamma = 5.0/3.0   #????
            ZstarPlus1C = 2.86
            Albedo = 23.920538030089528*log(1+0.20137080524063228*ZstarPlusC1)
            flelectro = exp(-1.936)
            fL = (1.0 - Albedo/100.0)*flelectro
            IstC = 60.0
            if Te>30.0
               Ttmp = Te
            else
                Ttmp = 30.0
            end
            loglaCSlow=log(2.0 * Ttmp/IstC * sqrt(ℯ*2.0)) 
            Blamdaq = 1/(ZC12*loglaCSlow)*(4/(2.5 + 2.2 * sqrt(ZstarPlus1C)))
            Gpr = (C0 * WC12^(2.0/3.0)*(gamma - 1)^(1.0/3.0)*(fL*ne)^(1.0/3.0)*pellet_radius^(4.0/3.0)*Te_interp^(11.0/6.0)*BLamdaq^(2.0/3.0))
            xiexp = 0.601
            lamdaa = 0.0933979540623963
            lamdab = -0.7127242270013098
            lamdac = -0.2437544205933372
            lamdad = -0.8534855445478313
            av = 10.420403555938629 * (Ttmp / 2000.0)^ lamdaa
            bv = 0.6879779829877795 * (Ttmp / 2000.0)^ lamdab
            cv = 1.5870910225610804 * (Ttmp / 2000.0)^ lamdac
            dv = 2.9695640286641840 * (Ttmp / 2000.0)^ lamdad
            fugCG = 0.777686
            CG = fugCG * av *log(1 + bv * (ne_interp / 1e14)^( 2.0 / 3.0) * ycm^( 2.0 / 3.0))/ log(cv + dv * (ne_interp / 1e14)^(2.0 / 3.0) * (pellet_radius^( 2.0 /3.0)))
            G = xiexp *CG*Gpr*(2.0/Bt)^Bt_exp
            ρ_mean = ???


            if pellet_radius<0.0
              drpdt=0.0
            else
              drpdt = -G/(4.0 *pi*ρ_mean*pellet_radius*pellet_radius)
              G = -drpdt*(4*pi*ρ_mean*pellet_radius*pellet_radius)
              drpdt *=1e-3
              GC = -fracC12*G*N_A*1e3/WC12*ZC12
            end
        end

        if model=="Li"
            Te_interp *=1.0e3
            ne_interp *=1.0e14

            rabl = 2.7e13*(Te_interp^1.562)*ne_interp^0.497*pellet_radius^1.497

            G = rabl*WLi/6.022e23
            # ρ_mean = ???
            if pellet_radius < 0.0
                drpdt = 0.0
              else
                drpdt = -G/(4.0 *pi*ρ_mean*pellet_radius*pellet_radius)
                G = -drpdt*(4*pi*ρ_mean*pellet_radius*pellet_radius)
                drpdt *=1e-3
                GLi = -fracLi*G*N_A*1e3/WLi*ZLi
              end

        end



        if model=="Parks"
           Temin = 1.0e-3
           if pellet_radius <0.0 | Te_interp <Temin
              drpdt = 0.0
           end
           
           Ihyd = 7.514
           loglamH = 1
           drpdt= -8.2e15  ### ???? constant?

        end
        """

    else
        println("Ablation is zero, pellet outside of plasma")


    end
    return G, Gd, Gt, GC, GNe, GLi, drpdt
end

function update_pellet_densities(drhodR::AbstractVector, drhodZ::AbstractVector, rho2ds::AbstractVector, iteration::Int, dt::Float64, model::String)

    if model == "2DGaussian"


    end

    if model == "RadialGaissian"

    end

end

#

mutable struct Pellet1{A,T}
    properties::IMAS.pellets__time_slice___pellet
    time::A
    t::T
    r::A
    z::A
    x::A
    y::A
    ρ::A
    Te::A
    ne::A
    radius::A
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

    # get plasma profiles in rho coordinates and set to zero outside of the plasma
    ρ = RHO_interpolant.(r, z)
    ne = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density).(ρ)
    ne[ρ.>1.0] .= 0.0
    Te = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(ρ)
    Te[ρ.>1.0] .= 0.0

    radii = cumsum([layer.thickness for layer in pelt.layer])
    @assert pelt.shape.size[1] == radii[end] "The layer's thickness don't add up to the total radius"
    radius = fill(radii[end], size(time))
    radius = time[end] .- time

    return Pellet1(pelt, time, time[1], r, z, x, y, ρ, Te, ne, radius)
end

@recipe function plot_pellet(pelt::Pellet1)
    @series begin
        label := "pellet $(IMAS.index(pelt.properties))"
        linewidth := pelt.radius ./ pelt.radius[1] * 5
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

    return pellet

    # # plot the pellet direction in time
    # if inputs.debug_trajectory
    #     plot(dd.equilibrium; cx=true)
    #     plot!(r, z)
    #     display(scatter!([r[1]], [z[1]]))

    #     rmin, rmax = extrema(dd.equilibrium.time_slice[].boundary.outline.r)
    #     plot(x, y)
    #     t = 0:0.01:2π
    #     plot!(rmin .* cos.(t), rmin .* sin.(t))
    #     display(plot!(rmax .* cos.(t), rmax .* sin.(t); aspect_ratio=:equal))
    # end



    # # plot the profiles that the pellet sees
    # if inputs.debug_plasma
    #     display(scatter(ρ, ne; alpha=0.1))
    #     display(scatter(ρ, Te; alpha=0.1))
    # end

    # n = 100
    # r_t, z_t = r[n], z[n]

    # # #------------ start OMFIT part with rpdot for each t-----------------------------------------
    # # r_t, z_t = PAM.pellet_position_t(starting_position, velocity_vector, 0.02, tinj)
    # # rho = RHO_interpolant.(r_t, z_t)

    # # if inputs.debug_trajectory
    # #     println("-------Pellet1 position at t--------")
    # #     plot(dd.equilibrium; cx=true)
    # #     plot!(r, z)
    # #     display(scatter!([r_t], [z_t]))

    # # end


    # # initial pellet radius
    # y0 = pelt[1].shape.size[]
    # nlayers = length(pelt[1].layer)  # number of layers in the pellet
    # # generate array of pellet layers
    # layers_radii0 = []
    # for i = 1:nlayers
    #     push!(layers_radii0, pelt[1].layer[1].thickness)
    # end

    # nComp = length(pelt[1].layer[1].species)
    # layer_index = PAM.get_ilayer(y0, layers_radii0, nComp)
    # y = y0  # will define y somethere later in the function
    # # apply ablation model only if the pellet inside the plasma rho <1.0 - no plasma inside the SOL in FUSE?

    # nComp = length(pelt[1].layer[layer_index].species)  #
    # labels = []

    # Wavg = 0
    # Wratio = 1
    # for i in 1:nComp
    #     label = pelt[1].layer[layer_index].species[i].label
    #     push!(labels, label)
    #     Wavg += pelt[1].layer[layer_index].species[i].fraction * pelt[1].layer[layer_index].species[i].a

    #     @eval global $(Symbol("frac", label)) = $(pelt[1].layer[layer_index].species[i].fraction)
    #     @eval global $(Symbol("W", label)) = $(pelt[1].layer[layer_index].species[i].a)
    #     @eval global $(Symbol("Z", label)) = $(pelt[1].layer[layer_index].species[i].z_n)

    # end

    # if ("D" in labels) & ("T" in labels)
    #     println("P.P. model for deuterium-tritium ablation")
    #     model = "DT"
    # elseif ("Ne" in labels) & ("D" in labels)
    #     println("P.P. model for neon-deuterium ablation")
    #     model = "NeD"
    # elseif (nComp == 1) & ("C" in labels)
    #     println("Model for carbon ablation")
    #     model = "C"

    # elseif (nComp == 1) & ("Li" in labels)
    #     println("Model for Litium ablation")
    #     model = "Li"

    # else
    #     println("Parks model for the ablation")
    #     model = "Parks"

    # end



    # drpdt = PAM.ablation_model(y0, rho, Bt, Bt_exp, layer_index, ne, Te, model)  # all of this hould be for particular time step (iteration) if update
    # # same, but for all pellets? y0 should be changed to current y

    # # update: get shift (drift?), update pellet density
    # # variable for the pellet density? 
    # np = PAM.update_pellet_densities()




end

end
