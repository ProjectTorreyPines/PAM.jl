module PAM

import IMAS
using Plots
import DifferentialEquations

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

mutable struct Pellet{A,T,N,S,B,X}
    properties::IMAS.pellets__time_slice___pellet
    Btdep::B
    update_plasma::B
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
    temp_drop::A
end

function Pellet(
    pelt::IMAS.pellets__time_slice___pellet,
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    time::Vector{Float64},
    surfaces::Vector{IMAS.FluxSurface},
    drift_model::Symbol,
    Bt_dependance::Bool,
    update_plasma::Bool
)

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
    velocity_vector = [X2 - X1, Y2 - Y1, Z2 - Z1] .* pelt.velocity_initial ./ (dist)

    # trajectory for all time steps
    x, y, z = pellet_position(starting_position, velocity_vector, time, time[1])
    r, z = PAM.projected_coordinate(x, y, z)

    # rho interpolant
    _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)

    # check the pellet position
    rho_start = RHO_interpolant.(R1, Z1)

    if (rho_start < 1)
        error("Pellet starting inside plasma at rho = $rho_start")
    end

    #Bt = abs(eqt.eqt.global_quantities.magnetic_axis.r .b_field_tor) .* eqt.global_quantities.magnetic_axis.r ./r
    Bt = abs(eqt.global_quantities.vacuum_toroidal_field.b0) .* eqt.global_quantities.vacuum_toroidal_field.r0 ./ r

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
    temp_drop = fill(0.0, size(time))
    R_drift = fill(0.0, size(time))
    density_source = fill(0.0, (length(time), length(surfaces)))
    return Pellet(
        pelt,
        Bt_dependance,
        update_plasma,
        drift_model,
        time,
        time[1],
        Bt,
        velocity_vector,
        r,
        R_drift,
        z,
        x,
        y,
        ρ,
        Te,
        ne,
        radius,
        ablation_rate,
        density_source,
        temp_drop
    )
end

"""
    get_ilayer(pelt::Pellet, k::Int)

Return the layer index based on current pellet radius and thikness of the pellet radii
"""
function get_ilayer(pelt::Pellet, k::Int)
    layers = cumsum([layer.thickness for layer in pelt.properties.layer])
    pos = pelt.radius[k] .- layers
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

function drift!(::Val{:none}, pelt::Pellet, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, k::Int)
    pelt.R_drift[k] = 0.0
    return
end

function drift!(::Val{:HPI2}, pelt::Pellet, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, k::Int)
    Raxis = eqt.global_quantities.magnetic_axis.r
    Zaxis = eqt.global_quantities.magnetic_axis.z

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

    vp = sqrt(pelt.velocity_vector[1]^2 + pelt.velocity_vector[2]^2 + pelt.velocity_vector[3]^2)

    rp = pelt.radius[1] * 1e3  # from m to mm
    ne0 = cp1d.electrons.density[1] * 1e-19 #to 1e19 m^-3
    te0 = cp1d.electrons.temperature[1] * 1e-3
    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    a0 = eqt.boundary.minor_radius
    kappa = eqt.boundary.elongation
    alpha = atan(pelt.z[1] - Zaxis, pelt.r[1] - Raxis)
    lam = minimum(pelt.ρ)
    
    dr_drift = real(c1 * (vp / 100)^c2 * rp^c3  * ne0^c4 * te0^c5 * (abs(abs(alpha) - c6) + c8)^c7  * (1 - lam)^c9 * a0^c10 * R0^c11 * Complex(B0)^c12 * kappa^c13)
    pelt.R_drift[k] = dr_drift

    return
end

function drift!(::Val{:Parks}, pelt::Pellet, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, k::Int)
    Raxis = eqt.global_quantities.magnetic_axis.r
    Zaxis = eqt.global_quantities.magnetic_axis.z
    _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)

    # ---- constants ----------------
    m_p = 1.67262192369e-27
    e = 1.602176634e-19
    mu_0 = 1.25663706212e-6

    mi = m_p * 2
    T0 = 2.0    # eV, temperature at the pellet surface (channel entrance)
    M0 = 1  # Mach number at the channel entrance, from 2007 Roman's paper, in 2000 paper M0=0.8
    c_perp = 0.02 # pellet cloud width factor taken as a constant from OMFIT
    #---------------------------------------------------------
    Ca_inf = pelt.Bt[k] / sqrt(mu_0 * pelt.ne[k] * mi)
    C_s = sqrt(e * T0 / mi)

    if pelt.Btdep
        attenuation_factor = pelt.Bt[k]^0.7
    else
        attenuation_factor = 2.0^0.7
    end

    rperp_calc = c_perp * sqrt(pelt.ablation_rate[k] / (pelt.ne[k] * pelt.Te[k]^1.5 * attenuation_factor))

    rperp = max(pelt.radius[1], rperp_calc)

    Lc0 = sqrt(Raxis * rperp)

    pfact_exp = 1.69423774
    Sigma0_exp = 0.85

    loglam = 23.5 - log((pelt.ne[k] * 1e-6)^0.5 / pelt.Te[k]^(5 / 6))
    tau_inf = 2.248882e16 * pelt.Te[k]^2 / loglam
    a1 = 0.01 * M0

    n0_0 = pelt.ablation_rate[k] / (2 * π * rperp^2 * C_s) / (a1 * (T0 / pelt.ne[k] / pelt.Te[k])^(pfact_exp) * (Lc0 / tau_inf)^(Sigma0_exp))

    n0 = n0_0^(1 / (1 + pfact_exp + Sigma0_exp))

    Sigma0 = n0 * Lc0 / tau_inf
    pe = pelt.Te[k] * pelt.ne[k] * e

    pe_rho = e .* cp1d.electrons.density .* cp1d.electrons.temperature

    function vRdot!(du, u, p, t)

        pfact_exp = 1.69423774
        Sigma0_exp = 0.85

        R0 = u[1]

        vR0 = u[2]

        Z0 = pelt.z[k]

        tbar = t / Lc0 * C_s

        rho_cloud = RHO_interpolant.(R0, Z0)[1]

        if rho_cloud > 2.0
            du[1] = 0.0
            du[2] = 0.0
        end

        peinf = IMAS.interp1d(cp1d.grid.rho_tor_norm, pe_rho).(rho_cloud)
        peinf = min(e * n0 * T0, max(pe, peinf))

        Pfact = e * n0 * T0 / peinf

        a1 = 0.52296257 / Lc0 / (Sigma0^Sigma0_exp * Pfact^pfact_exp)
        a2 = 17.19090852
        Lc = 1.0 + (Sigma0^Sigma0_exp * Pfact^pfact_exp) * tbar

        PSI = a1 * (exp(-a2 * (Lc - Lc0) / Lc0 / Pfact^(2 * pfact_exp) / Sigma0^(2 * Sigma0_exp)) - 1 / Pfact) * Lc / Lc0 + (1 - a1) * (1 - 1 / Pfact)
        PSI = max(0.0, PSI)

        if PSI <= 0.0
            du[1] = 0.0
            du[2] = 0.0
        else
            du[1] = vR0
            du[2] = -2 * pelt.Bt[k]^2 / Ca_inf / mu_0 * vR0 / (mi * n0 * Lc0) + 2 / R0 * PSI * (C_s^2 / Lc0)
        end
    end

    t_eval = pelt.time[k:end] .- pelt.time[k]

    t_span = [0, (t_eval[end] - t_eval[1])]

    u0 = [pelt.r[k], pelt.velocity_vector[1]]

    prob = DifferentialEquations.ODEProblem(vRdot!, u0, t_span)


    #----------------------- Comment ---------------------------------
    #  RadauIIA5(; autodiff=false) is the most accurate method for this problem and provide nice solution, but fails during the regression
    #  test as can't work without autodiff option
    #  Tsit5() or Vern7() are able to solve the system, but the solution is oscillating, but use it as it should work for the regression test 
    sol = DifferentialEquations.solve(prob, DifferentialEquations.RadauIIA5(; autodiff=false); verbose=false)
    #sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(); verbose=false)
    #------------------------------------------------------------------------------------------------------------
    #sol = DifferentialEquations.solve(prob, DifferentialEquations.Vern7(); verbose=false)

    dr_drift = sol[1, end] - pelt.r[k]

    pelt.R_drift[k] = dr_drift

    return
end

function dr_dt!(pelt::Pellet, k::Int)
    ilayer = get_ilayer(pelt, k)
    layer = pelt.properties.layer[ilayer]

    Te_func = IMAS.interp1d(pelt.time, pelt.Te)
    ne_func = IMAS.interp1d(pelt.time, pelt.ne)
    Bt_function = IMAS.interp1d(pelt.time, pelt.Bt)

    A = Float64[]
    Z = Float64[]
    fractions = Float64[]
    species_list = String[]

    for species in layer.species
        tmp = IMAS.ion_properties(Symbol(species.label))
        push!(A, tmp.a)
        push!(Z, tmp.z_n)
        push!(fractions, species.fraction)
        push!(species_list, species.label)

    end

    A_mean = sum(fractions .* A)
    @assert sum(fractions) == 1.0 "Species fractions don't sum up to 1.0"
    if pelt.Btdep
        Bt_exp = -0.35  #scaling of the magnetic field 
    else
        Bt_exp = 0
    end

    # ablation model for the DT pellet
    if ("D" in species_list) & ("T" in species_list)
        if species_list[1] == "D"
            AD = A[1]
            AT = A[2]
            FD = fractions[1]
            FT = fractions[2]
            ZD = Z[1]
            ZT = Z[2]

        else
            AD = A[2]
            AT = A[1]
            FD = fractions[2]
            FT = fractions[1]
            ZD = Z[2]
            ZT = Z[1]
        end
        # according equation #1 in J.McClenaghan at al Nucl.Fusion 63 (2023) 036015 (two coefficients are missing) and normalization to FUSE units

        ρ_zero = (1 - FD + FD * AD / AT) * ((1 - FD) / pellet_mass_density("T") + (FD * AD / AT) / pellet_mass_density("D"))^(-1) #[g cm^-3]
        Wratio = (1 - FD) * AT / AD + FD

        function dr_dt_DT(du, u, p, t)
            c0 = 8.358 * Wratio^0.6667 * (abs(Bt_function(t)) / 2.0)^Bt_exp
            if u[1] < 0
                u[1] = 0
            else
                du[1] = (-c0 / ρ_zero * (Te_func(t) * 1e-3)^(1.6667) * (ne_func(t) * 1e-20)^(0.3333) / (u[1] * 1e2)^0.6667) * 1e-2
            end
        end

        t_span = [pelt.time[k-1], pelt.time[k]]

        u0 = [pelt.radius[k-1]]

        prob = DifferentialEquations.ODEProblem(dr_dt_DT, u0, t_span)

        sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(); verbose=false)

        pelt.radius[k] = sol[1, end]

        #------- ablation rate calculations ---------------------------------
        c0 = 8.358 * Wratio^0.6667 * (abs(pelt.Bt[k]) / 2.0)^Bt_exp
        drdt = -c0 / ρ_zero * (pelt.Te[k] * 1e-3)^(1.6667) * (pelt.ne[k] * 1e-20)^(0.3333) / (pelt.radius[k] * 1e2)^0.6667
        G = -drdt * (4 * π * ρ_zero * (pelt.radius[k] * 1e2)^2)

        GD = G * 6.022e23 * FD / A_mean * ZD
        GT = G * 6.022e23 * FT / A_mean * ZT

        pelt.ablation_rate[k] = GT + GD

    elseif ("Ne20" in species_list) & ("D" in species_list)

        if species_list[1] == "Ne"
            ANe = A[1]
            AD = A[2]
            FNe = fractions[1]
            FD = fractions[2]
            ZNe = Z[1]
            ZD = Z[2]

        else
            ANe = A[2]
            AD = A[1]
            FD = fractions[1]
            FNe = fractions[2]
            ZNe = Z[2]
            ZD = Z[1]
        end

        Wratio = (1 - FD) * ANe / AD + FD

        ρ_zero = (1 - FD + FD * AD / ANe) * ((1 - FD) / pellet_mass_density("Ne") + (FD * AD / ANe) / pellet_mass_density("D"))^(-1) #[g cm^-3]

        X = FD / (2 - FD)
        AoX = 27.0 + tan(1.48 * X)

        function dr_dt_DNe(du, u, p, t)
            if u[1] < 0
                u[1] = 0
            else
                c0 = AoX / (4 * π) * (abs(Bt_function(t)) / 2.0)^Bt_exp

                du[1] = (-c0 / ρ_zero * (Te_func(t) * 1e-3)^(1.6667) * (ne_func(t) * 1e-20)^(0.3333) / (u[1] * 1e2)^0.6667) * 1e-2
            end
        end

        t_span = [pelt.time[k-1], pelt.time[k]]

        u0 = [pelt.radius[k-1]]

        prob = DifferentialEquations.ODEProblem(dr_dt_DNe, u0, t_span)

        sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(); verbose=false)

        pelt.radius[k] = sol[1, end]

        #---- ablation rate calculation ----------------------------------------------

        dr_dt = -c0 / ρ_zero * (pelt.Te[k] * 1e-3)^(1.6667) * (pelt.ne[k] * 1e-20)^(0.3333) / (pelt.radius[k] * 1e2)^0.6667
        G = -dr_dt * (4 * π * ρ_zero * (pelt.radius[k-1] * 1e2)^2)

        GD = G * 6.022e23 * FD / A_mean * ZD
        GNe = G * 6.022e23 * FNe / A_mean * ZNe

        pelt.ablation_rate[k] = GD + GNe

    elseif "C12" in species_list

        C0 = 8.146777e-9
        AC = A[1]
        ZC = Z[1]
        gamma = 5.0 / 3.0

        ZstarPlus1C = 2.86
        Albedo = 23.920538030089528 * log(1 + 0.20137080524063228 * ZstarPlus1C)
        flelectro = exp(-1.936)
        fL = (1.0 - Albedo / 100) * flelectro

        IstC = 60

        xiexp = 0.601
        lamdaa = 0.0933979540623963
        lamdab = -0.7127242270013098
        lamdac = -0.2437544205933372
        lamdad = -0.8534855445478313

        function dr_dt_C(du, u, p, t)
            if u[1] < 0
                u[1] = 0
            else
                if Te_func(t) > 30
                    Ttmp = Te_func(t)
                else
                    Ttmp = 30
                end

                loglamCSlow = log(2.0 * Ttmp / IstC * sqrt(2.718 * 2.0))

                BLamdaq = 1 / (ZC * loglamCSlow) * (4 / (2.5 + 2.2 * sqrt(ZstarPlus1C)))

                av = 10.420403555938629 * (Ttmp / 2000.0)^lamdaa
                bv = 0.6879779829877795 * (Ttmp / 2000.0)^lamdab
                cv = 1.5870910225610804 * (Ttmp / 2000.0)^lamdac
                dv = 2.9695640286641840 * (Ttmp / 2000.0)^lamdad
                fugCG = 0.777686
                Gpr =
                    C0 * AC^(2.0 / 3.0) * (gamma - 1)^(1.0 / 3.0) * (fL * ne_func(t) * 1e-6)^(1.0 / 3.0) * (u[1] * 1e2)^(4.0 / 3.0) * (Te_func(t))^(11.0 / 6.0) *
                    BLamdaq^(2.0 / 3.0)

                CG =
                    fugCG * av * log(1 + bv * (ne_func(t) * 1e-20)^(2.0 / 3.0) * (u[1] * 1e2)^(2.0 / 3.0)) /
                    log(cv + dv * (ne_func(t) * 1e-20)^(2.0 / 3.0) * (u[1] * 1e2)^(2.0 / 3.0))

                G = xiexp * CG * Gpr * (2.0 / Bt_function(t))^Bt_exp

                du[1] = (-G / (4.0 * π * pellet_mass_density("C") * (u[1] * 1e2)^2)) * 1e-2
            end
        end

        t_span = [pelt.time[k-1], pelt.time[k]]

        u0 = [pelt.radius[k-1]]

        prob = DifferentialEquations.ODEProblem(dr_dt_C, u0, t_span)

        sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(); verbose=false)

        pelt.radius[k] = sol[1, end]

        #---ablation calculations --------------------------
        if pelt.Te[k] > 30
            Ttmp = pelt.Te[k]
        else
            Ttmp = 30
        end

        loglamCSlow = log(2.0 * Ttmp / IstC * sqrt(2.718 * 2.0))

        BLamdaq = 1 / (ZC * loglamCSlow) * (4 / (2.5 + 2.2 * sqrt(ZstarPlus1C)))

        av = 10.420403555938629 * (Ttmp / 2000.0)^lamdaa
        bv = 0.6879779829877795 * (Ttmp / 2000.0)^lamdab
        cv = 1.5870910225610804 * (Ttmp / 2000.0)^lamdac
        dv = 2.9695640286641840 * (Ttmp / 2000.0)^lamdad
        fugCG = 0.777686

        Gpr =
            C0 * AC^(2.0 / 3.0) * (gamma - 1)^(1.0 / 3.0) * (fL * pelt.ne[k] * 1e-6)^(1.0 / 3.0) * (pelt.radius[k] * 1e2)^(4.0 / 3.0) * (pelt.Te[k])^(11.0 / 6.0) *
            BLamdaq^(2.0 / 3.0)

        CG =
            fugCG * av * log(1 + bv * (pelt.ne[k] * 1e-20)^(2.0 / 3.0) * (pelt.radius[k] * 1e2)^(2.0 / 3.0)) /
            log(cv + dv * (pelt.ne[k] * 1e-20)^(2.0 / 3.0) * (pelt.radius[k] * 1e2)^(2.0 / 3.0))

        G = xiexp * CG * Gpr * (2.0 / pelt.Bt[k])^Bt_exp

        pelt.ablation_rate[k] = G * 6.022e23 / AC * ZC

    else
        error("No ablation model is implemented for such combination of species")
    end

    return
end

function pellet_density(pelt::Pellet, surface::IMAS.FluxSurface, k::Int)
    # calculations of the pellet cloud size
    cloudFactor = 5

    cloudFactorR = cloudFactor
    cloudFactorZ = cloudFactor

    rcloudR = pelt.radius[k] * cloudFactorR
    rcloudZ = pelt.radius[k] * cloudFactorZ

    # Assume density source as a 2D Gaussian function
    nsource = exp.(-0.5 .* ((pelt.r[k] .- surface.r .+ 0.5 .* pelt.R_drift[k]) .^ 2 ./ (rcloudR .+ 0.25 .* pelt.R_drift[k]) .^ 2 .+ (pelt.z[k] .- surface.z) .^ 2 ./ rcloudZ .^ 2))

    # need to normilize on the surface area under the gaussian shape
    nsource ./= (2 * π)^2 * (rcloudR + 0.25 * pelt.R_drift[k]) * (pelt.r[k] + 0.5 * pelt.R_drift[k]) * rcloudZ

    nsource .*= pelt.ablation_rate[k]

    return IMAS.flux_surface_avg(nsource, surface)
end

function ablate!(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, pelt::Pellet, surfaces::Vector{IMAS.FluxSurface})
    pellet_source = zeros(length(pelt.time), length(surfaces))

    for k in 2:length(pelt.time)
        dt = pelt.time[k] - pelt.time[k-1]

        if pelt.ρ[k] > 1.0 && pelt.radius[k-1] > 0
            pelt.radius[k] = pelt.radius[k-1]

        else
            dr_dt!(pelt, k)
            if isnan(pelt.ablation_rate[k])

                pelt.ablation_rate[k] = 0.0
            end

            drift!(Val(pelt.drift_model), pelt, eqt, cp1d, k)

            for (ks, surface) in enumerate(surfaces)
                tmp = pellet_density(pelt, surface, k)
                if isnan(tmp)
                    pellet_source[k, ks] += 0.0
                else
                    pellet_source[k, ks] += tmp * dt
                end
            end

            #------calculate energy sink and update plasma -----------
            # pressure remains the same, causing a drop in T
            ne0 = pelt.ne[k]

            ne_new = pelt.ne[k] + sum(pellet_source[k, :]) * (pelt.time[k] - pelt.time[k-1])

            pelt.temp_drop[k] = ne0 / ne_new

            if pelt.update_plasma
                pelt.ne[k] = ne_new
                pelt.Te[k] *= pelt.temp_drop[k]
                if pelt.Te[k] < 1
                    pelt.Te[k] = 1
                end
            end

        end
        pelt.density_source = pellet_source
    end

    return
end

@recipe function plot_pellet_deposition(pelt::Pellet; plot_cloud=true)
    # deposition = abs.(IMAS.gradient(pelt.radius))

    deposition = sum(pelt.density_source; dims=2)

    @series begin
        label := "normilized desity source Pellet $(IMAS.index(pelt.properties))"
        linewidth := deposition ./ maximum(deposition) .* 5 .+ 0.01
        linealpha := 0.5
        pelt.r, pelt.z
    end
    if plot_cloud
        @series begin
            label := "pellet cloud trajectory Pellet $(IMAS.index(pelt.properties))"
            linewidth := 2
            arrow := true
            pelt.r + pelt.R_drift, pelt.z
        end
    end

    @series begin
        primary := false
        seriestype := :scatter
        markersize := pelt.radius[1] * 1000
        [pelt.r[1]], [pelt.z[1]]
    end
end

"""
    run_PAM(dd::IMAS.dd; t_start::Float64, t_finish::Float64, time_step::Float64, drift_model::Symbol, Bt_dependance::Bool, update_plasma::Bool)

Run Pellet Ablation Model and write data to 
"""
function run_PAM(dd::IMAS.dd; t_start::Float64, t_finish::Float64, time_step::Float64, drift_model::Symbol, Bt_dependance::Bool, update_plasma::Bool)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    # generate time array for the simulations
    time = collect(range(t_start, t_finish; step=time_step))

    # define flux surfaces, will be needed for the pellet source calculations
    surfaces = IMAS.trace_surfaces(eqt, IMAS.first_wall(dd.wall)...)

    # initialize the pellet structure 
    pellet = Pellet(dd.pellets.time_slice[].pellet[1], eqt, cp1d, time, surfaces, drift_model, Bt_dependance, update_plasma)

    ablate!(eqt, cp1d, pellet, surfaces)

    return pellet
end

export run_PAM

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
