import PAM
import IMAS
using Test

@testset "PAM" begin
    dd_D3D_json = IMAS.json2imas(joinpath(@__DIR__, "..", "examples", "template_D3D_1layer_2species.json"))
    dd_D3D = IMAS.hdf2imas(joinpath(@__DIR__, "..", "examples", "template_D3D_1layer_2species.h5"))

    @test dd_D3D_json == dd_D3D

    dd_D3D.pellets.time_slice[].pellet[1].velocity_initial = 200.0

    for drift_model in (:Parks, :HPI2, :none)
        println(drift_model)
        @test begin
            PAM.run_PAM(dd_D3D;
                t_start=0.0,
                t_finish=0.0045,
                time_step=0.0001,
                drift_model,
                Bt_dependance=true,
                update_plasma=true)
            true
        end
    end

end

@testset "Comparison with OMFIT PAM" begin
    dd_D3D = IMAS.hdf2imas(joinpath(@__DIR__, "..", "examples", "template_D3D_1layer_2species.h5"))
    dd_D3D.pellets.time_slice[].pellet[1].velocity_initial = 200.0;

    inputs=(
        t_start = 0.0,
        t_finish = 0.001,
        time_step = 0.0001, 
        drift_model= :HPI2,
        Bt_dependance=true,
        update_plasma = false,
    )

    pellet_PAM = PAM.run_PAM(dd_D3D; inputs...);

    # Reference from OMFIT PAM
    omfit_radius = 0.01 * [0.1, 0.1, 0.09907301, 0.09514626, 0.08962336, 0.08174432, 0.07028775, 0.05297672, 0.02178495, 0.0, 0.0]
    omfit_R = [2.3, 2.28, 2.26, 2.24, 2.22, 2.2, 2.18, 2.16, 2.14, 2.12, 2.1]
    omfit_ablation_D = [0.0, 0.0, 1.15829553e+23, 1.62737315e+23, 2.03252305e+23, 2.43448343e+23, 2.64645468e+23, 2.35555840e+23,
        9.11641845e+22, 0.0, 0.0]

    @test isapprox(pellet_PAM.radius, omfit_radius, rtol=5e-3)
    @test isapprox(pellet_PAM.r, omfit_R, rtol=5e-3)
    @test isapprox(pellet_PAM.ablation_rate, 2 * omfit_ablation_D, rtol=5e-3)

    # Regression test (with previous PAM version)
    @test isapprox(maximum(pellet_PAM.density_source), 8.79030213181948e19)
    @test isapprox(sum(pellet_PAM.density_source), 3.9712279520318335e21)

    @test isapprox(pellet_PAM.ablation_rate[4], 3.2596193988579004e23)
    @test isapprox(pellet_PAM.ablation_rate[9], 1.7870515483957858e23)
    @test isapprox(pellet_PAM.ne[4], 3.437934198621271e19)
    @test isapprox(pellet_PAM.Te[4],  517.7166642235971)
    @test isapprox(pellet_PAM.temp_drop[4],  0.9994708595114168)
    @test isapprox(pellet_PAM.radius[4],  0.0009515680832308055)
    @test isapprox(pellet_PAM.x[4], 2.2399999999999998)
end

@testset "I/O and Equality operators" begin
    pellet = PAM.hdf2pellet(joinpath(pkgdir(PAM), "examples", "pellet_D3D.h5"))

    pellet2 = deepcopy(pellet)

    @test pellet == pellet2
    @test isequal(pellet, pellet2)
    @test isapprox(pellet, pellet2)
    @test !diff(pellet, pellet2)

    # Modify a field slightly
    pellet2.Te[end] += 1.0e-10

    @test pellet != pellet2
    @test !isequal(pellet, pellet2)
    @test isapprox(pellet, pellet2)
    @test diff(pellet, pellet2)

    @test !isapprox(pellet, pellet2; atol=1e-12) # should fail for 1e-10 difference
    @test isapprox(pellet, pellet2; atol=1e-5) # should pass for 1e-10 difference
end
