import PAM
import IMAS
using Test

@testset "PAM" begin
    dd_D3D = IMAS.json2imas(joinpath(@__DIR__, "..", "examples", "template_D3D_1layer_2species.json"))

    dd_D3D.pellets.time_slice[].pellet[1].velocity_initial = 200.0

    @test begin
        try
            PAM.run_PAM(dd_D3D;
                t_start=0.0,
                t_finish=0.0045,
                time_step=0.0001,
                drift_model=:Parks,
                Bt_dependance=true,
                update_plasma=true)
            return true
        catch
            return false
        end
    end

end
