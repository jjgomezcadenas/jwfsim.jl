using Test
using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
#using Interpolations

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    fs, ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W, kW
using jwfsim


function laser_test()
    epl375l = jwfsim.epl375(1.0MHz)
    pq375l  = jwfsim.pqLDH375B(1.0MHz)
    pr = uconvert(NoUnits, pq375l.P/epl375l.P)
    @test pr > 12

    lsr = epl375l
    obj = jwfsim.Objective(0.5, 40.0, 5.0mm, 5.1mm, 7.8mm, 0.85)
    w0 = jwfsim.w0f(lsr.λ, obj.d, obj.f)
    glsr = jwfsim.GaussianLaser(lsr, w0)
    Igl = jwfsim.I(glsr)
    @test uconvert(kW/cm^2, Igl(0.0nm,0.0nm)) ≈ uconvert(kW/cm^2, 2*lsr.P/(π*w0)^2)
end

# function lmu_40x_nuv()
#     NA = 0.5
#     M  = 40
#     f  = 5.0mm
#     d  = 5.1mm
#     wd = 7.8mm
#     T  = 0.85
#     jwfsim.Objective(NA, M, f, d, wd, T)
# end

# @testset "glaser_test.jl" begin
    
#     epl375 = LaserLab.PulsedLaser(375.0nm, 140.0mW, 1.0MHz, 65.0ps)
#     lmu_40x_nuv = LaserLab.Objective("LMU-40X-NUV", 5.3mm, 5.1mm, 40.0)
#     w0 = lmu_40x_nuv.d /2.0
#     gepl375 = LaserLab.GaussianLaser(epl375, w0)
#     g2epl375 = LaserLab.propagate_paralell_beam(epl375, lmu_40x_nuv)
    
    
#     @test gepl375.w0 ≈ w0
#     @test gepl375.laser.P ≈ epl375.P
#     @test g2epl375.w0 ≈ epl375.λ/(π * lmu_40x_nuv.NA)
#     @test LaserLab.spot_size(g2epl375) ≈ 2 * g2epl375.w0
#     @test LaserLab.angular_divergence(g2epl375) ≈ 2 * g2epl375.θ0
#     @test LaserLab.depth_of_focus(g2epl375)  ≈ 2 * g2epl375.z0

#     fI = LaserLab.I(g2epl375) 
#     I0 = fI(0.0*μm, 0.0*μm )
#     @test g2epl375.I0 ≈ I0

#     aspot = π * (LaserLab.spot_size(g2epl375) / 2.0)^2 
#     @test aspot ≈ π * g2epl375.w0^2

#     P0 = I0 * aspot
#     @test LaserLab.n_photons(g2epl375.laser.λ, P0) ≈ LaserLab.n_photons(g2epl375.laser.λ, g2epl375.I0 * aspot)

#     nγ = uconvert(Hz, epl375.P / LaserLab.photon_energy(epl375.λ))
#     @test 2 * nγ / (π * gepl375.w0^2) == gepl375.γ0
# end
