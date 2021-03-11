using StochasticGroundMotionSimulation
using Test

@testset "StochasticGroundMotionSimulation.jl" begin

    @testset "Oscillator" begin
        ζ = 0.05
        f_n = 1.0
        sdof = Oscillator(f_n, ζ)

        @test f_n ≈ 1.0/period(sdof)

        @test transfer(0.5, sdof)^2 ≈ squared_transfer(0.5, sdof)

        fi = [ 0.5, 1.0, 2.0 ]
        Hfi = transfer(fi, sdof)
        transfer!(Hfi, 2 .* fi, sdof)
        @test Hfi ≈ transfer(2 .* fi, sdof)

        squared_transfer!(Hfi, fi, sdof)
        @test Hfi ≈ transfer(fi, sdof).^2

    end

    @testset "Site Response" begin

        f = 0.05
        Af0 = site_amplification(f)
        Af1 = site_amplification(f; amp_model=:AlAtik2021_cy14)
        Af3 = site_amplification(f, amp_model=:Unit)

        @test Af0 == Af1
        @test Af3 == 1.0

        f0 = 80.0
        f1 = 100.0
        Af0 = site_amplification(f0, amp_model=:Boore2016)
        Af1 = site_amplification(f1, amp_model=:Boore2016)

        @test Af0 == Af1

    end


end
