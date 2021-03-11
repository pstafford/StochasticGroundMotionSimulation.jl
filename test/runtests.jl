using StochasticGroundMotionSimulation
using Test
using BenchmarkTools

@testset "StochasticGroundMotionSimulation.jl" begin

    @testset "Performance" begin
        Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]
        m = 6.0
        r = 10.0
        fas = FASParams(100.0, [1.0, 50.0, Inf], [1.0, 0.5], 200.0, 0.4, 0.039 )

        @time Sai = rvt_response_spectrum(Ti, m, r, fas)
        #@btime Sai = rvt_response_spectrum( Ti, m, r, fas )
        @time rvt_response_spectrum_cy!(Sai, Ti, m, r, fas)

    end

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
