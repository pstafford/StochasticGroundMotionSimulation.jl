using StochasticGroundMotionSimulation
using Test
using BenchmarkTools

@testset "StochasticGroundMotionSimulation.jl" begin

    @testset "Performance" begin
        Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]
        m = 4.0+π
        r = 500.0+π
        fas = FASParams(100.0, [1.0, 50.0, Inf], [1.0, 0.5], 200.0, 0.4, 0.039 )

        @time Sai = rvt_response_spectrum(Ti, m, r, fas)
        # @btime Sai = rvt_response_spectrum(Ti, m, r, fas)
        # @btime rvt_response_spectrum_cy!(Sai, Ti, m, r, fas)
    
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

    @testset "Duration" begin

        # Boore & Thompson 2014
        fas = FASParams(100.0, [1.0, 50.0, Inf], [1.0, 0.5], 200.0, 0.4, 0.039)

        m = 6.0
        fa, fb, ε = corner_frequency(m, fas)
        Ds = 1.0 / fa
        @test boore_thompson_2014(m, 0.0, fas) ≈ Ds

        m = 6.0
        r = 7.0
        fa, fb, ε = corner_frequency(m, fas)
        Dur = 1.0 / fa + 2.4
        @test boore_thompson_2014(m, r, fas) ≈ Dur

        fa, fb, ε = corner_frequency(m, fas, fc_fun=:Atkinson_Silva_2000)
        Ds = 0.5 * ( 1.0 / fa + 1.0 / fb )
        @test boore_thompson_2014(m, 0.0, fas, fc_fun=:Atkinson_Silva_2000) ≈ Ds

        c11 = [ 8.4312e-01, -2.8671e-02, 2.0,  1.7316e+00,  1.1695e+00,  2.1671e+00,  9.6224e-01 ]
        c11f = boore_thompson_2012_coefs(1, 1)

        @test all(isapprox.(c11, c11f))

        sdof = Oscillator(1.0)

        Drms, Dex, Dratio = boore_thompson_2012(6.0, 10.0, fas, sdof)
        Dex0 = boore_thompson_2014(6.0, 10.0, fas)

        @test Dex == Dex0




    end

end
