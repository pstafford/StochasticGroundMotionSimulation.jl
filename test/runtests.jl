using StochasticGroundMotionSimulation
using Test
using Traceur
using ForwardDiff: Dual
using FastGaussQuadrature
using QuadGK
using LinearAlgebra


@testset "StochasticGroundMotionSimulation.jl" begin

    @testset "Performance" begin
        Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]
        m = 4.0+π
        r = 500.0+π

        src = SourceParameters(100.0)
        geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
        ane = AnelasticAttenuationParameters(200.0, 0.4)
        sat = NearSourceSaturationParameters(:BT15)
        path = PathParameters(geo, sat, ane)
        site = SiteParameters(0.039)

        fas = FourierParameters(src, path, site)
        rvt = RandomVibrationParameters(:CL56)


        sdof = Oscillator(1.0)
        # @trace Oscillator(1.0)
        @code_warntype Oscillator(1.0)
        @code_warntype period(sdof)
        @code_warntype transfer(1.0, sdof)

        # @trace rvt_response_spectral_ordinate(Ti[1], m, r, fas, rvt)
        @code_warntype rvt_response_spectral_ordinate(Ti[1], m, r, fas, rvt)
        @code_warntype rvt_response_spectrum(Ti, m, r, fas, rvt)
        @time Sai = rvt_response_spectrum(Ti, m, r, fas, rvt)
        # @btime Sai = rvt_response_spectrum(Ti, m, r, fas)
        # @btime rvt_response_spectrum_cy!(Sai, Ti, m, r, fas)

    end

    @testset "Source" begin

        m = 6.0
        Δσ = 100.0
        β = 3.5

        # @trace magnitude_to_moment(m)
        @code_warntype magnitude_to_moment(m)

        @code_warntype corner_frequency_brune(m, Δσ)
        @code_warntype corner_frequency_brune(m, Δσ, β)
        @code_warntype corner_frequency_atkinson_silva_2000(m)

        srcf = SourceParameters(Δσ)
        srcd = SourceParameters(Dual{Float64}(Δσ))

        @code_warntype corner_frequency(m, srcf)
        @code_warntype corner_frequency(m, srcd)
        @code_warntype corner_frequency(Dual(m), srcf)
        @code_warntype corner_frequency(Dual(m), srcd)

        # get_parametric_type(srcd)

        # @trace corner_frequency(m, srcf)
        # @trace corner_frequency(m, srcd)

        @time faf, fbf, fεf = corner_frequency(m, srcf)
        @time fad, fbd, fεd = corner_frequency(m, srcd)

        @test faf == fad.value

        srcf = SourceParameters(Δσ, :Atkinson_Silva_2000)
        srcd = SourceParameters(Dual{Float64}(Δσ), :Atkinson_Silva_2000)

        @time faf, fbf, fεf = corner_frequency(m, srcf)
        @time fad, fbd, fεd = corner_frequency(m, srcd)

        @test faf == fad.value
        @test fbf == fbd.value
        @test fεf == fεd.value

    end

    @testset "Path" begin

        Rrefi = [ 1.0, 50.0, Inf ]
        γi = [ 1.0, 0.5 ]
        geof = GeometricSpreadingParameters(Rrefi, γi)
        geod = GeometricSpreadingParameters(Rrefi, [ 0.5 ], [ Dual{Float64}(1.0) ], BitVector([1,0]), :Piecewise)

        sat = NearSourceSaturationParameters(:BT15)

        Q0 = 200.0
        anef = AnelasticAttenuationParameters(Q0)
        aned = AnelasticAttenuationParameters(Dual{Float64}(Q0))

        pathf = PathParameters(geof, sat, anef)
        pathd = PathParameters(geod, sat, aned)

        r = 10.0
        m = 6.0
        @code_warntype near_source_saturation(m, pathf.saturation)
        @code_warntype near_source_saturation(m, pathd.saturation)

        @code_warntype near_source_saturation(m, pathf)
        @code_warntype near_source_saturation(m, pathd)

        @code_warntype equivalent_point_source_distance(r, m, pathf)
        @code_warntype equivalent_point_source_distance(r, m, pathd)

        @code_warntype geometric_spreading_piecewise(r, geof)
        @code_warntype geometric_spreading_piecewise(r, geod)

        @code_warntype geometric_spreading_cy14(r, geof)
        @code_warntype geometric_spreading_cy14(r, geod)

        @code_warntype geometric_spreading(r, geof)
        @code_warntype geometric_spreading(r, geod)

        @time geometric_spreading(r, m, geof, sat)
        @time geometric_spreading(r, m, geod, sat)

        f = 1.0
        r = 100.0
        @code_warntype anelastic_attenuation(f, r, anef)
        @code_warntype anelastic_attenuation(f, r, aned)

        @time anelastic_attenuation(f, r, anef)
        @time anelastic_attenuation(f, r, aned)

    end

    @testset "Site" begin

        κ0f = 0.039
        κ0d = Dual{Float64}(κ0f)

        site0f = SiteParameters(κ0f)
        siteAf = SiteParameters(κ0f, :AlAtik2021_cy14)
        siteBf = SiteParameters(κ0f, :Boore2016)
        siteUf = SiteParameters(κ0f, :Unit)
        siteNf = SiteParameters(κ0f, :NaN)

        site0d = SiteParameters(κ0d)
        siteAd = SiteParameters(κ0d, :AlAtik2021_cy14)
        siteBd = SiteParameters(κ0d, :Boore2016)
        siteUd = SiteParameters(κ0d, :Unit)
        siteNd = SiteParameters(κ0d, :NaN)

        f = 0.05
        @code_warntype site_amplification(f, site0f)
        @code_warntype site_amplification(f, site0d)

        @time Af0f = site_amplification(f, site0f)
        @time Af1f = site_amplification(f, siteAf)
        @time Af2f = site_amplification(f, siteBf)
        @time Af3f = site_amplification(f, siteUf)
        @time Af4f = site_amplification(f, siteNf)

        @test Af0f == Af1f
        @test Af3f == 1.0
        @test isnan(Af4f)

        @time Af0d = site_amplification(f, site0d)
        @time Af1d = site_amplification(f, siteAd)
        @time Af2d = site_amplification(f, siteBd)
        @time Af3d = site_amplification(f, siteUd)
        @time Af4d = site_amplification(f, siteNd)

        @test Af0d == Af1d
        @test Af3d == 1.0
        @test isnan(Af4d)

        @test Af0f == Af0d

        f0 = 80.0
        f1 = 100.0
        Af0 = site_amplification(f0, siteBf)
        Af1 = site_amplification(f1, siteBf)

        @test Af0 == Af1

        f = 10.0
        @code_warntype kappa_filter(f, siteAf)
        @code_warntype kappa_filter(f, siteAd)
        @time kappa_filter(f, siteAf)
        @time kappa_filter(f, siteAd)

    end

    @testset "Oscillator" begin
        ζ = 0.05
        f_n = 1.0
        sdof = Oscillator(f_n, ζ)

        @test f_n ≈ 1.0/period(sdof)

        @test transfer(0.5, sdof)^2 ≈ squared_transfer(0.5, sdof)

        fi = [ 0.5, 1.0, 2.0 ]

        @code_warntype transfer(fi, sdof)
        # @trace transfer(fi, sdof)

        @time Hfi = transfer(fi, sdof)
        tfi = 2 * fi
        @time transfer!(Hfi, tfi, sdof)
        @test Hfi ≈ transfer(tfi, sdof)

        squared_transfer!(Hfi, fi, sdof)
        @test Hfi ≈ transfer(fi, sdof).^2

    end



    @testset "Duration" begin

        src = SourceParameters(100.0)
        geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
        ane = AnelasticAttenuationParameters(200.0, 0.4)
        sat = NearSourceSaturationParameters(:BT15)
        path = PathParameters(geo, sat, ane)
        site = SiteParameters(0.039)

        fas = FourierParameters(src, path, site)

        # Boore & Thompson 2014

        m = 6.0
        fa, fb, ε = corner_frequency(m, src)
        Ds = 1.0 / fa

        Δσf = 100.0
        Δσd = Dual{Float64}(Δσf)
        srcf = SourceParameters(Δσf)
        srcd = SourceParameters(Δσd)
        fasf = FourierParameters(srcf, path, site)
        fasd = FourierParameters(srcd, path, site)

        fa, fb, ε = corner_frequency(m, srcf)
        @code_warntype boore_thompson_2014(m, 0.0, srcf)
        fa, fb, ε = corner_frequency(m, srcd)
        @code_warntype boore_thompson_2014(m, 0.0, srcd)


        @code_warntype boore_thompson_2014(m, 0.0, fasf)
        @test boore_thompson_2014(m, 0.0, fasf) ≈ Ds
        @test boore_thompson_2014(m, 0.0, fasd) ≈ Ds

        m = 6.0
        r = 7.0
        fa, fb, ε = corner_frequency(m, fasf)
        Dur = 1.0 / fa + 2.4
        @test boore_thompson_2014(m, r, fasf) ≈ Dur
        fa, fb, ε = corner_frequency(m, fasd)
        Dur = 1.0 / fa + 2.4
        @test boore_thompson_2014(m, r, fasd) ≈ Dur


        srcAS = SourceParameters(Δσf, :Atkinson_Silva_2000)
        fa, fb, ε = corner_frequency(m, srcAS)
        Ds = 0.5 * ( 1.0 / fa + 1.0 / fb )
        @test boore_thompson_2014(m, 0.0, srcAS) ≈ Ds


        h = 0.1
        m1 = 8.0
        m2 = m1 + h
        r_ps1 = 1.0 + near_source_saturation(m1, fasf)
        r_ps2 = 1.0 + near_source_saturation(m2, fasf)

        Dex1 = boore_thompson_2014(m1, r_ps1, fasf)
        Dex2 = boore_thompson_2014(m2, r_ps2, fasf)

        log(Dex2/Dex1)/h

        using ForwardDiff
        d(x) = log(boore_thompson_2014(x[1], 1.0 + near_source_saturation(x[1], fasf), fasf))
        gd(x) = ForwardDiff.gradient(d, x)
        gd([8.0])[1]


        rvt = RandomVibrationParameters(:BT14)
        @code_warntype excitation_duration(m, r, fasf, rvt)
        @code_warntype excitation_duration(m, r, fasd, rvt)


        c11 = [ 8.4312e-01, -2.8671e-02, 2.0,  1.7316e+00,  1.1695e+00,  2.1671e+00,  9.6224e-01 ]
        c11f = boore_thompson_2012_coefs(1, 1)

        @code_warntype boore_thompson_2012_coefs(1, 1)

        @time c = boore_thompson_2012_coefs(1, 1)

        m = 8.0
        r = 1.0
        Dex = boore_thompson_2014(m, r, srcf)
        # get the oscillator period
        sdof = Oscillator(100.0)
        T_n = period(sdof)
        ζ = sdof.ζ_n
        # define the η parameter as T_n/Dex
        η = T_n / Dex

        @time boore_thompson_2012_base(η, c, ζ)
        @code_warntype boore_thompson_2012_base(η, c, ζ)

        Dex = boore_thompson_2014(m, r, srcd)
        # get the oscillator period
        sdof = Oscillator(1.0)
        T_n = period(sdof)
        ζ = sdof.ζ_n
        # define the η parameter as T_n/Dex
        η = T_n / Dex

        @time boore_thompson_2012_base(η, c, ζ)
        @code_warntype boore_thompson_2012_base(η, c, ζ)

        @time boore_thompson_2012(m, r, srcf, sdof, rvt)
        @code_warntype boore_thompson_2012(m, r, srcf, sdof, rvt)
        # @trace boore_thompson_2012(m, r, srcf, sdof, rvt)

        @time boore_thompson_2012(m, r, srcd, sdof, rvt)
        @code_warntype boore_thompson_2012(m, r, srcd, sdof, rvt)


        @test all(isapprox.(c11, c11f))

        sdof = Oscillator(1.0)

        Drms, Dex, Dratio = boore_thompson_2012(6.0, 10.0, fas, sdof, rvt)
        Dex0 = boore_thompson_2014(6.0, 10.0, fas)

        @test Dex == Dex0

        @code_warntype rms_duration(m, r, srcf, path, sdof, rvt)
        @code_warntype rms_duration(m, r, srcd, path, sdof, rvt)

        @code_warntype rms_duration(m, r, fasf, sdof, rvt)
        @code_warntype rms_duration(m, r, fasd, sdof, rvt)

        @time Drmsf, Dexf, Dratiof = rms_duration(m, r, fasf, sdof, rvt)
        @time Drmsd, Dexd, Dratiod = rms_duration(m, r, fasd, sdof, rvt)

        @test Drmsf == Drmsd.value
        @test Dexf == Dexd.value
        @test Dratiof == Dratiod.value

    end

    @testset "Fourier" begin

        Δσf = 100.0
        γ1f = 1.0
        γ2f = 0.5
        Q0f = 200.0
        ηf = 0.4
        κ0f = 0.039

        Δσd = Dual{Float64}(Δσf)
        γ1d = Dual{Float64}(γ1f)
        γ2d = Dual{Float64}(γ2f)
        Q0d = Dual{Float64}(Q0f)
        ηd = Dual{Float64}(ηf)
        κ0d = Dual{Float64}(κ0f)

        srcf = SourceParameters(Δσf)
        srcd = SourceParameters(Δσd)

        Rrefi = [1.0, 50.0, Inf]
        geof = GeometricSpreadingParameters(Rrefi, [γ1f, γ2f])
        geod = GeometricSpreadingParameters(Rrefi, [γ1d, γ2d])
        geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0,1]), :Piecewise)

        anef = AnelasticAttenuationParameters(Q0f, ηf)
        aned = AnelasticAttenuationParameters(Q0d, ηd)
        anem = AnelasticAttenuationParameters(Q0f, ηd)

        sat = NearSourceSaturationParameters(:BT15)

        pathf = PathParameters(geof, sat, anef)
        pathd = PathParameters(geod, sat, aned)
        pathm = PathParameters(geom, sat, anem)

        sitef = SiteParameters(κ0f)
        sited = SiteParameters(κ0d)

        fasf = FourierParameters(srcf, pathf, sitef)
        fasd = FourierParameters(srcd, pathd, sited)
        fasm = FourierParameters(srcf, pathm, sited)


        @code_warntype fourier_constant(srcf)
        @time fourier_constant(srcf)
        @time fourier_constant(fasf)

        @code_warntype fourier_constant(srcd)
        @time fourier_constant(srcd)
        @time fourier_constant(fasd)

        f = 1.0
        m = 6.0

        @code_warntype fourier_source_shape(f, m, srcf)
        @code_warntype fourier_source_shape(f, m, srcd)
        @time fourier_source_shape(f, m, srcf)
        @time fourier_source_shape(f, m, fasf)
        @time fourier_source_shape(f, m, srcd)
        @time fourier_source_shape(f, m, fasd)

        fa, fb, ε = corner_frequency(m, srcf)
        @code_warntype fourier_source_shape(f, fa, fb, ε, srcf.model)
        @time fourier_source_shape(f, fa, fb, ε, srcf.model)

        fa, fb, ε = corner_frequency(m, srcd)
        @code_warntype fourier_source_shape(f, fa, fb, ε, srcd.model)
        @time fourier_source_shape(f, fa, fb, ε, srcd.model)

        @code_warntype fourier_source(f, m, srcf)
        @code_warntype fourier_source(f, m, srcd)
        @time fourier_source(f, m, srcf)
        @time fourier_source(f, m, fasf)
        @time fourier_source(f, m, srcd)
        @time fourier_source(f, m, fasd)

        f = 10.0
        r = 100.0

        @code_warntype fourier_path(f, r, geof, anef)
        @code_warntype fourier_path(f, r, geod, aned)
        @code_warntype fourier_path(f, r, geom, anef)
        @code_warntype fourier_path(f, r, pathf)
        @code_warntype fourier_path(f, r, pathd)
        @code_warntype fourier_path(f, r, pathm)
        @code_warntype fourier_path(f, r, fasf)
        @code_warntype fourier_path(f, r, fasd)
        @code_warntype fourier_path(f, r, fasm)

        @time fourier_path(f, r, fasf)
        @time fourier_path(f, r, fasd)
        @time fourier_path(f, r, fasm)


        f = 10.0
        @code_warntype fourier_attenuation(f, r, anef, sitef)
        @code_warntype fourier_attenuation(f, r, aned, sited)
        @code_warntype fourier_attenuation(f, r, anef, sited)
        @code_warntype fourier_attenuation(f, r, pathf, sitef)
        @code_warntype fourier_attenuation(f, r, pathd, sited)
        @code_warntype fourier_attenuation(f, r, pathm, sited)
        @code_warntype fourier_attenuation(f, r, fasf)
        @code_warntype fourier_attenuation(f, r, fasd)
        @code_warntype fourier_attenuation(f, r, fasm)

        @time fourier_attenuation(f, r, fasf)
        @time fourier_attenuation(f, r, fasd)
        @time fourier_attenuation(f, r, fasm)

        @code_warntype fourier_site(f, sitef)
        @code_warntype fourier_site(f, sited)
        @time fourier_site(f, sitef)
        @time fourier_site(f, sited)

        @time fourier_site(f, fasf)
        @time fourier_site(f, fasd)
        @time fourier_site(f, fasm)

        f = 1.0
        m = 6.0
        r = 10.0
        r_psf = equivalent_point_source_distance(r, m, fasf)
        r_psd = equivalent_point_source_distance(r, m, fasd)
        r_psm = equivalent_point_source_distance(r, m, fasm)

        @code_warntype fourier_spectral_ordinate(f, m, r_psf, fasf)
        @code_warntype fourier_spectral_ordinate(f, m, r_psd, fasd)
        @code_warntype fourier_spectral_ordinate(f, m, r_psm, fasm)

        @time fourier_spectral_ordinate(f, m, r_psf, fasf)
        @time fourier_spectral_ordinate(f, m, r_psd, fasd)
        @time fourier_spectral_ordinate(f, m, r_psm, fasm)

        fi = [ 0.01, 0.1, 1.0, 10.0, 100.0 ]
        size(fi)
        length(fi)

        get_parametric_type(geof)
        get_parametric_type(srcf)
        get_parametric_type(anef)
        get_parametric_type(aned)
        get_parametric_type(fasm)

        @code_warntype fourier_spectrum(fi, m, r_psf, fasf)
        @code_warntype fourier_spectrum(fi, m, r_psf, fasd)
        @code_warntype fourier_spectrum(fi, m, r_psd, fasd)
        @code_warntype fourier_spectrum(fi, m, r_psd, fasm)

        @time Afif = fourier_spectrum(fi, m, r_psf, fasf)
        @time Afid = fourier_spectrum(fi, m, r_psf, fasd)
        @time Afim = fourier_spectrum(fi, m, r_psd, fasm)

        @code_warntype fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        @code_warntype fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        @code_warntype fourier_spectrum!(Afim, fi, m, r_psd, fasm)

        @time fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        @time fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        @time fourier_spectrum!(Afim, fi, m, r_psd, fasm)

        @code_warntype squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        @code_warntype squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        @code_warntype squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)

        @time squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        @time squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        @time squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)

        @code_warntype combined_kappa_frequency(r_psf, fasf)
        @code_warntype combined_kappa_frequency(r_psd, fasd)

    end

    @testset "RVT" begin

        @testset "Integration" begin

            n = 200
            @time xi, wi = gausslegendre(n)
            @time xi, wi = gausslaguerre(n)
            @time xi, wi = gausslobatto(n)


            Δσf = 100.0
            γ1f = 1.158
            γ2f = 0.5
            Q0f = 212.5
            ηf = 0.65
            κ0f = 0.038

            Δσd = Dual{Float64}(Δσf)
            γ1d = Dual{Float64}(γ1f)
            γ2d = Dual{Float64}(γ2f)
            Q0d = Dual{Float64}(Q0f)
            ηd = Dual{Float64}(ηf)
            κ0d = Dual{Float64}(κ0f)

            srcf = SourceParameters(Δσf)
            srcd = SourceParameters(Δσd)

            Rrefi = [1.0, 50.0, Inf]
            geof = GeometricSpreadingParameters(Rrefi, [γ1f, γ2f])
            geod = GeometricSpreadingParameters(Rrefi, [γ1d, γ2d])
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0,1]), :Piecewise)

            anef = AnelasticAttenuationParameters(Q0f, ηf)
            aned = AnelasticAttenuationParameters(Q0d, ηd)
            anem = AnelasticAttenuationParameters(Q0f, ηd)

            sat = NearSourceSaturationParameters(:BT15)

            pathf = PathParameters(geof, sat, anef)
            pathd = PathParameters(geod, sat, aned)
            pathm = PathParameters(geom, sat, anem)

            sitef = SiteParameters(κ0f, :Boore2016)
            sited = SiteParameters(κ0d, :Boore2016)

            fasf = FourierParameters(srcf, pathf, sitef)
            fasd = FourierParameters(srcd, pathd, sited)
            fasm = FourierParameters(srcf, pathm, sited)

            m = 5.5
            r = 10.88
            r_psf = equivalent_point_source_distance(r, m, fasf)
            r_psd = equivalent_point_source_distance(r, m, fasd)
            r_psm = equivalent_point_source_distance(r, m, fasm)

            sdof = Oscillator(100.0)


            function gauss_interval(n, fmin, fmax, integrand::Function)
                xi, wi = gausslegendre(n)
                ifi = @. integrand( (fmax-fmin)/2 * xi + (fmin+fmax)/2 )
                return (fmax-fmin)/2 * dot( wi, ifi )
            end

            function gauss_intervals(fun::Function, n, flim...)
                xi, wi = gausslegendre(n)
                ii = 0.0
                for i in 2:length(flim)
                    ii += (flim[i]-flim[i-1])/2 * dot( wi, fun.( (flim[i]-flim[i-1])/2 * xi .+ (flim[i]+flim[i-1])/2 ) )
                end
                return ii
            end

            function trapezoidal_rule(x::Vector, y::Vector)
                return (x[2] - x[1]) * ( sum(y) - (y[1] + y[end])/2 )
            end

            function trapezoidal(fun::Function, n, flim...)
                ii = 0.0
                for i in 2:length(flim)
                    xi = collect(range(flim[i-1], stop=flim[i], length=n))
                    yi = fun.(xi)
                    ii += trapezoidal_rule(xi, yi)
                end
                return ii
            end


            # Boore comparison (assume his are cgs units)
            # ps2db(f) = ( (2π * sdof.f_n) / ( (2π * f)^2 ) )^2 * 1e-4
            # ps2db(f) = ( 1.0 / ( 2π * f^2 * sdof.f_n ) )^2
            ps2db(f) = ( 1.0 / ( 2π*sdof.f_n ) )^2 * 1e4

            dbm0_integrand(f) = squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof) * ps2db(f)
            dbm0ln_integrand(lnf) = exp(lnf) * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm1_integrand(f) = (2π*f) * squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof) * ps2db(f)
            dbm1ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf)) * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm2_integrand(f) = (2π*f)^2 * squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof) * ps2db(f)
            dbm2ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^2  * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm4_integrand(f) = (2π*f)^4 * squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof) * ps2db(f)
            dbm4ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^4 * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            @time igk = 2*quadgk(dbm0_integrand, 0.0, Inf)[1]
            @time igk = 2*quadgk(dbm0_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*gauss_intervals(dbm0ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            @time igk = 2*quadgk(dbm1_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*gauss_intervals(dbm1ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            @time igk = 2*quadgk(dbm2_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*gauss_intervals(dbm2ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            @time igk = 2*quadgk(dbm4_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*gauss_intervals(dbm4ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)


            # using DifferentialEquations
            #
            # u0 = 0.0
            # tspan = (0.0, 300.0)
            # f0(u,p,t) = dbm0_integrand(t)
            # prob = ODEProblem(f0, u0, tspan)
            # sol = solve(prob, RK4())
            # 2*sol.u[end]


            m0_integrand(f) = squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof)
            m0ln_integrand(lnf) = exp(lnf) * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof)

            m1_integrand(f) = (2π*f) * squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof)
            m1ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf)) * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof)

            m2_integrand(f) = (2π*f)^2 * squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof)
            m2ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^2  * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof)

            m4_integrand(f) = (2π*f)^4 * squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * squared_transfer(f, sdof)
            m4ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^4 * squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * squared_transfer(exp(lnf), sdof)


            @time igk = quadgk(m0_integrand, exp(-7.0), exp(7.0))[1]
            @time igle = gauss_interval(2000, 0.0, 300.0, m0_integrand)
            @time igleln = gauss_interval(750, -7.0, 7.0, m0ln_integrand)

            @time iglelnm = gauss_intervals(m0ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            # @test igk ≈ igle rtol=1e-2
            @test igk ≈ igleln rtol=1e-4
            @test igk ≈ iglelnm rtol=1e-4

            lnfi = log.([ 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, sdof.f_n ])
            sort!(lnfi)

            @time igk = quadgk(m1_integrand, exp(-7.0), exp(7.0))[1]
            @time igle = gauss_interval(1500, 0.0, 300.0, m1_integrand)
            @time igleln = gauss_interval(750, -7.0, 7.0, m1ln_integrand)
            @time iglelnm = gauss_intervals(m1ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            @time iglelnm = gauss_intervals(m1ln_integrand, 30, lnfi...)

            @time itr = trapezoidal(m1ln_integrand, 60, lnfi...)

            # @test igk ≈ igle rtol=1e-2
            @test igk ≈ iglelnm rtol=1e-5
            @test igk ≈ itr rtol=1e-3


            @time igk = quadgk(m2_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = gauss_intervals(m2ln_integrand, 30, lnfi...)
            @time itr = trapezoidal(m2ln_integrand, 60, lnfi...)

            @test igk ≈ iglelnm rtol=1e-5
            @test igk ≈ itr rtol=1e-3

            @time igk = quadgk(m4_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = gauss_intervals(m4ln_integrand, 30, lnfi...)
            @time itr = trapezoidal(m4ln_integrand, 60, lnfi...)

            @test igk ≈ iglelnm rtol=1e-4
            @test igk ≈ itr rtol=1e-3




            integrand(x) = sin(x)

            intervals = 101
            x_min = 0.0
            x_max = 2π
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = simpsons_rule(xx, yy)
            ist = trapezoidal_rule(xx, yy)

            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk atol=10*eps()
            @test ist ≈ igk atol=10*eps()


            m = 6.0
            r = 100.0
            src = SourceParameters(100.0)
            geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
            sat = NearSourceSaturationParameters(:BT15)
            ane = AnelasticAttenuationParameters(200.0, 0.4, :Rrup)
            path = PathParameters(geo, sat, ane)
            site = SiteParameters(0.039)
            fas = FourierParameters(src, path, site)
            sdof = Oscillator(1.0)

            integrand(f) = squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = sdof.f_n/1.1
            x_max = sdof.f_n*1.1
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = simpsons_rule(xx, yy)
            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk rtol=1e-6
            @test isr ≈ igk atol=1e-6


            integrand(f) = squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = 100.0
            x_max = 200.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = simpsons_rule(xx, yy)
            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk rtol=1e-3
            @test isr ≈ igk atol=1e-6


            integrand(f) = (2π*f)^4 * squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = 100.0
            x_max = 200.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = simpsons_rule(xx, yy)
            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk rtol=1e-3
            @test isr ≈ igk atol=1e-6

            x_min = 300.0
            x_max = 500.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = simpsons_rule(xx, yy)
            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk rtol=1e-3
            @test isr ≈ igk atol=1e-6

        end

        @testset "Spectral Moments" begin

            Δσf = 100.0
            γ1f = 1.0
            γ2f = 0.5
            Q0f = 200.0
            ηf = 0.4
            κ0f = 0.039

            Δσd = Dual{Float64}(Δσf)
            γ1d = Dual{Float64}(γ1f)
            γ2d = Dual{Float64}(γ2f)
            Q0d = Dual{Float64}(Q0f)
            ηd = Dual{Float64}(ηf)
            κ0d = Dual{Float64}(κ0f)

            srcf = SourceParameters(Δσf)
            srcd = SourceParameters(Δσd)

            Rrefi = [1.0, 50.0, Inf]
            geof = GeometricSpreadingParameters(Rrefi, [γ1f, γ2f])
            geod = GeometricSpreadingParameters(Rrefi, [γ1d, γ2d])
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0,1]), :Piecewise)

            anef = AnelasticAttenuationParameters(Q0f, ηf)
            aned = AnelasticAttenuationParameters(Q0d, ηd)
            anem = AnelasticAttenuationParameters(Q0f, ηd)

            sat = NearSourceSaturationParameters(:BT15)

            pathf = PathParameters(geof, sat, anef)
            pathd = PathParameters(geod, sat, aned)
            pathm = PathParameters(geom, sat, anem)

            sitef = SiteParameters(κ0f)
            sited = SiteParameters(κ0d)

            fasf = FourierParameters(srcf, pathf, sitef)
            fasd = FourierParameters(srcd, pathd, sited)
            fasm = FourierParameters(srcf, pathm, sited)

            m = 6.0
            r = 10.0
            r_psf = equivalent_point_source_distance(r, m, fasf)
            r_psd = equivalent_point_source_distance(r, m, fasd)
            r_psm = equivalent_point_source_distance(r, m, fasm)

            sdof = Oscillator(0.1)

            order = 0
            @time spectral_moment(order, m, r_psf, fasf, sdof)
            @time spectral_moment(order, m, r_psd, fasd, sdof)
            @time spectral_moment(order, m, r_psm, fasm, sdof)

            @code_warntype spectral_moment(order, m, r_psf, fasf, sdof)
            @code_warntype spectral_moment(order, m, r_psd, fasd, sdof)
            @code_warntype spectral_moment(order, m, r_psm, fasm, sdof)

            @time spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            @time spectral_moments([0, 1, 2, 4], m, r_psd, fasd, sdof)
            @time spectral_moments([0, 1, 2, 4], m, r_psm, fasm, sdof)

            @code_warntype spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            @code_warntype spectral_moments([0, 1, 2, 4], m, r_psd, fasd, sdof)
            @code_warntype spectral_moments([0, 1, 2, 4], m, r_psm, fasm, sdof)

            @time smi = spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            @time smigk = spectral_moments_gk([0, 1, 2, 4], m, r_psf, fasf, sdof)

            @test all(isapprox.(smi, smigk, rtol=1e-3))
            # [ smi smigk ]

            sdof = Oscillator(1/3)
            @time smi = spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            @time smigk = spectral_moments_gk([0, 1, 2, 4], m, r_psf, fasf, sdof)

            @test all(isapprox.(smi, smigk, rtol=1e-3))

        end

        @testset "Peak Factor" begin

            Δσf = 100.0
            γ1f = 1.0
            γ2f = 0.5
            Q0f = 200.0
            ηf = 0.4
            κ0f = 0.039

            Δσd = Dual{Float64}(Δσf)
            γ1d = Dual{Float64}(γ1f)
            γ2d = Dual{Float64}(γ2f)
            Q0d = Dual{Float64}(Q0f)
            ηd = Dual{Float64}(ηf)
            κ0d = Dual{Float64}(κ0f)

            srcf = SourceParameters(Δσf)
            srcd = SourceParameters(Δσd)

            Rrefi = [1.0, 50.0, Inf]
            geof = GeometricSpreadingParameters(Rrefi, [γ1f, γ2f])
            geod = GeometricSpreadingParameters(Rrefi, [γ1d, γ2d])
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0,1]), :Piecewise)

            anef = AnelasticAttenuationParameters(Q0f, ηf)
            aned = AnelasticAttenuationParameters(Q0d, ηd)
            anem = AnelasticAttenuationParameters(Q0f, ηd)

            sat = NearSourceSaturationParameters(:BT15)

            pathf = PathParameters(geof, sat, anef)
            pathd = PathParameters(geod, sat, aned)
            pathm = PathParameters(geom, sat, anem)

            sitef = SiteParameters(κ0f)
            sited = SiteParameters(κ0d)

            fasf = FourierParameters(srcf, pathf, sitef)
            fasd = FourierParameters(srcd, pathd, sited)
            fasm = FourierParameters(srcf, pathm, sited)


            m = 7.0
            r = 1.0

            r_psf = equivalent_point_source_distance(r, m, fasf)
            r_psd = equivalent_point_source_distance(r, m, fasd)
            r_psm = equivalent_point_source_distance(r, m, fasm)

            sdof = Oscillator(1.0)

            @time pfps = peak_factor_dk80(m, r_psf, fasf, sdof)
            @time pfpsn = peak_factor_dk80(m, r_psf, fasf, sdof, nodes=30)
            @time pfgk = peak_factor_dk80_gk(m, r_psf, fasf, sdof)

            @test pfps ≈ pfgk rtol=1e-6
            @test pfpsn ≈ pfgk rtol=1e-6

            @time pfps = peak_factor_cl56(m, r_psf, fasf, sdof)
            @time pfpsn = peak_factor_cl56(m, r_psf, fasf, sdof, nodes=40)
            @time pfgk = peak_factor_cl56_gk(m, r_psf, fasf, sdof)

            @test pfps ≈ pfgk rtol=1e-5
            @test pfpsn ≈ pfgk rtol=1e-5


        end

        @testset "Response Spectra" begin

            Δσf = 100.0
            γ1f = 1.0
            γ2f = 0.5
            Q0f = 200.0
            ηf = 0.4
            κ0f = 0.039

            Δσd = Dual{Float64}(Δσf)
            γ1d = Dual{Float64}(γ1f)
            γ2d = Dual{Float64}(γ2f)
            Q0d = Dual{Float64}(Q0f)
            ηd = Dual{Float64}(ηf)
            κ0d = Dual{Float64}(κ0f)

            srcf = SourceParameters(Δσf)
            srcd = SourceParameters(Δσd)

            Rrefi = [1.0, 50.0, Inf]
            geof = GeometricSpreadingParameters(Rrefi, [γ1f, γ2f])
            geod = GeometricSpreadingParameters(Rrefi, [γ1d, γ2d])
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0,1]), :Piecewise)

            anef = AnelasticAttenuationParameters(Q0f, ηf)
            aned = AnelasticAttenuationParameters(Q0d, ηd)
            anem = AnelasticAttenuationParameters(Q0f, ηd)

            sat = NearSourceSaturationParameters(:BT15)

            pathf = PathParameters(geof, sat, anef)
            pathd = PathParameters(geod, sat, aned)
            pathm = PathParameters(geom, sat, anem)

            sitef = SiteParameters(κ0f)
            sited = SiteParameters(κ0d)

            fasf = FourierParameters(srcf, pathf, sitef)
            fasd = FourierParameters(srcd, pathd, sited)
            fasm = FourierParameters(srcf, pathm, sited)


            m = 7.0
            r = 1.0

            r_psf = equivalent_point_source_distance(r, m, fasf)
            r_psd = equivalent_point_source_distance(r, m, fasd)
            r_psm = equivalent_point_source_distance(r, m, fasm)

            # Ti = exp10.(range(-2.0, stop=1.0, length=31))
            Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]

            rvt = RandomVibrationParameters()

            @time Saif = rvt_response_spectrum(Ti, m, r_psf, fasf, rvt)
            # @btime Saif = rvt_response_spectrum(Ti, m, r_psf, fasf, rvt)
            @time Said = rvt_response_spectrum(Ti, m, r_psd, fasd, rvt)
            # @btime Said = rvt_response_spectrum(Ti, m, r_psd, fasd, rvt)
            @time Saim = rvt_response_spectrum(Ti, m, r_psm, fasm, rvt)

            @code_warntype rvt_response_spectrum(Ti, m, r_psf, fasf, rvt)
            @code_warntype rvt_response_spectrum(Ti, m, r_psd, fasd, rvt)
            @code_warntype rvt_response_spectrum(Ti, m, r_psm, fasm, rvt)


            function spectral_slope_rtp(x::Vector, r_rup, T, par::Vector)
                src = SourceParameters(exp(par[1]))
                geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [ par[2], 0.5 ], :CY14mod)
                heff = par[3] * exp( par[4]*m )
                sat = NearSourceSaturationParameters(heff, 1)
                ane = AnelasticAttenuationParameters(par[5], par[6], :Rrup)
                path = PathParameters(geo, sat, ane)
                site = SiteParameters(0.039)
                fas = FourierParameters(src, path, site)
                rvt = RandomVibrationParameters()

                r_ps = equivalent_point_source_distance(r_rup, x[1], fas)

                return log(rvt_response_spectral_ordinate(T, x[1], r_ps, fas, rvt))
            end

            r_rup = 1.0
            T = 0.01
            hα = 0.023
            hβ = 0.4
            γ = 1.2
            γ*hβ
            par = [ 5.25, γ, hα, hβ, 185.0, 0.66 ]


            x = [ 7.8 ]
            h = hα * exp( hβ * x[1] )
            r_ps = r_rup + h

            spectral_slope(x) = spectral_slope_rtp(x, r_rup, T, par)
            exp(spectral_slope(x))
            gspectral_slope(x) = ForwardDiff.gradient(spectral_slope, x)
            gspectral_slope(x)[1]

            lnSa1 = spectral_slope([7.0])
            lnSa2 = spectral_slope([8.40])
            (lnSa2 - lnSa1)/1.4



            r_rup = 1.0
            T = 0.01
            hα = 0.023
            hβ = 0.88
            γ = 1.2
            par = [ 5.25, γ, hα, hβ, 185.0, 0.66 ]

            m = 8.4

            src = SourceParameters(exp(par[1]))
            geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [ par[2], 0.5 ], :CY14mod)
            heff = hα * exp( hβ * m )
            sat = NearSourceSaturationParameters(heff, 1)
            ane = AnelasticAttenuationParameters(par[5], par[6], :Rrup)
            path = PathParameters(geo, sat, ane)
            fas = FourierParameters(src, path, site)
            # equivalent_point_source_distance(r_rup, m, fas)
            r_ps = r_rup + heff
            site = SiteParameters(0.039)
            rvt = RandomVibrationParameters()
            sdof = Oscillator(1.0/T)

            rvt_response_spectral_ordinate(T, m, r_ps, fas, rvt)

            fc, fb, fe = corner_frequency(m, fas)
            fk = combined_kappa_frequency(r_rup, ane, site)

            m0 = spectral_moment(0, m, r_ps, fas, sdof)

            Afc = fourier_spectral_ordinate(fc*1.5, m, r_ps, fas)
            Afk = fourier_spectral_ordinate(fk, m, r_ps, fas)

            fn = 100.0
            transfer(fn/3, Oscillator(fn))


            mi = collect(range(7.0, stop=8.5, step=0.1))
            pfi = zeros(length(mi))
            for i in 1:length(mi)
                heff = hα * exp( hβ * mi[i] )
                sat = NearSourceSaturationParameters(heff, 1)
                ane = AnelasticAttenuationParameters(par[5], par[6], :Rrup)
                path = PathParameters(geo, sat, ane)
                fas = FourierParameters(src, path, site)
                r_ps = r_rup + heff
                pfi[i] = peak_factor(mi[i], r_ps, fas, sdof, rvt)
            end

            [ mi pfi ]

        end

    end

end
