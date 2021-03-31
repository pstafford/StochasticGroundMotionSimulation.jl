using StochasticGroundMotionSimulation
using Test
using ForwardDiff
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
        # @code_warntype Oscillator(1.0)
        # @code_warntype period(sdof)
        # @code_warntype transfer(1.0, sdof)

        # @code_warntype rvt_response_spectral_ordinate(Ti[1], m, r, fas, rvt)
        # @code_warntype rvt_response_spectrum(Ti, m, r, fas, rvt)
        # @time Sai = rvt_response_spectrum(Ti, m, r, fas, rvt)

    end

    @testset "Source" begin

        m = 6.0
        Δσ = 100.0
        β = 3.5

        # @code_warntype magnitude_to_moment(m)
        #
        # @code_warntype corner_frequency_brune(m, Δσ)
        # @code_warntype corner_frequency_brune(m, Δσ, β)
        # @code_warntype corner_frequency_atkinson_silva_2000(m)

        srcf = SourceParameters(Δσ)
        srcd = SourceParameters(Dual{Float64}(Δσ))

        # @code_warntype corner_frequency(m, srcf)
        # @code_warntype corner_frequency(m, srcd)
        # @code_warntype corner_frequency(Dual(m), srcf)
        # @code_warntype corner_frequency(Dual(m), srcd)

        T = get_parametric_type(srcf)
        @test T == Float64
        T = get_parametric_type(srcd)
        @test T <: Dual


        faf, fbf, fεf = corner_frequency(m, srcf)
        fad, fbd, fεd = corner_frequency(m, srcd)
        # @time faf, fbf, fεf = corner_frequency(m, srcf)
        # @time fad, fbd, fεd = corner_frequency(m, srcd)

        @test faf == fad.value

        srcf = SourceParameters(Δσ, :Atkinson_Silva_2000)
        srcd = SourceParameters(Dual{Float64}(Δσ), :Atkinson_Silva_2000)

        faf, fbf, fεf = corner_frequency(m, srcf)
        fad, fbd, fεd = corner_frequency(m, srcd)
        # @time faf, fbf, fεf = corner_frequency(m, srcf)
        # @time fad, fbd, fεd = corner_frequency(m, srcd)

        @test faf == fad.value
        @test fbf == fbd.value
        @test fεf == fεd.value

    end

    @testset "Path" begin

        Rrefi = [ 1.0, 50.0, Inf ]
        γi = [ 1.0, 0.5 ]
        geof = GeometricSpreadingParameters(Rrefi, γi)
        geod = GeometricSpreadingParameters(Rrefi, [ 0.5 ], [ Dual{Float64}(1.0) ], BitVector([1,0]), :Piecewise)

        T = get_parametric_type(geof)
        @test T == Float64
        T = get_parametric_type(geod)
        @test T <: Dual


        sat = NearSourceSaturationParameters(:BT15)

        T = get_parametric_type(sat)
        @test T == Float64

        Q0 = 200.0
        anef = AnelasticAttenuationParameters(Q0)
        aned = AnelasticAttenuationParameters(Dual{Float64}(Q0))

        T = get_parametric_type(anef)
        @test T == Float64
        T = get_parametric_type(aned)
        @test T <: Dual

        pathf = PathParameters(geof, sat, anef)
        pathd = PathParameters(geod, sat, aned)

        T = get_parametric_type(pathf)
        @test T == Float64
        T = get_parametric_type(pathd)
        @test T <: Dual

        r = 10.0
        m = 6.0

        # @code_warntype near_source_saturation(m, pathf.saturation)
        # @code_warntype near_source_saturation(m, pathd.saturation)
        hf = near_source_saturation(m, pathf.saturation)
        hd = near_source_saturation(m, pathd.saturation)
        if get_parametric_type(pathd.saturation) <: Dual
            @test hf == hd.value
        else
            @test hf == hd
        end

        # @code_warntype near_source_saturation(m, pathf)
        # @code_warntype near_source_saturation(m, pathd)
        hf = near_source_saturation(m, pathf)
        hd = near_source_saturation(m, pathd)
        if get_parametric_type(pathd.saturation) <: Dual
            @test hf == hd.value
        else
            @test hf == hd
        end

        # @code_warntype equivalent_point_source_distance(r, m, pathf)
        # @code_warntype equivalent_point_source_distance(r, m, pathd)
        r_psf = equivalent_point_source_distance(r, m, pathf)
        r_psd = equivalent_point_source_distance(r, m, pathd)

        # @code_warntype geometric_spreading_piecewise(r, geof)
        # @code_warntype geometric_spreading_piecewise(r, geod)
        grf = geometric_spreading_piecewise(r, geof)
        grd = geometric_spreading_piecewise(r, geod)
        @test grf == grd.value

        # @code_warntype geometric_spreading_cy14(r, geof)
        # @code_warntype geometric_spreading_cy14(r, geod)
        grf = geometric_spreading_cy14(r, geof)
        grd = geometric_spreading_cy14(r, geod)
        @test grf == grd.value

        # @code_warntype geometric_spreading(r, geof)
        # @code_warntype geometric_spreading(r, geod)

        grf = geometric_spreading(r, m, geof, sat)
        grd = geometric_spreading(r, m, geod, sat)
        @test grf == grd.value

        f = 1.0
        r = 100.0
        # @code_warntype anelastic_attenuation(f, r, anef)
        # @code_warntype anelastic_attenuation(f, r, aned)
        qrf = anelastic_attenuation(f, r, anef)
        qrd = anelastic_attenuation(f, r, aned)
        @test qrf == qrd.value

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
        # @code_warntype site_amplification(f, site0f)
        # @code_warntype site_amplification(f, site0d)
        Sff = site_amplification(f, site0f)
        Sfd = site_amplification(f, site0d)
        @test Sff == Sfd

        Af0f = site_amplification(f, site0f)
        Af1f = site_amplification(f, siteAf)
        Af2f = site_amplification(f, siteBf)
        Af3f = site_amplification(f, siteUf)
        Af4f = site_amplification(f, siteNf)

        @test Af0f == Af1f
        @test Af3f == 1.0
        @test Af2f < Af1f
        @test isnan(Af4f)

        Af0d = site_amplification(f, site0d)
        Af1d = site_amplification(f, siteAd)
        Af2d = site_amplification(f, siteBd)
        Af3d = site_amplification(f, siteUd)
        Af4d = site_amplification(f, siteNd)

        @test Af0d == Af1d
        @test Af3d == 1.0
        @test Af2d < Af1d
        @test isnan(Af4d)

        @test Af0f == Af0d

        f0 = 80.0
        f1 = 100.0
        Af0 = site_amplification(f0, siteBf)
        Af1 = site_amplification(f1, siteBf)

        @test Af0 == Af1

        f = 10.0
        # @code_warntype kappa_filter(f, siteAf)
        # @code_warntype kappa_filter(f, siteAd)
        Kff = kappa_filter(f, siteAf)
        Kfd = kappa_filter(f, siteAd)
        @test Kff == Kfd.value

    end

    @testset "Oscillator" begin
        ζ = 0.05
        f_n = 1.0
        sdof = Oscillator(f_n, ζ)

        @test f_n ≈ 1.0/period(sdof)

        @test transfer(0.5, sdof)^2 ≈ squared_transfer(0.5, sdof)

        fi = [ 0.5, 1.0, 2.0 ]

        # @code_warntype transfer(fi, sdof)

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

        faf, fbf, εf = corner_frequency(m, srcf)
        # @code_warntype boore_thompson_2014(m, 0.0, srcf)
        Dsf = boore_thompson_2014(m, 0.0, srcf)
        @test Dsf ≈ 1.0/faf
        @test isnan(fbf)
        @test isnan(εf)
        fad, fbd, εd = corner_frequency(m, srcd)
        # @code_warntype boore_thompson_2014(m, 0.0, srcd)
        Dsd = boore_thompson_2014(m, 0.0, srcd)
        @test Dsd ≈ 1.0/fad
        @test isnan(fbd)
        @test isnan(εd)

        # @code_warntype boore_thompson_2014(m, 0.0, fasf)
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
        @test isnan(fb) == false
        @test isnan(ε) == false
        @test fa < fb

        # test gradient of BT14 duration model w.r.t. magnitude
        h = 0.05
        m1 = 8.0
        m2 = m1 + h
        r_ps1 = 1.0 + near_source_saturation(m1, fasf)
        r_ps2 = 1.0 + near_source_saturation(m2, fasf)
        Dex1 = boore_thompson_2014(m1, r_ps1, fasf)
        Dex2 = boore_thompson_2014(m2, r_ps2, fasf)
        fdg = log(Dex2/Dex1)/h

        d(x) = log(boore_thompson_2014(x[1], 1.0 + near_source_saturation(x[1], fasf), fasf))
        gd(x) = ForwardDiff.gradient(d, x)
        adg = gd([8.0])[1]

        @test fdg ≈ adg atol=1e-2


        rvt = RandomVibrationParameters(:BT14)
        # @code_warntype excitation_duration(m, r, fasf, rvt)
        # @code_warntype excitation_duration(m, r, fasd, rvt)
        Dexf = excitation_duration(m, r, fasf, rvt)
        Dexd = excitation_duration(m, r, fasd, rvt)
        @test Dexf == Dexd.value


        c11 = [ 8.4312e-01, -2.8671e-02, 2.0,  1.7316e+00,  1.1695e+00,  2.1671e+00,  9.6224e-01 ]
        c11f = boore_thompson_2012_coefs(1, 1)
        @test c11f[1] == c11[1]
        @test all(isapprox.(c11, c11f))

        # @code_warntype boore_thompson_2012_coefs(1, 1)

        m = 8.0
        r = 1.0
        Dex = boore_thompson_2014(m, r, srcf)
        # get the oscillator period
        sdof = Oscillator(100.0)
        T_n = period(sdof)
        ζ = sdof.ζ_n
        # define the η parameter as T_n/Dex
        η = T_n / Dex

        # @time c = boore_thompson_2012_coefs(1, 1)
        # @time boore_thompson_2012_base(η, c, ζ)
        #
        # @time boore_thompson_2012(m, r, srcf, sdof, rvt)
        # @code_warntype boore_thompson_2012(m, r, srcf, sdof, rvt)
        #
        # @time boore_thompson_2012(m, r, srcd, sdof, rvt)
        # @code_warntype boore_thompson_2012(m, r, srcd, sdof, rvt)

        sdof = Oscillator(1.0)

        Drms, Dex, Dratio = boore_thompson_2012(m, r, fas, sdof, rvt)
        Dex0 = boore_thompson_2014(m, r, fas)

        @test Dex == Dex0

        # @code_warntype rms_duration(m, r, srcf, path, sdof, rvt)
        # @code_warntype rms_duration(m, r, srcd, path, sdof, rvt)
        # @code_warntype rms_duration(m, r, fasf, sdof, rvt)
        # @code_warntype rms_duration(m, r, fasd, sdof, rvt)

        # @time Drmsf, Dexf, Dratiof = rms_duration(m, r, fasf, sdof, rvt)
        # @time Drmsd, Dexd, Dratiod = rms_duration(m, r, fasd, sdof, rvt)
        Drmsf, Dexf, Dratiof = rms_duration(m, r, fasf, sdof, rvt)
        Drmsd, Dexd, Dratiod = rms_duration(m, r, fasd, sdof, rvt)

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

        Tf = get_parametric_type(fasf)
        Td = get_parametric_type(fasd)
        Tm = get_parametric_type(fasm)
        @test Tf == Float64
        @test Td <: Dual
        @test Tm <: Dual
        @test Td == Tm

        # @code_warntype fourier_constant(srcf)
        Cfs = fourier_constant(srcf)
        Cff = fourier_constant(fasf)
        @test Cfs == Cff

        # @code_warntype fourier_constant(srcd)
        Cfsd = fourier_constant(srcd)
        Cffd = fourier_constant(fasd)
        @test Cfsd == Cffd

        f = 0.001
        m = 6.0

        # @code_warntype fourier_source_shape(f, m, srcf)
        # @code_warntype fourier_source_shape(f, m, srcd)
        Affs = fourier_source_shape(f, m, srcf)
        Afff = fourier_source_shape(f, m, fasf)
        Afds = fourier_source_shape(f, m, srcd)
        Afdf = fourier_source_shape(f, m, fasd)
        @test Affs == Afds.value
        @test Afff == Afdf.value
        @test Affs == Afff
        @test Afds == Afdf

        @test Afff ≈ 1.0 atol=1e-3

        fa, fb, ε = corner_frequency(m, srcf)
        # @code_warntype fourier_source_shape(f, fa, fb, ε, srcf.model)
        Afc = fourier_source_shape(f, fa, fb, ε, srcf.model)
        @test Afc ≈ 1.0 atol=1e-3

        fa, fb, ε = corner_frequency(m, srcd)
        # @code_warntype fourier_source_shape(f, fa, fb, ε, srcd.model)
        Afcd = fourier_source_shape(f, fa, fb, ε, srcd.model)
        @test Afcd ≈ 1.0 atol=1e-3

        # @code_warntype fourier_source(f, m, srcf)
        # @code_warntype fourier_source(f, m, srcd)
        @time fourier_source(f, m, srcf)
        @time fourier_source(f, m, fasf)
        @time fourier_source(f, m, srcd)
        @time fourier_source(f, m, fasd)

        f = 10.0
        r = 100.0

        # @code_warntype fourier_path(f, r, geof, anef)
        # @code_warntype fourier_path(f, r, geod, aned)
        # @code_warntype fourier_path(f, r, geom, anef)
        # @code_warntype fourier_path(f, r, pathf)
        # @code_warntype fourier_path(f, r, pathd)
        # @code_warntype fourier_path(f, r, pathm)
        # @code_warntype fourier_path(f, r, fasf)
        # @code_warntype fourier_path(f, r, fasd)
        # @code_warntype fourier_path(f, r, fasm)

        Pf = fourier_path(f, r, fasf)
        Pd = fourier_path(f, r, fasd)
        Pm = fourier_path(f, r, fasm)
        @test Pf == Pd.value
        @test Pd == Pm


        f = 10.0
        # @code_warntype fourier_attenuation(f, r, anef, sitef)
        # @code_warntype fourier_attenuation(f, r, aned, sited)
        # @code_warntype fourier_attenuation(f, r, anef, sited)
        # @code_warntype fourier_attenuation(f, r, pathf, sitef)
        # @code_warntype fourier_attenuation(f, r, pathd, sited)
        # @code_warntype fourier_attenuation(f, r, pathm, sited)
        # @code_warntype fourier_attenuation(f, r, fasf)
        # @code_warntype fourier_attenuation(f, r, fasd)
        # @code_warntype fourier_attenuation(f, r, fasm)

        Qf = fourier_attenuation(f, r, fasf)
        Qd = fourier_attenuation(f, r, fasd)
        Qm = fourier_attenuation(f, r, fasm)
        @test Qf == Qd.value
        @test Qd == Qm

        # @code_warntype fourier_site(f, sitef)
        # @code_warntype fourier_site(f, sited)
        # @time fourier_site(f, sitef)
        # @time fourier_site(f, sited)

        Sf = fourier_site(f, fasf)
        Sd = fourier_site(f, fasd)
        Sm = fourier_site(f, fasm)
        @test Sf == Sd.value
        @test Sd == Sm

        f = 1.0
        m = 6.0
        r = 10.0
        r_psf = equivalent_point_source_distance(r, m, fasf)
        r_psd = equivalent_point_source_distance(r, m, fasd)
        r_psm = equivalent_point_source_distance(r, m, fasm)

        # @code_warntype fourier_spectral_ordinate(f, m, r_psf, fasf)
        # @code_warntype fourier_spectral_ordinate(f, m, r_psd, fasd)
        # @code_warntype fourier_spectral_ordinate(f, m, r_psm, fasm)

        Af = fourier_spectral_ordinate(f, m, r_psf, fasf)
        Ad = fourier_spectral_ordinate(f, m, r_psd, fasd)
        Am = fourier_spectral_ordinate(f, m, r_psm, fasm)
        @test Af == Ad.value
        @test Ad == Am

        fi = [ 0.01, 0.1, 1.0, 10.0, 100.0 ]

        # @code_warntype fourier_spectrum(fi, m, r_psf, fasf)
        # @code_warntype fourier_spectrum(fi, m, r_psf, fasd)
        # @code_warntype fourier_spectrum(fi, m, r_psd, fasd)
        # @code_warntype fourier_spectrum(fi, m, r_psd, fasm)

        Afif = fourier_spectrum(fi, m, r_psf, fasf)
        Afid = fourier_spectrum(fi, m, r_psf, fasd)
        Afim = fourier_spectrum(fi, m, r_psd, fasm)
        for i in 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))

        # @code_warntype fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        # @code_warntype fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        # @code_warntype fourier_spectrum!(Afim, fi, m, r_psd, fasm)

        fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        fourier_spectrum!(Afim, fi, m, r_psd, fasm)
        for i in 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))

        # @code_warntype squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        # @code_warntype squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        # @code_warntype squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)

        squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)
        for i in 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))

        # @code_warntype combined_kappa_frequency(r_psf, fasf)
        # @code_warntype combined_kappa_frequency(r_psd, fasd)
        fkf = combined_kappa_frequency(r_psf, fasf)
        fkd = combined_kappa_frequency(r_psd, fasd)
        @test fkf == fkd.value

    end

    @testset "RVT" begin

        @testset "Integration" begin

            # n = 200
            # @time xi, wi = gausslegendre(n)
            # @time xi, wi = gausslaguerre(n)
            # @time xi, wi = gausslobatto(n)

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
            @test igk ≈ iglelnm rtol=1e-4

            @time igk = 2*quadgk(dbm1_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*gauss_intervals(dbm1ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol=1e-4

            @time igk = 2*quadgk(dbm2_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*gauss_intervals(dbm2ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol=1e-3

            @time igk = 2*quadgk(dbm4_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*gauss_intervals(dbm4ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol=1e-3


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
            m0f = spectral_moment(order, m, r_psf, fasf, sdof)
            m0d = spectral_moment(order, m, r_psd, fasd, sdof)
            m0m = spectral_moment(order, m, r_psm, fasm, sdof)
            @test m0f == m0d.value
            @test m0d == m0m

            # @code_warntype spectral_moment(order, m, r_psf, fasf, sdof)
            # @code_warntype spectral_moment(order, m, r_psd, fasd, sdof)
            # @code_warntype spectral_moment(order, m, r_psm, fasm, sdof)

            # @time spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            # @time spectral_moments([0, 1, 2, 4], m, r_psd, fasd, sdof)
            # @time spectral_moments([0, 1, 2, 4], m, r_psm, fasm, sdof)

            # @code_warntype spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            # @code_warntype spectral_moments([0, 1, 2, 4], m, r_psd, fasd, sdof)
            # @code_warntype spectral_moments([0, 1, 2, 4], m, r_psm, fasm, sdof)

            smi = spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            smigk = spectral_moments_gk([0, 1, 2, 4], m, r_psf, fasf, sdof)

            @test all(isapprox.(smi, smigk, rtol=1e-3))


            sdof = Oscillator(1/3)
            smi = spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            smigk = spectral_moments_gk([0, 1, 2, 4], m, r_psf, fasf, sdof)

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

            Saif = rvt_response_spectrum(Ti, m, r_psf, fasf, rvt)
            Said = rvt_response_spectrum(Ti, m, r_psd, fasd, rvt)
            Saim = rvt_response_spectrum(Ti, m, r_psm, fasm, rvt)
            for i in 1:length(Ti)
                @test Saif[i] ≈ Said[i].value
            end
            @test all(isapprox.(Said, Saim))

            # @code_warntype rvt_response_spectrum(Ti, m, r_psf, fasf, rvt)
            # @code_warntype rvt_response_spectrum(Ti, m, r_psd, fasd, rvt)
            # @code_warntype rvt_response_spectrum(Ti, m, r_psm, fasm, rvt)


            # function spectral_slope_rtp(x::Vector, r_rup, T, par::Vector)
            #     src = SourceParameters(exp(par[1]))
            #     geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [ par[2], 0.5 ], :CY14mod)
            #     heff = par[3] * exp( par[4]*m )
            #     sat = NearSourceSaturationParameters(heff, 1)
            #     ane = AnelasticAttenuationParameters(par[5], par[6], :Rrup)
            #     path = PathParameters(geo, sat, ane)
            #     site = SiteParameters(0.039)
            #     fas = FourierParameters(src, path, site)
            #     rvt = RandomVibrationParameters()
            #
            #     r_ps = equivalent_point_source_distance(r_rup, x[1], fas)
            #
            #     return log(rvt_response_spectral_ordinate(T, x[1], r_ps, fas, rvt))
            # end
            #
            # r_rup = 1.0
            # T = 0.01
            # hα = 0.023
            # hβ = 0.4
            # γ = 1.2
            # γ*hβ
            # par = [ 5.25, γ, hα, hβ, 185.0, 0.66 ]
            #
            #
            # x = [ 7.8 ]
            # h = hα * exp( hβ * x[1] )
            # r_ps = r_rup + h
            #
            # spectral_slope(x) = spectral_slope_rtp(x, r_rup, T, par)
            # exp(spectral_slope(x))
            # gspectral_slope(x) = ForwardDiff.gradient(spectral_slope, x)
            # gspectral_slope(x)[1]
            #
            # lnSa1 = spectral_slope([7.0])
            # lnSa2 = spectral_slope([8.40])
            # (lnSa2 - lnSa1)/1.4



            # r_rup = 1.0
            # T = 0.01
            # hα = 0.023
            # hβ = 0.88
            # γ = 1.2
            # par = [ 5.25, γ, hα, hβ, 185.0, 0.66 ]
            #
            # m = 8.4
            #
            # src = SourceParameters(exp(par[1]))
            # geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [ par[2], 0.5 ], :CY14mod)
            # heff = hα * exp( hβ * m )
            # sat = NearSourceSaturationParameters(heff, 1)
            # ane = AnelasticAttenuationParameters(par[5], par[6], :Rrup)
            # path = PathParameters(geo, sat, ane)
            # fas = FourierParameters(src, path, site)
            # # equivalent_point_source_distance(r_rup, m, fas)
            # r_ps = r_rup + heff
            # site = SiteParameters(0.039)
            # rvt = RandomVibrationParameters()
            # sdof = Oscillator(1.0/T)
            #
            # rvt_response_spectral_ordinate(T, m, r_ps, fas, rvt)
            #
            # fc, fb, fe = corner_frequency(m, fas)
            # fk = combined_kappa_frequency(r_rup, ane, site)
            #
            # m0 = spectral_moment(0, m, r_ps, fas, sdof)
            #
            # Afc = fourier_spectral_ordinate(fc*1.5, m, r_ps, fas)
            # Afk = fourier_spectral_ordinate(fk, m, r_ps, fas)
            #
            # fn = 100.0
            # transfer(fn/3, Oscillator(fn))
            #
            #
            # mi = collect(range(7.0, stop=8.5, step=0.1))
            # pfi = zeros(length(mi))
            # for i in 1:length(mi)
            #     heff = hα * exp( hβ * mi[i] )
            #     sat = NearSourceSaturationParameters(heff, 1)
            #     ane = AnelasticAttenuationParameters(par[5], par[6], :Rrup)
            #     path = PathParameters(geo, sat, ane)
            #     fas = FourierParameters(src, path, site)
            #     r_ps = r_rup + heff
            #     pfi[i] = peak_factor(mi[i], r_ps, fas, sdof, rvt)
            # end
            #
            # [ mi pfi ]

        end

    end

end
