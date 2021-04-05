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

        T = StochasticGroundMotionSimulation.get_parametric_type(srcf)
        @test T == Float64
        T = StochasticGroundMotionSimulation.get_parametric_type(srcd)
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

        srcd = SourceParameters(Dual(100.0), 3.5, 2.75)
        srcf = SourceParameters(100.0, 3.5, 2.75)
        @test srcd.Δσ.value == srcf.Δσ

        @test StochasticGroundMotionSimulation.magnitude_to_moment(6.0) == exp10(25.05)

    end

    @testset "Path" begin

        Rrefi = [ 1.0, 50.0, Inf ]
        γi = [ 1.0, 0.5 ]
        geof = GeometricSpreadingParameters(Rrefi, γi)
        geod = GeometricSpreadingParameters(Rrefi, [ 0.5 ], [ Dual{Float64}(1.0) ], BitVector([1,0]), :Piecewise)

        T = StochasticGroundMotionSimulation.get_parametric_type(geof)
        @test T == Float64
        T = StochasticGroundMotionSimulation.get_parametric_type(geod)
        @test T <: Dual

        sat = NearSourceSaturationParameters(:BT15)

        T = StochasticGroundMotionSimulation.get_parametric_type(sat)
        @test T == Float64

        Q0 = 200.0
        anef = AnelasticAttenuationParameters(Q0)
        aned = AnelasticAttenuationParameters(Dual{Float64}(Q0))

        T = StochasticGroundMotionSimulation.get_parametric_type(anef)
        @test T == Float64
        T = StochasticGroundMotionSimulation.get_parametric_type(aned)
        @test T <: Dual

        pathf = PathParameters(geof, sat, anef)
        pathd = PathParameters(geod, sat, aned)

        T = StochasticGroundMotionSimulation.get_parametric_type(pathf)
        @test T == Float64
        T = StochasticGroundMotionSimulation.get_parametric_type(pathd)
        @test T <: Dual

        r = 10.0
        m = 6.0

        geo_cy14 = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5], :CY14)
        geo_cy14mod = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5], :CY14mod)
        geo_null = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5], :Null)

        sat = NearSourceSaturationParameters(:BT15)
        r_ps = equivalent_point_source_distance(r, m, sat)

        gr_cy14 = geometric_spreading(r_ps, m, geo_cy14, sat)
        gr_cy14mod = geometric_spreading(r_ps, m, geo_cy14mod, sat)
        @test gr_cy14mod < gr_cy14
        @test isnan(geometric_spreading(r_ps, m, geo_null, sat))

        grp = geometric_spreading(r_ps, pathf)
        fasf = FourierParameters(SourceParameters(50.0), pathf)
        grf = geometric_spreading(r_ps, fasf)
        @test grp == grf

        geop = GeometricSpreadingParameters([1.0, Inf], [1.0], Vector{Float64}(), BitVector(undef,0), :Piecewise )
        @test StochasticGroundMotionSimulation.geometric_spreading_piecewise(r_ps, geop) == 1.0
        @test StochasticGroundMotionSimulation.geometric_spreading_piecewise(Dual(r_ps), geop) == 1.0
        fasp = FourierParameters(SourceParameters(50.0), PathParameters(geop, sat, anef))
        @test StochasticGroundMotionSimulation.geometric_spreading_piecewise(r_ps, fasp) == 1.0
        @test StochasticGroundMotionSimulation.geometric_spreading_piecewise(Dual(r_ps), fasp) == 1.0

        geo_cy14d = GeometricSpreadingParameters([1.0, 50.0, Inf], [Dual(1.0), Dual(1.0)], :CY14)
        sat = NearSourceSaturationParameters(:None)
        r_ps = equivalent_point_source_distance(1.0, -5.0, sat)
        fas_cy14d = FourierParameters(SourceParameters(50.0), PathParameters(geo_cy14d, sat, anef))

        @test StochasticGroundMotionSimulation.geometric_spreading_cy14(r_ps, fas_cy14d).value ≈ 1.0

        geo_cy14d = GeometricSpreadingParameters([1.0, 50.0, Inf], [Dual(1.0), Dual(1.0)], :CY14mod)
        sat = NearSourceSaturationParameters(:None)
        r_ps = equivalent_point_source_distance(1.0, -5.0, sat)
        fas_cy14d = FourierParameters(SourceParameters(50.0), PathParameters(geo_cy14d, sat, anef))

        @test StochasticGroundMotionSimulation.geometric_spreading_cy14mod(r_ps, -5.0, fas_cy14d).value ≈ 1.0

        # @code_warntype near_source_saturation(m, pathf.saturation)
        # @code_warntype near_source_saturation(m, pathd.saturation)
        hf = near_source_saturation(m, pathf.saturation)
        hd = near_source_saturation(m, pathd.saturation)
        if StochasticGroundMotionSimulation.get_parametric_type(pathd.saturation) <: Dual
            @test hf == hd.value
        else
            @test hf == hd
        end

        src = SourceParameters(100.0)
        geo = GeometricSpreadingParameters([1.0, Inf], [ 1.0 ])
        ane = AnelasticAttenuationParameters(200.0, 0.5)
        sat_ya = NearSourceSaturationParameters(:YA15)
        sat_cy = NearSourceSaturationParameters(:CY14)
        sat_none = NearSourceSaturationParameters(:None)
        sat_con = NearSourceSaturationParameters(5.0)
        sat_var = NearSourceSaturationParameters(Dual(5.0))
        sat_null = NearSourceSaturationParameters(:Null)

        path_ya = PathParameters(geo, sat_ya, ane)
        path_cy = PathParameters(geo, sat_cy, ane)
        path_none = PathParameters(geo, sat_none, ane)
        path_con = PathParameters(geo, sat_con, ane)
        path_var = PathParameters(geo, sat_var, ane)
        path_null = PathParameters(geo, sat_null, ane)

        h_ya = near_source_saturation(5.0, sat_ya)
        h_cy = near_source_saturation(5.0, sat_cy)
        h_none = near_source_saturation(5.0, sat_none)
        h_con = near_source_saturation(5.0, sat_con)
        h_var = near_source_saturation(5.0, sat_var)
        h_null = near_source_saturation(5.0, sat_null)

        # @code_warntype near_source_saturation(m, pathf)
        # @code_warntype near_source_saturation(m, pathd)
        hf = near_source_saturation(m, pathf)
        hd = near_source_saturation(m, pathd)
        if StochasticGroundMotionSimulation.get_parametric_type(pathd.saturation) <: Dual
            @test hf == hd.value
        else
            @test hf == hd
        end

        # @code_warntype equivalent_point_source_distance(r, m, pathf)
        # @code_warntype equivalent_point_source_distance(r, m, pathd)
        r_psf = equivalent_point_source_distance(r, m, pathf)
        r_psd = equivalent_point_source_distance(r, m, pathd)

        # @code_warntype StochasticGroundMotionSimulation.geometric_spreading_piecewise(r, geof)
        # @code_warntype StochasticGroundMotionSimulation.geometric_spreading_piecewise(r, geod)
        grf = StochasticGroundMotionSimulation.geometric_spreading_piecewise(r, geof)
        grd = StochasticGroundMotionSimulation.geometric_spreading_piecewise(r, geod)
        @test grf == grd.value

        # @code_warntype StochasticGroundMotionSimulation.geometric_spreading_cy14(r, geof)
        # @code_warntype StochasticGroundMotionSimulation.geometric_spreading_cy14(r, geod)
        grf = StochasticGroundMotionSimulation.geometric_spreading_cy14(r, geof)
        grd = StochasticGroundMotionSimulation.geometric_spreading_cy14(r, geod)
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


        geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [ 1.0, 0.5 ])
        geoc = GeometricSpreadingParameters([1.0, 50.0, Inf], [ 1.0, 0.5 ], :CY14)
        @test geo.model == :Piecewise
        @test geoc.model == :CY14
        geod = GeometricSpreadingParameters([1.0, 50.0, Inf], [Dual(1.0), Dual(0.5)])
        @test geo.γconi[1] == geod.γvari[1].value

        sat = NearSourceSaturationParameters(1, :BT15)
        sat = NearSourceSaturationParameters([5.5, 7.0, Inf], [4.0, 6.0], :ConstantConstrained)
        sat = NearSourceSaturationParameters([5.5, 7.0, Inf], [Dual(4.0), Dual(6.0)], :ConstantVariable)
        sat = NearSourceSaturationParameters([5.5, 7.0, Inf], [4.0, 6.0])
        sat = NearSourceSaturationParameters([5.5, 7.0, Inf], [Dual(4.0), Dual(6.0)])
        sat = NearSourceSaturationParameters(5.0)
        satd = NearSourceSaturationParameters(Dual(5.0))
        sat = NearSourceSaturationParameters(5.0, 2)
        sat = NearSourceSaturationParameters(Dual(5.0), 2)

        ane = AnelasticAttenuationParameters(200.0, 0.5, 3.5)
        @test ane.rmetric == :Rps

        path = PathParameters(geo, ane)
        @test path.saturation.model == :None

        path = PathParameters(geo, satd, ane)
        @test StochasticGroundMotionSimulation.get_parametric_type(path) <: Dual

        fas = FourierParameters(SourceParameters(100.0), path)
        @test fas.site.model == :Unit

        q_r = anelastic_attenuation(1.0, 10.0, fas)

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


        @test StochasticGroundMotionSimulation.boore_2016_generic_amplification(0.015) == 1.01
        @test isnan(StochasticGroundMotionSimulation.boore_2016_generic_amplification(NaN))

        numf = length(StochasticGroundMotionSimulation.fii_b16)
        for i in 1:numf
            fi = StochasticGroundMotionSimulation.fii_b16[i]
            @test StochasticGroundMotionSimulation.boore_2016_generic_amplification(fi) == StochasticGroundMotionSimulation.Aii_b16[i]
        end

        @test isnan(StochasticGroundMotionSimulation.alatik_2021_cy14_inverted_amplification_seg(NaN))

        numf = length(StochasticGroundMotionSimulation.fii_aa21_cy14)
        for i in 1:numf
            fi = StochasticGroundMotionSimulation.fii_aa21_cy14[i]
            @test StochasticGroundMotionSimulation.alatik_2021_cy14_inverted_amplification_seg(fi) == StochasticGroundMotionSimulation.Aii_aa21_cy14[i]
        end

        @test isnan(StochasticGroundMotionSimulation.alatik_2021_cy14_inverted_amplification(NaN))

        numf = length(StochasticGroundMotionSimulation.fii_aa21_cy14)
        for i in 1:numf
            fi = StochasticGroundMotionSimulation.fii_aa21_cy14[i]
            @test StochasticGroundMotionSimulation.alatik_2021_cy14_inverted_amplification(fi) == StochasticGroundMotionSimulation.Aii_aa21_cy14[i]
        end

        numf = length(StochasticGroundMotionSimulation.fii_b16)
        for i in 1:numf
            fi = StochasticGroundMotionSimulation.fii_b16[i]
            Afi_all = StochasticGroundMotionSimulation.alatik_2021_cy14_inverted_amplification(fi)
            Afi_seg = StochasticGroundMotionSimulation.alatik_2021_cy14_inverted_amplification_seg(fi)
            @test Afi_all ≈ Afi_seg
        end

    end

    @testset "Oscillator" begin
        ζ = 0.05
        f_n = 1.0
        sdof = Oscillator(f_n, ζ)

        @test f_n ≈ 1.0/period(sdof)

        @test transfer(0.5, sdof)^2 ≈ StochasticGroundMotionSimulation.squared_transfer(0.5, sdof)

        fi = [ 0.5, 1.0, 2.0 ]

        # @code_warntype transfer(fi, sdof)

        Hfi = transfer(fi, sdof)
        tfi = 2 * fi
        transfer!(Hfi, tfi, sdof)
        @test Hfi ≈ transfer(tfi, sdof)

        StochasticGroundMotionSimulation.squared_transfer!(Hfi, fi, sdof)
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
        # @code_warntype StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, srcf)
        Dsf = StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, srcf)
        @test Dsf ≈ 1.0/faf
        @test isnan(fbf)
        @test isnan(εf)
        fad, fbd, εd = corner_frequency(m, srcd)
        # @code_warntype StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, srcd)
        Dsd = StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, srcd)
        @test Dsd ≈ 1.0/fad
        @test isnan(fbd)
        @test isnan(εd)

        # @code_warntype StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, fasf)
        @test StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, fasf) ≈ Ds
        @test StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, fasd) ≈ Ds

        m = 6.0
        r = 7.0
        fa, fb, ε = corner_frequency(m, fasf)
        Dur = 1.0 / fa + 2.4
        @test StochasticGroundMotionSimulation.boore_thompson_2014(m, r, fasf) ≈ Dur
        fa, fb, ε = corner_frequency(m, fasd)
        Dur = 1.0 / fa + 2.4
        @test StochasticGroundMotionSimulation.boore_thompson_2014(m, r, fasd) ≈ Dur


        srcAS = SourceParameters(Δσf, :Atkinson_Silva_2000)
        fa, fb, ε = corner_frequency(m, srcAS)
        Ds = 0.5 * ( 1.0 / fa + 1.0 / fb )
        @test StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, srcAS) ≈ Ds
        @test isnan(fb) == false
        @test isnan(ε) == false
        @test fa < fb

        # test gradient of BT14 duration model w.r.t. magnitude
        h = 0.05
        m1 = 8.0
        m2 = m1 + h
        r_ps1 = 1.0 + near_source_saturation(m1, fasf)
        r_ps2 = 1.0 + near_source_saturation(m2, fasf)
        Dex1 = StochasticGroundMotionSimulation.boore_thompson_2014(m1, r_ps1, fasf)
        Dex2 = StochasticGroundMotionSimulation.boore_thompson_2014(m2, r_ps2, fasf)
        fdg = log(Dex2/Dex1)/h

        d(x) = log(StochasticGroundMotionSimulation.boore_thompson_2014(x[1], 1.0 + near_source_saturation(x[1], fasf), fasf))
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
        c11f = StochasticGroundMotionSimulation.StochasticGroundMotionSimulation.boore_thompson_2012_coefs(1, 1)
        @test c11f[1] == c11[1]
        @test all(isapprox.(c11, c11f))

        # @code_warntype StochasticGroundMotionSimulation.StochasticGroundMotionSimulation.boore_thompson_2012_coefs(1, 1)

        m = 8.0
        r = 1.0
        Dex = StochasticGroundMotionSimulation.boore_thompson_2014(m, r, srcf)
        # get the oscillator period
        sdof = Oscillator(100.0)
        T_n = period(sdof)
        ζ = sdof.ζ_n
        # define the η parameter as T_n/Dex
        η = T_n / Dex

        # @time c = StochasticGroundMotionSimulation.StochasticGroundMotionSimulation.boore_thompson_2012_coefs(1, 1)
        # @time StochasticGroundMotionSimulation.boore_thompson_2012_base(η, c, ζ)
        #
        # @time StochasticGroundMotionSimulation.boore_thompson_2012(m, r, srcf, sdof, rvt)
        # @code_warntype StochasticGroundMotionSimulation.boore_thompson_2012(m, r, srcf, sdof, rvt)
        #
        # @time StochasticGroundMotionSimulation.boore_thompson_2012(m, r, srcd, sdof, rvt)
        # @code_warntype StochasticGroundMotionSimulation.boore_thompson_2012(m, r, srcd, sdof, rvt)

        sdof = Oscillator(1.0)

        Drms, Dex, Dratio = StochasticGroundMotionSimulation.boore_thompson_2012(m, r, fas, sdof, rvt)
        Dex0 = StochasticGroundMotionSimulation.boore_thompson_2014(m, r, fas)

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

        # Boore & Thompson 2014
        m = 6.0
        fa, fb, ε = corner_frequency(m, src)
        Ds = 1.0 / fa

        Dex270 = Ds + 34.2
        Dex300 = Dex270 + 0.156*30.0
        @test Dex270 ≈ StochasticGroundMotionSimulation.boore_thompson_2014(m, 270.0, src)
        @test Dex300 ≈ StochasticGroundMotionSimulation.boore_thompson_2014(m, 300.0, src)

        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2014(m, -1.0, src))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2014(m, -1.0, srcd))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2014(m, Dual(-1.0), src))

        @test isnan(excitation_duration(m, -1.0, src, rvt))
        @test isnan(excitation_duration(m, -1.0, srcd, rvt))
        @test isnan(excitation_duration(m, Dual(-1.0), src, rvt))

        c = StochasticGroundMotionSimulation.StochasticGroundMotionSimulation.boore_thompson_2012_coefs(1, 1, region=:ENA)
        idx = 1
        @test all(isapprox(c, StochasticGroundMotionSimulation.coefs_ena_bt12[idx,3:9]))

        d1a, d2a, d3a = StochasticGroundMotionSimulation.boore_thompson_2012(6.1234, 2.0, src, sdof, rvt)
        d1b, d2b, d3b = StochasticGroundMotionSimulation.boore_thompson_2012(6.1234, 2.1234, src, sdof, rvt)
        d1c, d2c, d3c = StochasticGroundMotionSimulation.boore_thompson_2012(6.0, 2.1234, src, sdof, rvt)
        d1d, d2d, d3d = StochasticGroundMotionSimulation.boore_thompson_2012(6.0, 2.0, src, sdof, rvt)
        @test d1a < d1b
        @test d2a < d2b
        @test d3a > d3b
        @test d1a > d1c
        @test d1d < d1a


        c = StochasticGroundMotionSimulation.boore_thompson_2015_coefs(1, 1, region=:ENA)
        idx = 1
        @test all(isapprox(c, StochasticGroundMotionSimulation.coefs_ena_bt15[idx,3:9]))

        d1a, d2a, d3a = StochasticGroundMotionSimulation.boore_thompson_2015(6.1234, 2.0, src, sdof, rvt)
        d1b, d2b, d3b = StochasticGroundMotionSimulation.boore_thompson_2015(6.1234, 2.1234, src, sdof, rvt)
        d1c, d2c, d3c = StochasticGroundMotionSimulation.boore_thompson_2015(6.0, 2.1234, src, sdof, rvt)
        d1d, d2d, d3d = StochasticGroundMotionSimulation.boore_thompson_2015(6.0, 2.0, src, sdof, rvt)
        @test d1a < d1b
        @test d2a < d2b
        @test d3a > d3b
        @test d1a > d1c
        @test d1d < d1a

        d1as, d2as, d3as = StochasticGroundMotionSimulation.boore_thompson_2015(6.1234, 2.0, src, sdof, rvt)
        d1af, d2af, d3af = StochasticGroundMotionSimulation.boore_thompson_2015(6.1234, 2.0, fas, sdof, rvt)
        @test d1as == d1af
        @test d2as == d2af
        @test d3as == d3af

        rvt = RandomVibrationParameters(:PS, :PS, :PS, :PS)
        d1, d2, d3 = rms_duration(6.0, 1.0, fas, sdof, rvt)
        @test isnan(d1)
        @test isnan(d2)
        @test isnan(d3)

        srcf = SourceParameters(100.0)
        srcd = SourceParameters(Dual(100.0))
        rvt = RandomVibrationParameters(:PS, :PS, :PS, :PS)

        @test isnan(excitation_duration(6.0, 10.0, srcf, rvt))
        @test isnan(excitation_duration(6.0, 10.0, srcd, rvt))
        @test isnan(excitation_duration(6.0, Dual(10.0), srcf, rvt))

        D50 = excitation_duration(6.0, 50.0, srcf, RandomVibrationParameters())
        D150 = excitation_duration(6.0, 150.0, srcf, RandomVibrationParameters())
        @test D150 > D50

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

        Tf = StochasticGroundMotionSimulation.get_parametric_type(fasf)
        Td = StochasticGroundMotionSimulation.get_parametric_type(fasd)
        Tm = StochasticGroundMotionSimulation.get_parametric_type(fasm)
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

        src_a = SourceParameters(100.0, :Atkinson_Silva_2000)
        Af_a = fourier_source_shape(f, m, src_a)
        @test Af_a ≈ 1.0 atol=1e-3
        src_n = SourceParameters(100.0, :Null)
        Af_n = fourier_source_shape(f, m, src_n)
        src_b = SourceParameters(100.0)
        Af_b = fourier_source_shape(f, m, src_b)
        @test Af_n == Af_b

        fa, fb, ε = corner_frequency(m, src_a)
        Af_a = fourier_source_shape(f, fa, fb, ε, src_a.model)
        @test Af_a ≈ 1.0 atol=1e-3
        fa, fb, ε = corner_frequency(m, src_b)
        Af_n = fourier_source_shape(f, fa, fb, ε, src_n.model)
        @test Af_n ≈ Af_b


        # @code_warntype fourier_source(f, m, srcf)
        # @code_warntype fourier_source(f, m, srcd)
        Afs = fourier_source(f, m, srcf)
        Aff = fourier_source(f, m, fasf)
        @test Afs == Aff
        # @time fourier_source(f, m, srcd)
        # @time fourier_source(f, m, fasd)

        f = 10.0
        r = 100.0
        m = 6.0

        ane = AnelasticAttenuationParameters(200.0, 0.5, :Rrup)
        Pfr = fourier_path(f, r, m, geof, sat, ane)
        path = PathParameters(geof, sat, ane)
        Pfp = fourier_path(f, r, m, path)
        fas = FourierParameters(SourceParameters(100.0), path)
        Pff = fourier_path(f, r, m, fas)
        @test Pfr == Pfp
        @test Pfr == Pff


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

        @test fourier_attenuation(-1.0, r, fasf) == 1.0

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

        fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        fourier_spectrum!(Afim, fi, m, r_psd, fasm)
        for i in 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))

        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)
        for i in 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))



        ane = AnelasticAttenuationParameters(200.0, 0.0, :Rrup)
        path = PathParameters(geof, sat, ane)
        fas = FourierParameters(SourceParameters(100.0), path)

        m = Dual(6.0)
        r_rup = 10.0
        r_ps = equivalent_point_source_distance(r_rup, m, fas)

        Afid = fourier_spectrum(Vector{Float64}(), m, r_ps, fas)
        @test eltype(Afid) <: Dual
        Afid = fourier_spectrum(fi, m, r_ps, fas)

        sqAfid = StochasticGroundMotionSimulation.squared_fourier_spectrum(fi, m, r_ps, fas)
        @test any(isapprox.(sqAfid, Afid.^2))

        # @code_warntype fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        # @code_warntype fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        # @code_warntype fourier_spectrum!(Afim, fi, m, r_psd, fasm)


        Afid = fourier_spectrum(Vector{Float64}(), m, r_ps, fas)
        fourier_spectrum!(Afid, Vector{Float64}(), m, r_ps, fas)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, Vector{Float64}(), m, r_ps, fas)

        Afid = fourier_spectrum(fi, m, r_ps, fas)
        fourier_spectrum!(Afid, fi, m, r_ps, fas)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_ps, fas)

        # @code_warntype StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        # @code_warntype StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        # @code_warntype StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)


        # @code_warntype combined_kappa_frequency(r_psf, fasf)
        # @code_warntype combined_kappa_frequency(r_psd, fasd)
        fkf = combined_kappa_frequency(r_psf, 0.5, fasf)
        fkd = combined_kappa_frequency(r_psd, 0.5, fasd)
        @test fkf == fkd.value

        fkf = combined_kappa_frequency(r_psf, 0.5, fasf)
        fkfd = combined_kappa_frequency(Dual(r_psf), 0.5, fasf)
        @test fkf == fkfd.value

        fkf0 = combined_kappa_frequency(r_psf, 0.5, fas)
        fkf1 = combined_kappa_frequency(r_psf, 0.5, fasf)
        @test fkf0 > fkf1

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


            # Boore comparison (assume his are cgs units)
            # ps2db(f) = ( (2π * sdof.f_n) / ( (2π * f)^2 ) )^2 * 1e-4
            # ps2db(f) = ( 1.0 / ( 2π * f^2 * sdof.f_n ) )^2
            ps2db(f) = ( 1.0 / ( 2π*sdof.f_n ) )^2 * 1e4

            dbm0_integrand(f) = StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm0ln_integrand(lnf) = exp(lnf) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm1_integrand(f) = (2π*f) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm1ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf)) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm2_integrand(f) = (2π*f)^2 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm2ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^2  * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm4_integrand(f) = (2π*f)^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm4ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            @time igk = 2*quadgk(dbm0_integrand, 0.0, Inf)[1]
            @time igk = 2*quadgk(dbm0_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*StochasticGroundMotionSimulation.gauss_intervals(dbm0ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol=1e-4

            @time igk = 2*quadgk(dbm1_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*StochasticGroundMotionSimulation.gauss_intervals(dbm1ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol=1e-4

            @time igk = 2*quadgk(dbm2_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*StochasticGroundMotionSimulation.gauss_intervals(dbm2ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol=1e-3

            @time igk = 2*quadgk(dbm4_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2*StochasticGroundMotionSimulation.gauss_intervals(dbm4ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol=1e-3


            # using DifferentialEquations
            #
            # u0 = 0.0
            # tspan = (0.0, 300.0)
            # f0(u,p,t) = dbm0_integrand(t)
            # prob = ODEProblem(f0, u0, tspan)
            # sol = solve(prob, RK4())
            # 2*sol.u[end]


            m0_integrand(f) = StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m0ln_integrand(lnf) = exp(lnf) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)

            m1_integrand(f) = (2π*f) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m1ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf)) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)

            m2_integrand(f) = (2π*f)^2 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m2ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^2  * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)

            m4_integrand(f) = (2π*f)^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m4ln_integrand(lnf) = exp(lnf) * (2π*exp(lnf))^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)


            @time igk = quadgk(m0_integrand, exp(-7.0), exp(7.0))[1]
            @time igle = StochasticGroundMotionSimulation.gauss_interval(m0_integrand, 2000, 0.0, 300.0)
            @time igleln = StochasticGroundMotionSimulation.gauss_interval(m0ln_integrand, 750, -7.0, 7.0)

            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m0ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            # @test igk ≈ igle rtol=1e-2
            @test igk ≈ igleln rtol=1e-4
            @test igk ≈ iglelnm rtol=1e-4

            lnfi = log.([ 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, sdof.f_n ])
            sort!(lnfi)

            @time igk = quadgk(m1_integrand, exp(-7.0), exp(7.0))[1]
            @time igle = StochasticGroundMotionSimulation.gauss_interval(m1_integrand, 1500, 0.0, 300.0)
            @time igleln = StochasticGroundMotionSimulation.gauss_interval(m1ln_integrand, 750, -7.0, 7.0)
            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m1ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m1ln_integrand, 30, lnfi...)

            @time itr = StochasticGroundMotionSimulation.trapezoidal(m1ln_integrand, 60, lnfi...)

            # @test igk ≈ igle rtol=1e-2
            @test igk ≈ iglelnm rtol=1e-5
            @test igk ≈ itr rtol=1e-3


            @time igk = quadgk(m2_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m2ln_integrand, 30, lnfi...)
            @time itr = StochasticGroundMotionSimulation.trapezoidal(m2ln_integrand, 60, lnfi...)

            @test igk ≈ iglelnm rtol=1e-5
            @test igk ≈ itr rtol=1e-3

            @time igk = quadgk(m4_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m4ln_integrand, 30, lnfi...)
            @time itr = StochasticGroundMotionSimulation.trapezoidal(m4ln_integrand, 60, lnfi...)

            @test igk ≈ iglelnm rtol=1e-4
            @test igk ≈ itr rtol=1e-3




            integrand(x) = sin(x)

            intervals = 101
            x_min = 0.0
            x_max = 2π
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            ist = StochasticGroundMotionSimulation.trapezoidal_rule(xx, yy)

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

            integrand(f) = StochasticGroundMotionSimulation.squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = sdof.f_n/1.1
            x_max = sdof.f_n*1.1
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk rtol=1e-6
            @test isr ≈ igk atol=1e-6


            integrand(f) = StochasticGroundMotionSimulation.squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = 100.0
            x_max = 200.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk rtol=1e-3
            @test isr ≈ igk atol=1e-6


            integrand(f) = (2π*f)^4 * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = 100.0
            x_max = 200.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk rtol=1e-3
            @test isr ≈ igk atol=1e-6

            x_min = 300.0
            x_max = 500.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
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
            smigk = StochasticGroundMotionSimulation.spectral_moments_gk([0, 1, 2, 4], m, r_psf, fasf, sdof)

            @test all(isapprox.(smi, smigk, rtol=1e-3))


            sdof = Oscillator(1/3)
            smi = spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            smigk = StochasticGroundMotionSimulation.spectral_moments_gk([0, 1, 2, 4], m, r_psf, fasf, sdof)

            @test all(isapprox.(smi, smigk, rtol=1e-3))

            m = Dual(6.0)
            smdi = spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            smfi = spectral_moments([0, 1, 2, 4], m.value, r_psf, fasf, sdof)
            @test any(isapprox.(smdi, smfi))

            m = Dual(6.0)
            smdi = spectral_moment(0, m, r_psf, fasf, sdof)
            smfi = spectral_moment(0, m.value, r_psf, fasf, sdof)
            @test smdi ≈ smfi

            rvt = RandomVibrationParameters()
            Ti = [ 0.01, 0.1, 1.0 ]
            m = Dual(6.0)
            Sadi = rvt_response_spectrum(Ti, m, 10.0, fasf, rvt)
            Safi = rvt_response_spectrum(Ti, m.value, 10.0, fasf, rvt)
            for i in 1:length(Ti)
                @test Sadi[i].value ≈ Safi[i]
            end

            rvt_response_spectrum!(Sadi, Ti, m, 10.0, fasf, rvt)
            rvt_response_spectrum!(Safi, Ti, m.value, 10.0, fasf, rvt)
            for i in 1:length(Ti)
                @test Sadi[i].value ≈ Safi[i]
            end


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

            @time pfps = StochasticGroundMotionSimulation.peak_factor_dk80(m, r_psf, fasf, sdof)
            @time pfpsn = StochasticGroundMotionSimulation.peak_factor_dk80(m, r_psf, fasf, sdof, nodes=30)
            @time pfgk = StochasticGroundMotionSimulation.peak_factor_dk80_gk(m, r_psf, fasf, sdof)

            @test pfps ≈ pfgk rtol=1e-6
            @test pfpsn ≈ pfgk rtol=1e-6

            @time pfps = StochasticGroundMotionSimulation.peak_factor_cl56(m, r_psf, fasf, sdof)
            @time pfpsn = StochasticGroundMotionSimulation.peak_factor_cl56(m, r_psf, fasf, sdof, nodes=40)
            @time pfgk = StochasticGroundMotionSimulation.peak_factor_cl56_gk(m, r_psf, fasf, sdof)

            @test pfps ≈ pfgk rtol=1e-5
            @test pfpsn ≈ pfgk rtol=1e-5

            @test StochasticGroundMotionSimulation.vanmarcke_cdf(-1.0, 10.0, 0.5) == 0.0

            rvt = RandomVibrationParameters(:CL56)
            pf = peak_factor(6.0, 10.0, fasf, sdof, rvt)
            rvt = RandomVibrationParameters(:DK80)
            pf = peak_factor(6.0, 10.0, fasf, sdof, rvt)
            rvt = RandomVibrationParameters(:PS)
            pf = peak_factor(6.0, 10.0, fasf, sdof, rvt)
            @test isnan(pf)

            m = 6.0
            r_rup = 10.0
            r_ps = equivalent_point_source_distance(r_rup, m, fasf)
            Dex = excitation_duration(m, r_ps, fasf, rvt)
            m0 = spectral_moment(0, m, r_ps, fasf, sdof)
            rvt = RandomVibrationParameters(:PS)
            pf = peak_factor(m, r_ps, Dex, m0, fasf, sdof, rvt)
            @test isnan(pf)
            pf = peak_factor(Dual(m), r_ps, Dex, m0, fasf, sdof, rvt)
            @test isnan(pf)
            pf = peak_factor(m, r_ps, Dex, Dual(m0), fasf, sdof, rvt)
            @test isnan(pf)
            rvt = RandomVibrationParameters(:CL56)
            pf = peak_factor(m, r_ps, Dex, m0, fasf, sdof, rvt)
            @test isnan(pf) == false
            @test pf == peak_factor(m, r_ps, fasf, sdof, RandomVibrationParameters(:CL56))

            pfi = StochasticGroundMotionSimulation.peak_factor_integrand_cl56(0.0, 10.0, 10.0)
            @test pfi ≈ 1.0
            pfi = StochasticGroundMotionSimulation.peak_factor_integrand_cl56(Inf, 10.0, 10.0)
            @test pfi ≈ 0.0

            pf0 = StochasticGroundMotionSimulation.peak_factor_cl56(10.0, 10.0)
            pf1 = StochasticGroundMotionSimulation.peak_factor_cl56(10.0, 10.0, nodes=50)
            @test pf0 ≈ pf1

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