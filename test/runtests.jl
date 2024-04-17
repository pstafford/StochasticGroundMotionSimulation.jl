using StochasticGroundMotionSimulation
using Test
using ForwardDiff
using ForwardDiff: Dual
using FastGaussQuadrature
using QuadGK
using LinearAlgebra
using StaticArrays
# using Distributions
# using BenchmarkTools


@testset "StochasticGroundMotionSimulation.jl" begin

    @testset "Performance" begin

        # @testset "Allocation Testing" begin
        #     src = SourceParameters(100.0)
        #     geo = GeometricSpreadingParameters([1.0, Inf], [1.0], :Piecewise)
        #     sat = NearSourceSaturationParameters(:BT15)
        #     ane = AnelasticAttenuationParameters(180.0, 0.45, :Rps)
        #     path = PathParameters(geo, sat, ane)
        #     site = SiteParameters(0.03)
        #     fas = FourierParameters(src, path, site)
        #     rvt = RandomVibrationParameters()


        #     function run_sims(T, num_sims, fas, rvt)
        #         md = Uniform(2.0, 8.0)
        #         rd = Uniform(1.0, 100.0)
        #         mi = rand(md, num_sims)
        #         ri = rand(rd, num_sims)
        #         Sai = zeros(num_sims)
        #         for i in 1:num_sims
        #             Sai[i] = rvt_response_spectral_ordinate(T, mi[i], ri[i], fas, rvt)
        #         end
        #         return sum(Sai)
        #     end

        #     T = 0.123
        #     num_sims = 100
        #     @benchmark run_sims(T, num_sims, fas, rvt)
        #     num_sims = 1_000
        #     @benchmark run_sims(T, num_sims, fas, rvt)
        # end

        Ti = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0]
        m = 4.0 + π
        r = 500.0 + π

        src = SourceParameters(100.0)
        geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
        ane = AnelasticAttenuationParameters(200.0, 0.4)
        sat = NearSourceSaturationParameters(:BT15)
        path = PathParameters(geo, sat, ane)
        site = SiteParameters(0.039)

        fas = FourierParameters(src, path, site)
        rvt = RandomVibrationParameters(:DK80)


        sdof = Oscillator(1.0)
        # @code_warntype Oscillator(1.0)
        # @code_warntype period(sdof)
        # @code_warntype transfer(1.0, sdof)

        # @code_warntype rvt_response_spectral_ordinate(Ti[1], m, r, fas, rvt)
        # @code_warntype rvt_response_spectrum(Ti, m, r, fas, rvt)
        # @time Sai = rvt_response_spectrum(Ti, m, r, fas, rvt)
        # @profile Sai = rvt_response_spectrum(Ti, m, r, fas, rvt)

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

        @test StochasticGroundMotionSimulation.magnitude_to_moment(6.0) ≈ exp10(25.05)


        @testset "Beresnev (2019)" begin
            Δσ = 100.0
            n = 1.0
            srcb = SourceParameters(Δσ, n)
            srcω = SourceParameters(Δσ)

            f = 25.0
            m = 5.0
            Afb = fourier_source_shape(f, m, srcb)
            Afω = fourier_source_shape(f, m, srcω)

            @test Afb ≈ Afω

            srcb1p5 = SourceParameters(Δσ, 1.5)
            Afb1p5 = fourier_source_shape(f, m, srcb1p5)

            @test Afb > Afb1p5

            Δσd = Dual{Float64}(Δσ)
            nd = Dual{Float64}(n)

            srcdd = SourceParameters(Δσd, nd)
            srcdf = SourceParameters(Δσd, n)

            @test fourier_source_shape(f, m, srcdd) ≈ fourier_source_shape(f, m, srcb)
            @test fourier_source_shape(f, m, srcdf) ≈ fourier_source_shape(f, m, srcb)

        end

    end

    @testset "Path" begin

        Rrefi = [1.0, 50.0, Inf]
        γi = [1.0, 0.5]

        @testset "Geometric Spreading Constructors" begin
            # test floating point instantiation
            @test typeof(GeometricSpreadingParameters(Rrefi, γi)) <: GeometricSpreadingParameters
            # test Dual instantiation
            @test typeof(GeometricSpreadingParameters(Rrefi, [0.5], [Dual{Float64}(1.0)], BitVector([1, 0]), :Piecewise)) <: GeometricSpreadingParameters

            geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
            geoc = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5], :CY14)
            @test geo.model == :Piecewise
            @test geoc.model == :CY14
            geod = GeometricSpreadingParameters([1.0, 50.0, Inf], [Dual(1.0), Dual(0.5)])
            @test geo.γconi[1] == geod.γvari[1].value

            geo_mix = GeometricSpreadingParameters([1.0, Dual{Float64}(50.0), Inf], [1.0, 0.5])
            @test geo_mix.model == :Piecewise
        end

        geof = GeometricSpreadingParameters(Rrefi, γi)
        geod = GeometricSpreadingParameters(Rrefi, [0.5], [Dual{Float64}(1.0)], BitVector([1, 0]), :Piecewise)
        geom = GeometricSpreadingParameters([1.0, Dual{Float64}(50.0), Inf], [1.0, 0.5])

        @testset "Geometric Spreading Types" begin
            T = StochasticGroundMotionSimulation.get_parametric_type(geof)
            @test T == Float64
            T = StochasticGroundMotionSimulation.get_parametric_type(geod)
            @test T <: Dual
            T = StochasticGroundMotionSimulation.get_parametric_type(geom)
            @test T <: Dual
        end

        @testset "Near Source Saturation Constructors" begin
            # test standard instantiation with known models
            @test typeof(NearSourceSaturationParameters(:BT15)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters(:YA15)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters(:CY14)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters(:CY14mod)) <: NearSourceSaturationParameters

            @test typeof(NearSourceSaturationParameters(1, :BT15)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters([5.5, 7.0, Inf], [4.0, 6.0], :ConstantConstrained)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters([5.5, 7.0, Inf], [Dual(4.0), Dual(6.0)], :ConstantVariable)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters([5.5, 7.0, Inf], [4.0, 6.0])) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters([5.5, 7.0, Inf], [Dual(4.0), Dual(6.0)])) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters(5.0)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters(Dual(5.0))) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters(5.0, 2)) <: NearSourceSaturationParameters
            @test typeof(NearSourceSaturationParameters(Dual(5.0), 2)) <: NearSourceSaturationParameters

        end

        sat = NearSourceSaturationParameters(:BT15)
        satd = NearSourceSaturationParameters(Dual{Float64}(5.0))

        @testset "Near Source Saturation Types" begin
            T = StochasticGroundMotionSimulation.get_parametric_type(sat)
            @test T == Float64
        end


        @testset "Anelastic Attenuation Constructors" begin
            # test floating point instantiation
            @test typeof(AnelasticAttenuationParameters(200.0)) <: AnelasticAttenuationParameters
            # test Dual instantiation
            @test typeof(AnelasticAttenuationParameters(Dual{Float64}(200.0))) <: AnelasticAttenuationParameters

            ane = AnelasticAttenuationParameters(200.0, 0.5, 3.5)
            @test ane.rmetric == :Rps
            ane = AnelasticAttenuationParameters(200.0, 0.5, 3.5, :Rrup)
            @test ane.rmetric == :Rrup

            # test the segmented versions
            # scalar inputs (internally mapped to vectors)
            @test typeof(AnelasticAttenuationParameters(200.0, 0.5, :Rrup)) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters(200.0, Dual{Float64}(0.5), :Rrup)) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters(Dual{Float64}(200.0), 0.5, :Rrup)) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters(Dual{Float64}(200.0), Dual{Float64}(0.5), :Rrup)) <: AnelasticAttenuationParameters
            # vector inputs
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0])) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0], :Rrup)) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [Dual{Float64}(200.0), Dual{Float64}(200.0)])) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [Dual{Float64}(200.0), Dual{Float64}(200.0)], :Rrup)) <: AnelasticAttenuationParameters

            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0], [0.5, 0.5])) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0], [0.5, 0.5], :Rrup)) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [Dual{Float64}(200.0), Dual{Float64}(200.0)], [Dual{Float64}(0.5), Dual{Float64}(0.5)])) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [Dual{Float64}(200.0), Dual{Float64}(200.0)], [Dual{Float64}(0.5), Dual{Float64}(0.5)], :Rrup)) <: AnelasticAttenuationParameters

            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0], [Dual{Float64}(0.5), Dual{Float64}(0.5)])) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0], [Dual{Float64}(0.5), Dual{Float64}(0.5)], :Rrup)) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [Dual{Float64}(200.0), Dual{Float64}(200.0)], [0.5, 0.5])) <: AnelasticAttenuationParameters
            @test typeof(AnelasticAttenuationParameters([0.0, 80.0, Inf], [Dual{Float64}(200.0), Dual{Float64}(200.0)], [0.5, 0.5], :Rrup)) <: AnelasticAttenuationParameters

        end

        Q0 = 200.0
        anef = AnelasticAttenuationParameters(Q0)
        aned = AnelasticAttenuationParameters(Dual{Float64}(Q0))

        @testset "Anelastic Attenuation Types" begin
            T = StochasticGroundMotionSimulation.get_parametric_type(anef)
            @test T == Float64
            T = StochasticGroundMotionSimulation.get_parametric_type(aned)
            @test T <: Dual

            # artificial test, but for ensuring complete coverage
            anet = AnelasticAttenuationParameters([Dual{Float64}(0.0), Dual{Float64}(Inf)], [Dual{Float64}(200.0)], Vector{Dual{Float64,Float64,0}}(), zeros(Dual{Float64,Float64,0}, 1), zeros(Dual{Float64,Float64,0}, 1), [3.5], BitVector(ones(1)), BitVector(ones(1)), :Rps)
            T = StochasticGroundMotionSimulation.get_parametric_type(anet)
            @test T <: Dual
        end

        @testset "Path Constructors" begin
            # test floating point components
            @test typeof(PathParameters(geof, sat, anef)) <: PathParameters
            # test Dual components
            @test typeof(PathParameters(geod, sat, aned)) <: PathParameters

            path = PathParameters(geof, anef)
            @test path.saturation.model == :None

            path = PathParameters(geof, satd, anef)
            @test StochasticGroundMotionSimulation.get_parametric_type(path) <: Dual

        end

        pathf = PathParameters(geof, sat, anef)
        pathd = PathParameters(geod, sat, aned)

        @testset "Path Types" begin
            T = StochasticGroundMotionSimulation.get_parametric_type(pathf)
            @test T == Float64
            T = StochasticGroundMotionSimulation.get_parametric_type(pathd)
            @test T <: Dual

        end

        r = 10.0
        m = 6.0

        @testset "Geometric Spreading Functionality" begin

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

            geop = GeometricSpreadingParameters([1.0, Inf], [1.0], Vector{Float64}(), BitVector(undef, 0), :Piecewise)
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
        end

        @testset "Near Source Saturation Functionality" begin

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
            geo = GeometricSpreadingParameters([1.0, Inf], [1.0])
            ane = AnelasticAttenuationParameters(200.0, 0.5)
            sat_ya = NearSourceSaturationParameters(:YA15)
            sat_cy = NearSourceSaturationParameters(:CY14)
            sat_sea = NearSourceSaturationParameters(:SEA21)
            sat_none = NearSourceSaturationParameters(:None)
            sat_con = NearSourceSaturationParameters(5.0)
            sat_var = NearSourceSaturationParameters(Dual(5.0))
            sat_null = NearSourceSaturationParameters(:Null)

            path_ya = PathParameters(geo, sat_ya, ane)
            path_cy = PathParameters(geo, sat_cy, ane)
            path_sea = PathParameters(geo, sat_sea, ane)
            path_none = PathParameters(geo, sat_none, ane)
            path_con = PathParameters(geo, sat_con, ane)
            path_var = PathParameters(geo, sat_var, ane)
            path_null = PathParameters(geo, sat_null, ane)

            h_ya = near_source_saturation(5.0, sat_ya)
            h_cy = near_source_saturation(5.0, sat_cy)
            h_sea = near_source_saturation(5.0, sat_sea)
            h_none = near_source_saturation(5.0, sat_none)
            h_con = near_source_saturation(5.0, sat_con)
            h_var = near_source_saturation(5.0, sat_var)
            h_null = near_source_saturation(5.0, sat_null)

            @test h_con == h_var.value
        end


        f = 1.0
        r = 100.0

        @testset "Anelastic Attenuation Functionality" begin
            # @code_warntype anelastic_attenuation(f, r, anef)
            # @code_warntype anelastic_attenuation(f, r, aned)
            qrf = anelastic_attenuation(f, r, anef)
            qrd = anelastic_attenuation(f, r, aned)
            @test qrf ≈ qrd.value

            fas = FourierParameters(SourceParameters(100.0), pathf)
            q_r = anelastic_attenuation(f, r, fas)
            @test qrf ≈ q_r

            f = [0.01, 0.1, 1.0, 10.0, 100.0]
            nf = length(f)
            qrf = anelastic_attenuation(f, r, anef)
            qrd = anelastic_attenuation(f, r, aned)
            @test qrf ≈ map(q -> q.value, qrd)

            fasf = FourierParameters(SourceParameters(100.0), pathf)
            fasd = FourierParameters(SourceParameters(100.0), pathd)
            q_r = anelastic_attenuation(f, r, fas)
            @test qrf ≈ q_r

            Aff = ones(eltype(qrf), nf)
            Afd = ones(eltype(qrd), nf)
            StochasticGroundMotionSimulation.apply_anelastic_attenuation!(Aff, f, r, anef)
            StochasticGroundMotionSimulation.apply_anelastic_attenuation!(Afd, f, r, aned)
            @test Aff ≈ map(a -> a.value, Afd)

            Aff = ones(eltype(qrf), nf)
            Afd = ones(eltype(qrd), nf)
            StochasticGroundMotionSimulation.apply_anelastic_attenuation!(Aff, f, r, fasf)
            StochasticGroundMotionSimulation.apply_anelastic_attenuation!(Afd, f, r, fasd)
            @test Aff ≈ map(a -> a.value, Afd)

            Qff = ones(nf)
            Qfd = ones(eltype(qrd), nf)
            StochasticGroundMotionSimulation.anelastic_attenuation!(Qff, f, r, anef)
            StochasticGroundMotionSimulation.anelastic_attenuation!(Qfd, f, r, aned)
            @test Qff ≈ map(a -> a.value, Qfd)

        end

        @testset "Anelastic Attenuation Segmentation" begin
            # test the vector descriptions of anelastic attenuation
            ane_vec = AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0], [0.4, 0.4])
            ane_con = AnelasticAttenuationParameters(200.0, 0.4)

            @test anelastic_attenuation(5.0, 50.0, ane_vec) ≈ anelastic_attenuation(5.0, 50.0, ane_con)
            @test anelastic_attenuation(5.0, 150.0, ane_vec) ≈ anelastic_attenuation(5.0, 150.0, ane_con)

            ane_vecd = AnelasticAttenuationParameters([0.0, 80.0, Inf], [Dual(200.0), Dual(200.0)], [0.4, 0.4])
            ane_cond = AnelasticAttenuationParameters(Dual(200.0), 0.4)

            @test anelastic_attenuation(5.0, 50.0, ane_vecd) ≈ anelastic_attenuation(5.0, 50.0, ane_cond)
            @test anelastic_attenuation(5.0, 150.0, ane_vecd) ≈ anelastic_attenuation(5.0, 150.0, ane_cond)

            ane_vecd = AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0, 200.0], [Dual(0.4), Dual(0.4)])
            ane_cond = AnelasticAttenuationParameters(200.0, Dual(0.4))

            @test anelastic_attenuation(5.0, 50.0, ane_vecd) ≈ anelastic_attenuation(5.0, 50.0, ane_cond)
            @test anelastic_attenuation(5.0, 150.0, ane_vecd) ≈ anelastic_attenuation(5.0, 150.0, ane_cond)

            ane_vecd = AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0], [Dual(200.0)], [0.4], [Dual(0.4)], 3.5 * ones(2), BitVector([1, 0]), BitVector([0, 1]), :Rrup)
            ane_cond = AnelasticAttenuationParameters(200.0, 0.4)

            @test anelastic_attenuation(5.0, 50.0, ane_vecd) ≈ anelastic_attenuation(5.0, 50.0, ane_cond)
            @test anelastic_attenuation(5.0, 150.0, ane_vecd) ≈ anelastic_attenuation(5.0, 150.0, ane_cond)

            ane_vecd = AnelasticAttenuationParameters([0.0, 80.0, Inf], [200.0], [Dual(200.0)], [0.4], [Dual(0.4)], 3.5 * ones(2), BitVector([0, 1]), BitVector([1, 0]), :Rrup)
            ane_cond = AnelasticAttenuationParameters(200.0, 0.4)

            @test anelastic_attenuation(5.0, 50.0, ane_vecd) ≈ anelastic_attenuation(5.0, 50.0, ane_cond)
            @test anelastic_attenuation(5.0, 150.0, ane_vecd) ≈ anelastic_attenuation(5.0, 150.0, ane_cond)

            ane_inf = AnelasticAttenuationParameters([0.0, 80.0, Inf], [Inf], [Dual{Float64}(200.0)], [0.0], [Dual{Float64}(0.5)], 3.5 * ones(2), BitVector([0, 1]), BitVector([0, 1]), :Rrup)
            @test anelastic_attenuation(5.0, 50.0, ane_inf) == 1.0

            f_vec = [0.1, 1.0, 10.0, 100.0]
            nf = length(f_vec)
            q_vec = anelastic_attenuation(f_vec, 200.0, ane_vec)
            q_con = anelastic_attenuation(f_vec, 200.0, ane_con)
            @test q_vec ≈ q_con

            Afv = ones(nf)
            Afc = ones(nf)
            ζ0f = 0.039
            ηf = 0.75
            site = SiteParameters(ζ0f, ηf)
            StochasticGroundMotionSimulation.apply_fourier_path_and_site_attenuation!(Afv, f_vec, 200.0, ane_vec, site)
            StochasticGroundMotionSimulation.apply_fourier_path_and_site_attenuation!(Afc, f_vec, 200.0, ane_con, site)
            @test Afv ≈ Afc

            Afv = ones(nf)
            Afc = ones(nf)
            StochasticGroundMotionSimulation.apply_anelastic_attenuation!(Afv, f_vec, 200.0, ane_vec)
            StochasticGroundMotionSimulation.apply_anelastic_attenuation!(Afc, f_vec, 200.0, ane_con)
            @test Afv ≈ Afc

            Kfv = ones(nf)
            Kfc = ones(nf)
            StochasticGroundMotionSimulation.anelastic_attenuation!(Kfv, f_vec, 200.0, ane_vec)
            StochasticGroundMotionSimulation.anelastic_attenuation!(Kfc, f_vec, 200.0, ane_con)
            @test Kfv ≈ Kfc
            @test Kfv ≈ Afv

            Kfv = anelastic_attenuation(f_vec, 200.0, ane_vecd)
            Kfc = anelastic_attenuation(f_vec, 200.0, ane_cond)
            StochasticGroundMotionSimulation.anelastic_attenuation!(Kfv, f_vec, 200.0, ane_vecd)
            StochasticGroundMotionSimulation.anelastic_attenuation!(Kfc, f_vec, 200.0, ane_cond)
            @test Kfv ≈ Kfc
            @test Kfv ≈ Afv

        end

    end

    @testset "Site" begin

        κ0f = 0.039
        κ0d = Dual{Float64}(κ0f)
        ζ0f = 0.039
        ζ0d = Dual{Float64}(ζ0f)
        ηf = 0.75
        ηd = Dual{Float64}(ηf)

        @testset "Site Constructors" begin
            @test typeof(SiteAmpUnit()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpConstant(2.0)) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpBoore2016_760()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_ask14_620()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_ask14_760()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_ask14_1100()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_bssa14_620()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_bssa14_760()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_bssa14_1100()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_cb14_620()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_cb14_760()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_cb14_1100()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_cy14_620()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_cy14_760()) <: StochasticGroundMotionSimulation.SiteAmplification
            @test typeof(SiteAmpAlAtikAbrahamson2021_cy14_1100()) <: StochasticGroundMotionSimulation.SiteAmplification

            @test typeof(SiteParameters(κ0f)) <: SiteParameters

            @test typeof(SiteParameters(κ0f, SiteAmpAlAtikAbrahamson2021_cy14_760())) <: SiteParameters
            @test typeof(SiteParameters(κ0f, SiteAmpBoore2016_760())) <: SiteParameters
            @test typeof(SiteParameters(κ0f, SiteAmpUnit())) <: SiteParameters
            @test typeof(SiteParameters(κ0f, SiteAmpConstant(2.0))) <: SiteParameters

            @test typeof(SiteParameters(κ0d)) <: SiteParameters
            @test typeof(SiteParameters(κ0d, SiteAmpAlAtikAbrahamson2021_cy14_760())) <: SiteParameters
            @test typeof(SiteParameters(κ0d, SiteAmpBoore2016_760())) <: SiteParameters
            @test typeof(SiteParameters(κ0d, SiteAmpUnit())) <: SiteParameters

            @test typeof(SiteParameters(ζ0f, ηf)) <: SiteParameters
            @test typeof(SiteParameters(ζ0d, ηf)) <: SiteParameters
            @test typeof(SiteParameters(ζ0f, ηd)) <: SiteParameters
            @test typeof(SiteParameters(ζ0d, ηd)) <: SiteParameters
            @test typeof(SiteParameters(ζ0f, ηf, SiteAmpAlAtikAbrahamson2021_cy14_760())) <: SiteParameters
            @test typeof(SiteParameters(ζ0d, ηf, SiteAmpAlAtikAbrahamson2021_cy14_760())) <: SiteParameters
            @test typeof(SiteParameters(ζ0f, ηd, SiteAmpAlAtikAbrahamson2021_cy14_760())) <: SiteParameters
            @test typeof(SiteParameters(ζ0d, ηd, SiteAmpAlAtikAbrahamson2021_cy14_760())) <: SiteParameters

        end

        site0f = SiteParameters(κ0f)
        siteAf = SiteParameters(κ0f, SiteAmpAlAtikAbrahamson2021_cy14_760())
        siteBf = SiteParameters(κ0f, SiteAmpBoore2016_760())
        siteUf = SiteParameters(κ0f, SiteAmpUnit())
        siteCf = SiteParameters(κ0f, SiteAmpConstant(2.0))

        site0d = SiteParameters(κ0d)
        siteAd = SiteParameters(κ0d, SiteAmpAlAtikAbrahamson2021_cy14_760())
        siteBd = SiteParameters(κ0d, SiteAmpBoore2016_760())
        siteUd = SiteParameters(κ0d, SiteAmpUnit())
        siteCd = SiteParameters(κ0d, SiteAmpConstant(2.0))

        f = 0.05

        @testset "Site Amplification" begin
            # @code_warntype site_amplification(f, site0f)
            # @code_warntype site_amplification(f, site0d)
            Sff = site_amplification(f, site0f)
            Sfd = site_amplification(f, site0d)
            @test Sff == Sfd


            Af0f = site_amplification(f, site0f)
            Af1f = site_amplification(f, siteAf)
            Af2f = site_amplification(f, siteBf)
            Af3f = site_amplification(f, siteUf)
            Af4f = site_amplification(f, siteCf)

            @test Af0f == Af1f
            @test Af3f == 1.0
            @test Af2f < Af1f
            @test Af4f == 2.0

            Af0d = site_amplification(f, site0d)
            Af1d = site_amplification(f, siteAd)
            Af2d = site_amplification(f, siteBd)
            Af3d = site_amplification(f, siteUd)
            Af4d = site_amplification(f, siteCd)

            @test Af0d == Af1d
            @test Af3d == 1.0
            @test Af2d < Af1d
            @test Af4d == 2.0

            @test Af0f == Af0d
            @test Af4f == Af4d

            f0 = 80.0
            f1 = 100.0
            Af0 = site_amplification(f0, siteBf)
            Af1 = site_amplification(f1, siteBf)

            @test Af0 ≈ Af1 atol=1e-2

            ft = 10.0
            Af620 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_ask14_620()))
            Af760 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_ask14_760()))
            Af1100 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_ask14_1100()))

            @test Af620 > Af760
            @test Af760 > Af1100

            Af620 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_bssa14_620()))
            Af760 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_bssa14_760()))
            Af1100 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_bssa14_1100()))

            @test Af620 > Af760
            @test Af760 > Af1100

            Af620 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_cb14_620()))
            Af760 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_cb14_760()))
            Af1100 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_cb14_1100()))

            @test Af620 > Af760
            @test Af760 > Af1100

            Af620 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_cy14_620()))
            Af760 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_cy14_760()))
            Af1100 = site_amplification(ft, SiteParameters(0.039, SiteAmpAlAtikAbrahamson2021_cy14_1100()))

            @test Af620 > Af760
            @test Af760 > Af1100
        end


        @testset "Impedance Functions" begin
            sa = SiteAmpBoore2016_760()
            @test sa.amplification(0.015) ≈ 1.01 rtol = 1e-5

            fi = @SVector [0.001, 0.010, 0.015, 0.021, 0.031, 0.045, 0.065, 0.095, 0.138, 0.200, 0.291, 0.423, 0.615, 0.894, 1.301, 1.892, 2.751, 4.000, 5.817, 8.459, 12.301, 17.889, 26.014, 37.830, 55.012, 80.000, 1e3]
            Ai = @SVector [1.00, 1.00, 1.01, 1.02, 1.02, 1.04, 1.06, 1.09, 1.13, 1.18, 1.25, 1.32, 1.41, 1.51, 1.64, 1.80, 1.99, 2.18, 2.38, 2.56, 2.75, 2.95, 3.17, 3.42, 3.68, 3.96, 3.96]

            numf = length(fi) - 1
            for i = 1:numf
                f = fi[i]
                @test sa.amplification(f) ≈ Ai[i] rtol = 1e-2
            end


            sa = SiteAmpAlAtikAbrahamson2021_cy14_760()

            fi = @SVector [0.001, 0.100000005278119, 0.102329305685959, 0.104712901088038, 0.10715191734136, 0.1096478161355, 0.112201844312325, 0.114815399493054, 0.117489811269133, 0.120226426884127, 0.123026912673429, 0.125892543180806, 0.128824943888933, 0.131825701931171, 0.134896292542321, 0.138038419685919, 0.141253726031804, 0.144543991682274, 0.147910841016824, 0.15135612531762, 0.154881699105396, 0.158489298275321, 0.162180994731809, 0.165958688177051, 0.169824406249264, 0.173780124732544, 0.177827987234069, 0.181970104503309, 0.18620874216892, 0.190546088276735, 0.194984413397179, 0.199526234462532, 0.204173818378493, 0.208929629962364, 0.213796203501538, 0.218776216272469, 0.223872092751581, 0.229086773295876, 0.234422913440041, 0.239883346374602, 0.245470959191782, 0.251188647299664, 0.257039583965111, 0.263026839625322, 0.269153386068061, 0.275422892536035, 0.281838357206786, 0.288403136014915, 0.29512088643951, 0.301995147118928, 0.309029636270837, 0.316227726764557, 0.323593743035287, 0.331131074274485, 0.338844138033324, 0.346736738117816, 0.354813374379999, 0.363078069427477, 0.371535318437207, 0.38018929961147, 0.389044995694541, 0.398107060813636, 0.407380379735427, 0.416869369537323, 0.426579727361764, 0.436515887692042, 0.446683606986212, 0.457088019319652, 0.46773532542682, 0.478629981874554, 0.489778948828998, 0.501187304961605, 0.512861186273828, 0.524807402402299, 0.537031596210482, 0.549541299892976, 0.562341327015613, 0.575439822942147, 0.588843516357273, 0.602559763585899, 0.61659488244124, 0.630957589108888, 0.645654602589703, 0.660693172670039, 0.676083352446191, 0.69183061540227, 0.707945609900711, 0.724435639005472, 0.741310319666438, 0.758577791401586, 0.776247032532349, 0.794327972885161, 0.812829888313304, 0.83176450189098, 0.85113745411528, 0.870963690655165, 0.891251760400126, 0.912010408577856, 0.933255245295474, 0.954991724274509, 0.977236931585538, 0.99999931915226, 1.0232926583713, 1.04712871370939, 1.07151972559382, 1.09647839861261, 1.12201798920987, 1.14815239015486, 1.17489619113164, 1.20226482282661, 1.230270509361, 1.25892568257472, 1.28825221490343, 1.31825970420816, 1.34896234511787, 1.38038569480746, 1.4125355875376, 1.44544035634817, 1.4791067005516, 1.51355991769356, 1.54881419338307, 1.58489135570109, 1.62180741098713, 1.65958724224312, 1.69824098493074, 1.73780538614063, 1.77828384029801, 1.81969768325487, 1.86209132642581, 1.90545951786295, 1.94984000723927, 1.99526323502638, 2.04173666608573, 2.08929373051144, 2.13795720949869, 2.18776679563686, 2.23871864940605, 2.29087274897834, 2.34422665758295, 2.398830455509, 2.45470143491059, 2.51187774934974, 2.57040187591233, 2.63027473614656, 2.69154228983537, 2.75423051595883, 2.81839385215467, 2.88403692114349, 2.95122086420034, 3.01995215562392, 3.09030065678915, 3.16227642295472, 3.23592346631478, 3.31132891490246, 3.3884299172534, 3.46735946910046, 3.54813821074709, 3.63078533725135, 3.71536846377799, 3.80191135853523, 3.89043708383884, 3.98108381905827, 4.07382969500055, 4.16870567544569, 4.26581070164795, 4.36518697265844, 4.46680411703705, 4.57085377834242, 4.67739610111622, 4.78632242414544, 4.89776678860413, 5.0118843229469, 5.12865332406043, 5.24803935339088, 5.37032511094351, 5.49537847413284, 5.62339713074269, 5.75435457691296, 5.88847905118778, 6.02561155450687, 6.16599641750948, 6.30960714497232, 6.45659960120723, 6.60689861529496, 6.76084973072587, 6.91824357308313, 7.07938541613807, 7.24442301919855, 7.41307142632331, 7.58568641738688, 7.76244066564601, 7.94326569146627, 8.12834448596361, 8.31759182767968, 8.51151566867657, 8.70973318814878, 8.91245786665494, 9.12030564789694, 9.33246351716575, 9.54996668260048, 9.77234156401459, 9.99990153906584, 10.2330151684908, 10.4711008559889, 10.7150364469821, 10.964750490257, 11.2201347329491, 11.4816670597186, 11.7492496241653, 12.022742813086, 12.3026981000034, 12.5889920759082, 12.8822769465087, 13.182459265256, 13.489400932834, 13.8038791736706, 14.1258004265005, 14.4539459356652, 14.7913578984967, 15.1357324856635, 15.4880193842308, 15.8493619845298, 16.2183086824161, 16.5960326767835, 16.9823769767776, 17.3787778053556, 17.7834601578259, 18.1979135417807, 18.6200418512582, 19.0554824815177, 19.49807748989, 19.9517012657342, 20.4186610376656, 20.8940733113534, 21.379997361205, 21.8789412825561, 22.387937087537, 22.9096392147486, 23.4439269263859, 23.9870856553391, 24.46058662679, 25.0340588613218, 25.6152039490899, 26.2158596298214, 26.8275446219964, 27.4542925051998, 28.0959151088508, 28.7521309339848, 29.4225501552548, 30.1066616578723, 30.8101318295223, 31.5331607803149, 32.2688473500071, 33.0235676906087, 33.7894041350745, 34.5814827413102, 35.3920282683592, 36.2114034095699, 37.0570384215831, 37.9296731960298, 38.808356830246, 39.7246511532133, 40.6446342457211, 41.5904771517877, 42.5624316552899, 43.5606664842275, 44.58525356566, 45.6203356517073, 46.679926731282, 47.7811828082431, 48.889095202211, 50.0392749799897, 51.1923784631563, 52.3876989427786, 53.6271327545174, 54.8639723203735, 56.1438945185322, 57.4685744162854, 58.8111011775592, 60.1686734628619, 61.6017896398734, 63.0158092367494, 64.5080917430415, 66.010239337733, 67.5570266517399, 69.1074610525002, 70.7437395379651, 72.3800126937281, 74.0584033514694, 75.8305646487854, 77.5949512656221, 79.4004309784528, 81.2461188644527, 83.1308139099541, 85.0529665055861, 87.0822300450569, 89.0763954752011, 91.1807407221922, 93.3202913166248, 95.4916631049442, 97.6908631563263, 100.012267904509, 1e3]
            Ai = @SVector [1.0, 1.26474365754693, 1.27187168763244, 1.27901461204394, 1.2861668468288, 1.29332769860997, 1.30049545713136, 1.30766852765824, 1.3148448358383, 1.32202240367609, 1.32920111176275, 1.33638102751286, 1.34356422475535, 1.35075365610266, 1.35795288790118, 1.36516683588672, 1.37240078818257, 1.37965956984633, 1.38694651866824, 1.3942644665929, 1.40161534480099, 1.40900004959024, 1.41641948596057, 1.42387361692424, 1.43136209798135, 1.43888401061133, 1.44643774756037, 1.45402034276764, 1.46162798960857, 1.46925573360015, 1.47689789847254, 1.48454848588452, 1.49220043010574, 1.49984643882526, 1.50747891120157, 1.51509009761108, 1.52267156234887, 1.53021542332592, 1.53771329159611, 1.54515767577195, 1.5525423017734, 1.55986193063354, 1.56711260297083, 1.57429107685139, 1.58139470109949, 1.58842230933339, 1.59537254872675, 1.60224497526539, 1.6090399193016, 1.61575793863017, 1.62240013364471, 1.62896756722738, 1.63546241542315, 1.64188624895108, 1.64824169962538, 1.65453105215902, 1.660757253143, 1.66692298962119, 1.67303139815318, 1.67908545861599, 1.68508877287354, 1.69104474095905, 1.69695691445104, 1.70282847890391, 1.70866337029487, 1.71446466460501, 1.72023619436562, 1.72598137429174, 1.73170412937894, 1.73740741092509, 1.74309547823134, 1.74877153446106, 1.7544392973187, 1.76010278844154, 1.7657654338474, 1.77143153644529, 1.77710410272141, 1.7827876649385, 1.78848607318854, 1.79420347226643, 1.79994358407784, 1.80571108575488, 1.81150988535825, 1.81734412132426, 1.82321904347219, 1.82913817252675, 1.83510670026067, 1.84112756121076, 1.84720347424249, 1.85333562714735, 1.85952470296974, 1.8657708979849, 1.87207351949043, 1.87843213317722, 1.8848436987637, 1.89130781117208, 1.89782159094484, 1.90438189213204, 1.91098741512956, 1.91763345905683, 1.92431870117879, 1.93103904302077, 1.9377917893536, 1.94457361968632, 1.95138121657727, 1.95821137709387, 1.96506101732005, 1.97192719683423, 1.97880718423323, 1.98569844861012, 1.99259770558709, 1.99950182429154, 2.00641002878188, 2.01331868946611, 2.02022543975041, 2.02713050388327, 2.03402982723964, 2.04092437374658, 2.04781056500155, 2.05468885904031, 2.0615574700588, 2.06841626488153, 2.0752640907423, 2.08210150758018, 2.08892649642953, 2.09574183898203, 2.10254459906843, 2.10933506179724, 2.11611724811337, 2.12288728136106, 2.12964821625687, 2.13640189502178, 2.14314680483325, 2.14988534511825, 2.15661842855298, 2.16334930355606, 2.17007541648867, 2.17680275350116, 2.18352920748937, 2.19025940643645, 2.19699395722506, 2.2037360384332, 2.21048929651311, 2.21725267854449, 2.22403035334881, 2.23082429847101, 2.23763960582685, 2.24447608988636, 2.25133959947153, 2.25823040581599, 2.26515517696813, 2.27211465284156, 2.279113030692, 2.28615863788022, 2.29324589179053, 2.30038750095143, 2.30758593730311, 2.31484372550238, 2.32216775190769, 2.32956116012962, 2.33702714372695, 2.34457871599795, 2.35221548594714, 2.35994167884065, 2.36776709937289, 2.37569708603435, 2.38373136578869, 2.39188740715789, 2.4001723655651, 2.40858055742088, 2.41712532444938, 2.42582176653978, 2.43467176969305, 2.44367644631041, 2.45286124704223, 2.46222059889919, 2.47177408512471, 2.48152482557902, 2.49149534670758, 2.50167973724668, 2.5121025424396, 2.52276859362371, 2.53368621015369, 2.54471244064013, 2.55579887509383, 2.56691556730856, 2.57807828112741, 2.58929124755218, 2.60052938345122, 2.61181075197994, 2.6231406506567, 2.63450877844716, 2.64592079024619, 2.65736566297975, 2.66886801235028, 2.68039925680379, 2.6919660621504, 2.70359735310195, 2.71524178942724, 2.72695021810889, 2.73869095536778, 2.75047462230107, 2.76231402600025, 2.77417366427405, 2.78609125336844, 2.79805689620922, 2.81005922472766, 2.82211450958847, 2.83421194016232, 2.84633914825058, 2.85851444304052, 2.87072623771478, 2.88299593842877, 2.89531316315953, 2.90766594630464, 2.92007909343311, 2.93254220452969, 2.94500235859685, 2.95756836203094, 2.97014763145451, 2.98276880579325, 2.99546605253435, 3.00818177435952, 3.02095009321979, 3.03375911466721, 3.04664949155654, 3.05955682058631, 3.07252220143771, 3.08547407825541, 3.09857806889149, 3.11164189991262, 3.1247740454825, 3.13803326445637, 3.15127392356254, 3.16454764394741, 3.17791564782097, 3.19129151012199, 3.20473834661461, 3.21824566667958, 3.23171364822652, 3.24322176393691, 3.25693690766748, 3.27057303493195, 3.28439587776679, 3.2982023981508, 3.3120769920673, 3.32600819520374, 3.33998266349085, 3.35398523152793, 3.36799876948903, 3.38213175244099, 3.39637913151413, 3.41059795760726, 3.42490457070842, 3.43914277447917, 3.45358581167583, 3.46808249458969, 3.4824554766011, 3.49700354471659, 3.51172864029507, 3.52627127838015, 3.54114568729358, 3.55579375948015, 3.57056350893963, 3.58544956703106, 3.60044533617299, 3.61554285093282, 3.63050282206819, 3.64552230305604, 3.66083358332543, 3.67594250063318, 3.69132730927766, 3.70645512356833, 3.72183502705042, 3.73747813266146, 3.75278864636778, 3.76832751756716, 3.78410258957844, 3.79978518957619, 3.81533904720756, 3.83144482482593, 3.84703224899097, 3.86316870335808, 3.87910227339303, 3.89519546271496, 3.91101665356257, 3.92739373576442, 3.94345781904698, 3.95961844868298, 3.97635732852111, 3.99270630339808, 4.00911573850308, 4.02556924273202, 4.04204836767728, 4.05853233081845, 4.07560361892911, 4.09205960420983, 4.10909253284836, 4.12608173443567, 4.14299493817811, 4.15979664012907, 4.17719144318364, 4.17719144318364]

            numf = length(fi) - 1
            for i = 1:numf
                f = fi[i]
                @test sa.amplification(f) ≈ Ai[i] rtol = 1e-2
            end

            f_hi = 500.0
            f_max = 999.99
            mms = [SiteAmpUnit(),
                SiteAmpBoore2016_760(),
                SiteAmpAlAtikAbrahamson2021_ask14_620(),
                SiteAmpAlAtikAbrahamson2021_ask14_760(),
                SiteAmpAlAtikAbrahamson2021_ask14_1100(),
                SiteAmpAlAtikAbrahamson2021_bssa14_620(),
                SiteAmpAlAtikAbrahamson2021_bssa14_760(),
                SiteAmpAlAtikAbrahamson2021_bssa14_1100(),
                SiteAmpAlAtikAbrahamson2021_cb14_620(),
                SiteAmpAlAtikAbrahamson2021_cb14_760(),
                SiteAmpAlAtikAbrahamson2021_cb14_1100(),
                SiteAmpAlAtikAbrahamson2021_cy14_620(),
                SiteAmpAlAtikAbrahamson2021_cy14_760(),
                SiteAmpAlAtikAbrahamson2021_cy14_1100()]
            for mm in mms
                @test site_amplification(f_hi, mm) ≈ site_amplification(f_max, mm)
            end

            # @code_warntype site_amplification(f_hi, mms[3])
        end

        @testset "Kappa Filter" begin
            f = 10.0
            # @code_warntype kappa_filter(f, siteAf)
            # @code_warntype kappa_filter(f, siteAd)
            Kff = kappa_filter(f, siteAf)
            Kfd = kappa_filter(f, siteAd)
            @test Kff == Kfd.value

            nf = 100
            fi = exp.(range(log(1e-2), stop=log(1e2), length=nf))
            Kf0fi = kappa_filter(fi, siteAf)
            Kf0di = kappa_filter(fi, siteAd)
            for i in 1:nf
                @test Kf0fi[i] == Kf0di[i].value
            end
            
            Affi = ones(eltype(Kf0fi), nf)
            Afdi = ones(eltype(Kf0di), nf)
            StochasticGroundMotionSimulation.apply_kappa_filter!(Affi, fi, siteAf)
            StochasticGroundMotionSimulation.apply_kappa_filter!(Afdi, fi, siteAd)
            for i in 1:nf
                @test Affi[i] == Afdi[i].value
            end

        end

        @testset "Zeta Filter" begin
            f = 10.0
            # @code_warntype kappa_filter(f, siteAf)
            # @code_warntype kappa_filter(f, siteAd)
            siteAzf = SiteParameters(ζ0f, ηf)
            siteAzd = SiteParameters(ζ0d, ηd)
            Kff = kappa_filter(f, siteAzf)
            Kfd = kappa_filter(f, siteAzd)
            @test Kff == Kfd.value

            nf = 100
            fi = exp.(range(log(1e-2), stop=log(1e2), length=nf))
            Kf0fi = kappa_filter(fi, siteAzf)
            Kf0di = kappa_filter(fi, siteAzd)
            for i in 1:nf
                @test Kf0fi[i] == Kf0di[i].value
            end

            Affi = ones(eltype(Kf0fi), nf)
            Afdi = ones(eltype(Kf0di), nf)
            StochasticGroundMotionSimulation.apply_kappa_filter!(Affi, fi, siteAzf)
            StochasticGroundMotionSimulation.apply_kappa_filter!(Afdi, fi, siteAzd)
            for i in 1:nf
                @test Affi[i] == Afdi[i].value
            end

        end

    end


    @testset "Oscillator" begin
        ζ = 0.05
        f_n = 1.0
        sdof = Oscillator(f_n, ζ)

        @test f_n ≈ 1.0 / period(sdof)

        @test transfer(0.5, sdof)^2 ≈ StochasticGroundMotionSimulation.squared_transfer(0.5, sdof)

        fi = [0.5, 1.0, 2.0]

        # @code_warntype transfer(fi, sdof)

        Hfi = transfer(fi, sdof)
        tfi = 2 * fi
        transfer!(Hfi, tfi, sdof)
        @test Hfi ≈ transfer(tfi, sdof)

        StochasticGroundMotionSimulation.squared_transfer!(Hfi, fi, sdof)
        @test Hfi ≈ transfer(fi, sdof) .^ 2

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
        @test Dsf ≈ 1.0 / faf
        @test isnan(fbf)
        @test isnan(εf)
        fad, fbd, εd = corner_frequency(m, srcd)
        # @code_warntype StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, srcd)
        Dsd = StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, srcd)
        @test Dsd ≈ 1.0 / fad
        @test isnan(fbd)
        @test isnan(εd)

        # @code_warntype StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, fasf)
        @test StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, fasf) ≈ Ds
        @test StochasticGroundMotionSimulation.boore_thompson_2014(m, 0.0, fasd) ≈ Ds
        @test StochasticGroundMotionSimulation.boore_thompson_2015(m, 0.0, fasf, :ACR) ≈ Ds
        @test StochasticGroundMotionSimulation.boore_thompson_2015(m, 0.0, fasd, :ACR) ≈ Ds

        rti = range(0.1, stop=700.0, step=1.0)
        m = -8.0
        rvt_acr = RandomVibrationParameters(:DK80, :BT14, :BT15, :ACR)
        rvt_scr = RandomVibrationParameters(:DK80, :BT15, :BT15, :SCR)
        for i in 1:lastindex(rti)
            if (rti[i] < 17.0) | (rti[i] > 281.0)
                @test excitation_duration(m, rti[i], src, rvt_scr) .<= excitation_duration(m, rti[i], src, rvt_acr)
            else
                @test excitation_duration(m, rti[i], src, rvt_scr) .> excitation_duration(m, rti[i], src, rvt_acr)
            end
        end

        @test excitation_duration(m, 10.0, src, rvt_acr) == excitation_duration(m, 10.0, fas, rvt_acr)

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
        Ds = 0.5 * (1.0 / fa + 1.0 / fb)
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
        fdg = log(Dex2 / Dex1) / h

        d(x) = log(StochasticGroundMotionSimulation.boore_thompson_2014(x[1], 1.0 + near_source_saturation(x[1], fasf), fasf))
        gd(x) = ForwardDiff.gradient(d, x)
        adg = gd([8.0])[1]

        @test fdg ≈ adg atol = 1e-2


        rvt = RandomVibrationParameters()
        # @code_warntype excitation_duration(m, r, fasf, rvt)
        # @code_warntype excitation_duration(m, r, fasd, rvt)
        Dexf = excitation_duration(m, r, fasf, rvt)
        Dexd = excitation_duration(m, r, fasd, rvt)
        @test Dexf == Dexd.value

        rvt = RandomVibrationParameters(:DK80, :BT15, :BT15, :SCR)
        Dexf = excitation_duration(m, r, fasf, rvt)
        Dexd = excitation_duration(m, r, fasd, rvt)
        @test Dexf == Dexd.value


        c11 = [8.4312e-01, -2.8671e-02, 2.0, 1.7316e+00, 1.1695e+00, 2.1671e+00, 9.6224e-01]
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

        @test Dex != Dex0


        # test combinations of :pf_method and :dur_rms methods
        rvt0 = RandomVibrationParameters()
        rvt1 = RandomVibrationParameters(:DK80)
        rvt2 = RandomVibrationParameters(:DK80, :SCR)
        rvt3 = RandomVibrationParameters(:CL56)
        rvt4 = RandomVibrationParameters(:CL56, :SCR)

        @test rvt0 == rvt1
        @test rvt0 != rvt2
        @test rvt0.pf_method == rvt2.pf_method
        @test rvt0.dur_ex == :BT15
        @test rvt0.dur_rms == :BT15
        @test rvt3 != rvt4
        @test rvt3.pf_method != rvt0.pf_method
        @test rvt3.dur_region == :ACR
        @test rvt3.pf_method == rvt4.pf_method
        @test rvt3.dur_rms == :BT12



        # check that the excitation durations are matching
        rvt = RandomVibrationParameters(:CL56, :BT14, :BT12, :ACR)
        Drms, Dex, Dratio = StochasticGroundMotionSimulation.boore_thompson_2012(m, r, fas, sdof, rvt)
        Dex0 = StochasticGroundMotionSimulation.boore_thompson_2014(m, r, fas)

        @test Dex == Dex0

        # confirm that incorrect combinations of peak factor and rms duration gives NaN result
        rvt = RandomVibrationParameters(:DK80, :BT14, :BT12, :ACR)
        Drms, Dex, Dratio = rms_duration(m, r, fas, sdof, rvt)

        @test isnan(Dex)

        # check that rms and excitation duration wrapper functions work as intended
        Dex0 = excitation_duration(m, r, fas, rvt0)
        Drms, Dex1, Dratio = rms_duration(m, r, fas, sdof, rvt0)

        @test Dex0 == Dex1

        Dex0 = excitation_duration(m, r, fas, rvt2)
        Drms, Dex1, Dratio = rms_duration(m, r, fas, sdof, rvt0)

        @test Dex0 != Dex1


        # @code_warntype rms_duration(m, r, srcf, path, sdof, rvt)
        # @code_warntype rms_duration(m, r, srcd, path, sdof, rvt)
        # @code_warntype rms_duration(m, r, fasf, sdof, rvt)
        # @code_warntype rms_duration(m, r, fasd, sdof, rvt)

        # @time Drmsf, Dexf, Dratiof = rms_duration(m, r, fasf, sdof, rvt)
        # @time Drmsd, Dexd, Dratiod = rms_duration(m, r, fasd, sdof, rvt)
        rvt = RandomVibrationParameters()
        Drmsf, Dexf, Dratiof = rms_duration(m, r, fasf, sdof, rvt)
        Drmsd, Dexd, Dratiod = rms_duration(m, r, fasd, sdof, rvt)

        @test Drmsf == Drmsd.value
        @test Dexf == Dexd.value
        @test Dratiof == Dratiod.value


        rvt = RandomVibrationParameters(:CL56, :BT14, :BT12, :ACR)
        Drmsf, Dexf, Dratiof = rms_duration(m, r, fasf, sdof, rvt)
        Drmsf1, Dexf1, Dratiof1 = StochasticGroundMotionSimulation.boore_thompson_2012(m, r, fasf, sdof, rvt)

        @test Drmsf == Drmsf1
        @test Dexf == Dexf1
        @test Dratiof == Dratiof1

        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, srcf, :ACR))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, srcd, :ACR))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, SourceParameters(1), :ACR))

        # Boore & Thompson 2014
        m = 6.0
        fa, fb, ε = corner_frequency(m, src)
        Ds = 1.0 / fa

        Dex270 = Ds + 34.2
        Dex300 = Dex270 + 0.156 * 30.0
        @test Dex270 ≈ StochasticGroundMotionSimulation.boore_thompson_2014(m, 270.0, src)
        @test Dex300 ≈ StochasticGroundMotionSimulation.boore_thompson_2014(m, 300.0, src)

        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2014(m, -1.0, src))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2014(m, -1.0, srcd))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2014(m, Dual(-1.0), src))

        @test isnan(excitation_duration(m, -1.0, src, rvt))
        @test isnan(excitation_duration(m, -1.0, srcd, rvt))
        @test isnan(excitation_duration(m, Dual(-1.0), src, rvt))

        c = StochasticGroundMotionSimulation.StochasticGroundMotionSimulation.boore_thompson_2012_coefs(1, 1, region=:SCR)
        idx = 1
        @test all(isapprox(c, StochasticGroundMotionSimulation.coefs_ena_bt12[idx, 3:9]))

        d1a, d2a, d3a = StochasticGroundMotionSimulation.boore_thompson_2012(6.1234, 2.0, src, sdof, rvt)
        d1b, d2b, d3b = StochasticGroundMotionSimulation.boore_thompson_2012(6.1234, 2.1234, src, sdof, rvt)
        d1c, d2c, d3c = StochasticGroundMotionSimulation.boore_thompson_2012(6.0, 2.1234, src, sdof, rvt)
        d1d, d2d, d3d = StochasticGroundMotionSimulation.boore_thompson_2012(6.0, 2.0, src, sdof, rvt)
        @test d1a < d1b
        @test d2a < d2b
        @test d3a > d3b
        @test d1a > d1c
        @test d1d < d1a


        c = StochasticGroundMotionSimulation.boore_thompson_2015_coefs(1, 1, region=:SCR)
        idx = 1
        @test all(isapprox(c, StochasticGroundMotionSimulation.coefs_ena_bt15[idx, 3:9]))

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

        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2014_path_duration(-1.0))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015_path_duration_acr(-1.0))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015_path_duration_scr(-1.0))

        rpi = range(0.0, stop=1000.0, step=1.0)
        DexAi = zeros(length(rpi))
        DexSi = zeros(length(rpi))
        DpAi = zeros(length(rpi))
        DpSi = zeros(length(rpi))
        rvtA = RandomVibrationParameters(:DK80, :ACR)
        rvtS = RandomVibrationParameters(:DK80, :SCR)
        src = SourceParameters(100.0)
        m = 2.0
        fc, d1, d2 = corner_frequency(m, src)
        Ds = 1.0 / fc

        for (i, r) in enumerate(rpi)
            DexAi[i] = excitation_duration(m, r, src, rvtA)
            DexSi[i] = excitation_duration(m, r, src, rvtS)
            DpAi[i] = StochasticGroundMotionSimulation.boore_thompson_2015_path_duration_acr(r)
            DpSi[i] = StochasticGroundMotionSimulation.boore_thompson_2015_path_duration_scr(r)
        end

        @test all(isapprox.(DexAi .- Ds, DpAi))
        @test all(isapprox.(DexSi .- Ds, DpSi))

        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, rvtA))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, rvtS))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, :PJS))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, :PJS))

        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, fas, rvtA))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, fas, rvtS))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, fas, :PJS))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, fas, :PJS))

        src = SourceParameters(100.0, :Atkinson_Silva_2000)

        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, rvtA))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, rvtS))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, :PJS))
        @test isnan(StochasticGroundMotionSimulation.boore_thompson_2015(m, -1.0, src, :PJS))


        @testset "Edwards (2023) duration" begin
            src = SourceParameters(50.0)
            geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.7], :Piecewise)
            sat = NearSourceSaturationParameters(:BT15)
            ane = AnelasticAttenuationParameters(3000.0, 0.0)
            path = PathParameters(geo, sat, ane)
            site = SiteParameters(0.03, SiteAmpUnit())
            fas = FourierParameters(src, path, site)

            rvt = RandomVibrationParameters(:CL56, :BE23, :LP99, :ACR)

            sdof = Oscillator(1.0)

            m = 6.0
            r_rup = 10.0
            r_ps = equivalent_point_source_distance(r_rup, m, fas)
            Dex = excitation_duration(m, r_ps, fas, rvt)
            (Drms, Dex0, Drat) = rms_duration(m, r_ps, fas, sdof, rvt)
            @test Dex == Dex0

            # @code_warntype rms_duration(m, r_ps, fas, sdof, rvt)
        end

        @testset "UK duration models" begin
            src = SourceParameters(50.0)
            geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.7], :Piecewise)
            sat = NearSourceSaturationParameters(:BT15)
            ane = AnelasticAttenuationParameters(3000.0, 0.0)
            path = PathParameters(geo, sat, ane)
            site = SiteParameters(0.03, SiteAmpUnit())
            fas = FourierParameters(src, path, site)

            rvt_free = RandomVibrationParameters(:DK80, :UKfree, :BT15, :ACR)
            rvt_fixed = RandomVibrationParameters(:DK80, :UKfixed, :BT15, :ACR)

            sdof = Oscillator(1.0)

            m = 1.0
            r_rup = 40.0
            r_ps = equivalent_point_source_distance(r_rup, m, fas)
            Dex_free = excitation_duration(m, r_ps, fas, rvt_free)
            (Drms_free, Dex0_free, Drat_free) = rms_duration(m, r_ps, fas, sdof, rvt_free)
            @test Dex_free == Dex0_free

            Dex_fixed = excitation_duration(m, r_ps, fas, rvt_fixed)
            (Drms_fixed, Dex0_fixed, Drat_fixed) = rms_duration(m, r_ps, fas, sdof, rvt_fixed)
            @test Dex_fixed == Dex0_fixed

            # @code_warntype StochasticGroundMotionSimulation.uk_path_duration_free(r_ps)

            m = 0.0
            for r_rup in range(-4.0, stop=600.0, step=5.0)
                r_ps = equivalent_point_source_distance(r_rup, m, fas)
                # StochasticGroundMotionSimulation.uk_path_duration_free(r_ps)
                # StochasticGroundMotionSimulation.uk_duration_free(m, r_ps, src)
                Dex_free = excitation_duration(m, r_ps, fas, rvt_free)
                (Drms_free, Dex0_free, Drat_free) = rms_duration(m, r_ps, fas, sdof, rvt_free)
                @test Dex_free == Dex0_free

                Dex_fixed = excitation_duration(m, r_ps, fas, rvt_fixed)
                (Drms_fixed, Dex0_fixed, Drat_fixed) = rms_duration(m, r_ps, fas, sdof, rvt_fixed)
                @test Dex_fixed == Dex0_fixed
            end
        end

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
        geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0, 1]), :Piecewise)

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

        @testset "Fourier Source Shape" begin
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

            @test Afff ≈ 1.0 atol = 1e-3

            fa, fb, ε = corner_frequency(m, srcf)
            # @code_warntype fourier_source_shape(f, fa, fb, ε, srcf.model)
            Afc = fourier_source_shape(f, fa, fb, ε, srcf)
            @test Afc ≈ 1.0 atol = 1e-3

            fa, fb, ε = corner_frequency(m, srcd)
            # @code_warntype fourier_source_shape(f, fa, fb, ε, srcd.model)
            Afcd = fourier_source_shape(f, fa, fb, ε, srcd)
            @test Afcd ≈ 1.0 atol = 1e-3

            src_a = SourceParameters(100.0, :Atkinson_Silva_2000)
            Af_a = fourier_source_shape(f, m, src_a)
            @test Af_a ≈ 1.0 atol = 1e-3
            src_n = SourceParameters(100.0, :Null)
            Af_n = fourier_source_shape(f, m, src_n)
            src_b = SourceParameters(100.0)
            Af_b = fourier_source_shape(f, m, src_b)
            @test Af_n == Af_b

            fa, fb, ε = corner_frequency(m, src_a)
            Af_a = fourier_source_shape(f, fa, fb, ε, src_a)
            @test Af_a ≈ 1.0 atol = 1e-3
            fa, fb, ε = corner_frequency(m, src_b)
            Af_n = fourier_source_shape(f, fa, fb, ε, src_n)
            @test Af_n ≈ Af_b


            f = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
            nf = length(f)

            # @code_warntype fourier_source_shape(f, m, srcf)
            # @code_warntype fourier_source_shape(f, m, srcd)
            Affs = fourier_source_shape(f, m, srcf)
            Afff = fourier_source_shape(f, m, fasf)
            Afds = fourier_source_shape(f, m, srcd)
            Afdf = fourier_source_shape(f, m, fasd)
            for i in 1:nf
                @test Affs[i] == Afds[i].value
                @test Afff[i] == Afdf[i].value
                @test Affs[i] == Afff[i]
                @test Afds[i] == Afdf[i]
            end

            fa, fb, ε = corner_frequency(m, srcf)
            # @code_warntype fourier_source_shape(f, fa, fb, ε, srcf.model)
            Afc = fourier_source_shape(f, fa, fb, ε, srcf)
            @test Afc[1] ≈ 1.0 atol = 1e-3

            fa, fb, ε = corner_frequency(m, srcd)
            # @code_warntype fourier_source_shape(f, fa, fb, ε, srcd.model)
            Afcd = fourier_source_shape(f, fa, fb, ε, srcd)
            @test Afcd[1] ≈ 1.0 atol = 1e-3

            src_a = SourceParameters(100.0, :Atkinson_Silva_2000)
            Af_a = fourier_source_shape(f, m, src_a)
            @test Af_a[1] ≈ 1.0 atol = 1e-3
            src_n = SourceParameters(100.0, :Null)
            Af_n = fourier_source_shape(f, m, src_n)
            src_b = SourceParameters(100.0)
            Af_b = fourier_source_shape(f, m, src_b)
            @test Af_n == Af_b

            fa, fb, ε = corner_frequency(m, src_a)
            Af_a = fourier_source_shape(f, fa, fb, ε, src_a)
            @test Af_a[1] ≈ 1.0 atol = 1e-3
            fa, fb, ε = corner_frequency(m, src_b)
            Af_n = fourier_source_shape(f, fa, fb, ε, src_n)
            @test Af_n ≈ Af_b

            src = SourceParameters(100.0, :Beresnev_2019)
            fa, fb, ε = corner_frequency(m, src)
            Af = fourier_source_shape(f, fa, fb, ε, src)
            @test Af[1] ≈ 1.0 atol = 1e-3

            fasbf = FourierParameters(src, pathf, sitef)
            fa, fb, ε = corner_frequency(m, fasbf)
            Af = fourier_source_shape(f, fa, fb, ε, fasbf)
            @test Af[1] ≈ 1.0 atol = 1e-3

        end

        @testset "Fourier Source" begin
            f = 0.001
            m = 6.0
            # @code_warntype fourier_source(f, m, srcf)
            # @code_warntype fourier_source(f, m, srcd)
            Afs = fourier_source(f, m, srcf)
            Aff = fourier_source(f, m, fasf)
            @test Afs == Aff
            # @time fourier_source(f, m, srcd)
            # @time fourier_source(f, m, fasd)
        end

        @testset "Fourier Path" begin
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
        end


        @testset "Fourier attenuation" begin
            f = 10.0
            r = 200.0
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

            f = [0.01, 0.1, 1.0, 10.0, 100.0]
            nf = length(f)
            Qf = fourier_attenuation(f, r, fasf)
            Qd = fourier_attenuation(f, r, fasd)
            Qm = fourier_attenuation(f, r, fasm)
            @test Qf == map(q -> q.value, Qd)
            @test Qd == Qm

            Aff = ones(eltype(Qf), nf)
            Afd = ones(eltype(Qd), nf)
            Afm = ones(eltype(Qm), nf)
            StochasticGroundMotionSimulation.apply_fourier_attenuation!(Aff, f, r, fasf)
            StochasticGroundMotionSimulation.apply_fourier_attenuation!(Afd, f, r, fasd)
            StochasticGroundMotionSimulation.apply_fourier_attenuation!(Afm, f, r, fasm)
            @test Aff == map(a -> a.value, Afd)
            @test Aff == map(a -> a.value, Afm)

        end


        @testset "Fourier site" begin
            # @code_warntype fourier_site(f, sitef)
            # @code_warntype fourier_site(f, sited)
            # @time fourier_site(f, sitef)
            # @time fourier_site(f, sited)
            f = 0.5

            Sf = fourier_site(f, fasf)
            Sd = fourier_site(f, fasd)
            Sm = fourier_site(f, fasm)
            @test Sf == Sd.value
            @test Sd == Sm
        end

        @testset "Point source distance" begin
            f = 1.0
            m = 6.0
            r = 10.0
            r_psf = equivalent_point_source_distance(r, m, fasf)
            r_psd = equivalent_point_source_distance(r, m, fasd)
            r_psm = equivalent_point_source_distance(r, m, fasm)
            @test r_psf ≈ r_psd
            @test r_psf ≈ r_psm
        end

        # @code_warntype fourier_spectral_ordinate(f, m, r_psf, fasf)
        # @code_warntype fourier_spectral_ordinate(f, m, r_psd, fasd)
        # @code_warntype fourier_spectral_ordinate(f, m, r_psm, fasm)

        f = 1.0
        m = 6.0
        r = 10.0
        r_psf = equivalent_point_source_distance(r, m, fasf)
        r_psd = equivalent_point_source_distance(r, m, fasd)
        r_psm = equivalent_point_source_distance(r, m, fasm)

        Af = fourier_spectral_ordinate(f, m, r_psf, fasf)
        Ad = fourier_spectral_ordinate(f, m, r_psd, fasd)
        Am = fourier_spectral_ordinate(f, m, r_psm, fasm)
        @test Af == Ad.value
        @test Ad == Am


        # test Beresnev source spectrum
        srcb1p0 = SourceParameters(100.0, 1.0)
        srcb1p5 = SourceParameters(100.0, 1.5)

        fasb1p0 = FourierParameters(srcb1p0, pathf, sitef)
        fasb1p5 = FourierParameters(srcb1p5, pathf, sitef)

        f = 10.0
        m = 6.0
        r = 10.0
        r_ps = equivalent_point_source_distance(r, m, fasb1p0)

        Afb1p0 = fourier_spectral_ordinate(f, m, r_ps, fasb1p0)
        Afb1p5 = fourier_spectral_ordinate(f, m, r_ps, fasb1p5)

        @test Afb1p0 > Afb1p5

        f = 1e-3
        Afb1p0 = fourier_spectral_ordinate(f, m, r_ps, fasb1p0)
        Afb1p5 = fourier_spectral_ordinate(f, m, r_ps, fasb1p5)

        @test Afb1p0 ≈ Afb1p5 rtol = 1e-3


        fi = [0.01, 0.1, 1.0, 10.0, 100.0]

        # @code_warntype fourier_spectrum(fi, m, r_psf, fasf)
        # @code_warntype fourier_spectrum(fi, m, r_psf, fasd)
        # @code_warntype fourier_spectrum(fi, m, r_psd, fasd)
        # @code_warntype fourier_spectrum(fi, m, r_psd, fasm)

        Afif = fourier_spectrum(fi, m, r_psf, fasf)
        Afid = fourier_spectrum(fi, m, r_psf, fasd)
        Afim = fourier_spectrum(fi, m, r_psd, fasm)
        for i = 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))

        fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        fourier_spectrum!(Afim, fi, m, r_psd, fasm)
        for i = 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))

        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)
        for i = 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))


        

        anef = AnelasticAttenuationParameters(Q0f, ηf, :Rrup)
        aned = AnelasticAttenuationParameters(Q0d, ηd, :Rrup)
        anem = AnelasticAttenuationParameters(Q0f, ηd, :Rrup)

        pathf = PathParameters(geof, sat, anef)
        pathd = PathParameters(geod, sat, aned)
        pathm = PathParameters(geom, sat, anem)

        fasf = FourierParameters(srcf, pathf, sitef)
        fasd = FourierParameters(srcd, pathd, sited)
        fasm = FourierParameters(srcf, pathm, sited)

        
        r_psf = equivalent_point_source_distance(r, m, fasf)
        r_psd = equivalent_point_source_distance(r, m, fasd)
        r_psm = equivalent_point_source_distance(r, m, fasm)

        Afif = fourier_spectrum(fi, m, r_psf, fasf)
        Afid = fourier_spectrum(fi, m, r_psf, fasd)
        Afim = fourier_spectrum(fi, m, r_psd, fasm)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afim, fi, m, r_psd, fasm)
        for i = 1:length(fi)
            @test Afif[i] == Afid[i].value
        end
        @test all(isapprox.(Afid, Afim))

        # test parallel threading of fourier spectrum computation
        # fi = exp10.(range(-2.0, stop=2.0, length=31))
        # Afif = fourier_spectrum(fi, m, r_psf, fasf)

        # @benchmark StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        # @benchmark StochasticGroundMotionSimulation.squared_fourier_spectrum_par!(Afif, fi, m, r_psf, fasf)


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
        @test any(isapprox.(sqAfid, Afid .^ 2))

        # @code_warntype fourier_spectrum!(Afif, fi, m, r_psf, fasf)
        # @code_warntype fourier_spectrum!(Afid, fi, m, r_psf, fasd)
        # @code_warntype fourier_spectrum!(Afim, fi, m, r_psd, fasm)


        Afid = fourier_spectrum(Vector{Float64}(), m, r_ps, fas)
        fourier_spectrum!(Afid, Vector{Float64}(), m, r_ps, fas)
        StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, Vector{Float64}(), m, r_ps, fas)

        Afid = fourier_spectrum(fi, m, r_ps, fas)
        fourier_spectrum!(Afid, fi, m, r_ps, fas)
        StochasticGroundMotionSimulation.get_parametric_type(fas)
        StochasticGroundMotionSimulation.get_parametric_type(fas.path.anelastic)

        fi = exp.(range(log(1e-2), stop=log(1e2), length=1000))
        Afid = fourier_spectrum(fi, m, r_ps, fas)
        # @benchmark StochasticGroundMotionSimulation.fourier_spectrum!(Afid, fi, m, r_ps, fas)
        # @benchmark StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_ps, fas)
        # @ballocated StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_ps, fas)
        # @profview @benchmark StochasticGroundMotionSimulation.fourier_spectrum!(Afid, fi, m, r_ps, fas)
        # @profview @benchmark StochasticGroundMotionSimulation.squared_fourier_spectrum!(Afid, fi, m, r_ps, fas)
        # @benchmark Afsqid = StochasticGroundMotionSimulation.squared_fourier_spectrum(fi, m, r_ps, fas)

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
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0, 1]), :Piecewise)

            anef = AnelasticAttenuationParameters(Q0f, ηf)
            aned = AnelasticAttenuationParameters(Q0d, ηd)
            anem = AnelasticAttenuationParameters(Q0f, ηd)

            sat = NearSourceSaturationParameters(:BT15)

            pathf = PathParameters(geof, sat, anef)
            pathd = PathParameters(geod, sat, aned)
            pathm = PathParameters(geom, sat, anem)

            sitef = SiteParameters(κ0f, SiteAmpBoore2016_760())
            sited = SiteParameters(κ0d, SiteAmpBoore2016_760())

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
            ps2db(f) = (1.0 / (2π * sdof.f_n))^2 * 1e4

            dbm0_integrand(f) = StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm0ln_integrand(lnf) = exp(lnf) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm1_integrand(f) = (2π * f) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm1ln_integrand(lnf) = exp(lnf) * (2π * exp(lnf)) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm2_integrand(f) = (2π * f)^2 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm2ln_integrand(lnf) = exp(lnf) * (2π * exp(lnf))^2 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            dbm4_integrand(f) = (2π * f)^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * ps2db(f)
            dbm4ln_integrand(lnf) = exp(lnf) * (2π * exp(lnf))^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof) * ps2db(exp(lnf))

            @time igk = 2 * quadgk(dbm0_integrand, 0.0, Inf)[1]
            @time igk = 2 * quadgk(dbm0_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2 * StochasticGroundMotionSimulation.gauss_intervals(dbm0ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol = 1e-4

            @time igk = 2 * quadgk(dbm1_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2 * StochasticGroundMotionSimulation.gauss_intervals(dbm1ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol = 1e-4

            @time igk = 2 * quadgk(dbm2_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2 * StochasticGroundMotionSimulation.gauss_intervals(dbm2ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol = 1e-3

            @time igk = 2 * quadgk(dbm4_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = 2 * StochasticGroundMotionSimulation.gauss_intervals(dbm4ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)
            @test igk ≈ iglelnm rtol = 1e-3


            m0_integrand(f) = StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m0ln_integrand(lnf) = exp(lnf) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)

            m1_integrand(f) = (2π * f) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m1ln_integrand(lnf) = exp(lnf) * (2π * exp(lnf)) * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)

            m2_integrand(f) = (2π * f)^2 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m2ln_integrand(lnf) = exp(lnf) * (2π * exp(lnf))^2 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)

            m4_integrand(f) = (2π * f)^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(f, m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(f, sdof)
            m4ln_integrand(lnf) = exp(lnf) * (2π * exp(lnf))^4 * StochasticGroundMotionSimulation.squared_fourier_spectral_ordinate(exp(lnf), m, r_psf, fasf) * StochasticGroundMotionSimulation.squared_transfer(exp(lnf), sdof)


            @time igk = quadgk(m0_integrand, exp(-7.0), exp(7.0))[1]
            @time igle = StochasticGroundMotionSimulation.gauss_interval(m0_integrand, 2000, 0.0, 300.0)
            @time igleln = StochasticGroundMotionSimulation.gauss_interval(m0ln_integrand, 750, -7.0, 7.0)

            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m0ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            # @test igk ≈ igle rtol=1e-2
            @test igk ≈ igleln rtol = 1e-4
            @test igk ≈ iglelnm rtol = 1e-4

            lnfi = log.([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, sdof.f_n])
            sort!(lnfi)

            @time igk = quadgk(m1_integrand, exp(-7.0), exp(7.0))[1]
            @time igle = StochasticGroundMotionSimulation.gauss_interval(m1_integrand, 1500, 0.0, 300.0)
            @time igleln = StochasticGroundMotionSimulation.gauss_interval(m1ln_integrand, 750, -7.0, 7.0)
            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m1ln_integrand, 250, -7.0, log(sdof.f_n), 7.0)

            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m1ln_integrand, 30, lnfi...)

            @time itr = StochasticGroundMotionSimulation.trapezoidal(m1ln_integrand, 60, lnfi...)

            # @test igk ≈ igle rtol=1e-2
            @test igk ≈ iglelnm rtol = 1e-4
            @test igk ≈ itr rtol = 1e-3


            @time igk = quadgk(m2_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m2ln_integrand, 30, lnfi...)
            @time itr = StochasticGroundMotionSimulation.trapezoidal(m2ln_integrand, 60, lnfi...)

            @test igk ≈ iglelnm rtol = 1e-4
            @test igk ≈ itr rtol = 1e-3

            @time igk = quadgk(m4_integrand, exp(-7.0), exp(7.0))[1]
            @time iglelnm = StochasticGroundMotionSimulation.gauss_intervals(m4ln_integrand, 30, lnfi...)
            @time itr = StochasticGroundMotionSimulation.trapezoidal(m4ln_integrand, 60, lnfi...)

            @test igk ≈ iglelnm rtol = 1e-4
            @test igk ≈ itr rtol = 1e-3


            integrand(x) = sin(x)

            intervals = 101
            x_min = 0.0
            x_max = 2π
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            ist = StochasticGroundMotionSimulation.trapezoidal_rule(xx, yy)

            igk = quadgk(integrand, x_min, x_max)[1]

            @test isr ≈ igk atol = 10 * eps()
            @test ist ≈ igk atol = 10 * eps()


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

            integrand1(f) = StochasticGroundMotionSimulation.squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = sdof.f_n / 1.1
            x_max = sdof.f_n * 1.1
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand1.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            igk = quadgk(integrand1, x_min, x_max)[1]

            @test isr ≈ igk rtol = 1e-6
            @test isr ≈ igk atol = 1e-6


            integrand2(f) = StochasticGroundMotionSimulation.squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = 100.0
            x_max = 200.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand2.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            igk = quadgk(integrand2, x_min, x_max)[1]

            @test isr ≈ igk rtol = 1e-3
            @test isr ≈ igk atol = 1e-6


            integrand3(f) = (2π * f)^4 * StochasticGroundMotionSimulation.squared_transfer(f, sdof) * fourier_spectral_ordinate(f, m, r, fas)^2

            intervals = 101
            x_min = 100.0
            x_max = 200.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand3.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            igk = quadgk(integrand3, x_min, x_max)[1]

            @test isr ≈ igk rtol = 1e-3
            @test isr ≈ igk atol = 1e-6

            x_min = 300.0
            x_max = 500.0
            xx = collect(range(x_min, stop=x_max, length=intervals))
            yy = integrand3.(xx)

            isr = StochasticGroundMotionSimulation.simpsons_rule(xx, yy)
            igk = quadgk(integrand3, x_min, x_max)[1]

            @test isr ≈ igk rtol = 1e-3
            @test isr ≈ igk atol = 1e-6

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
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0, 1]), :Piecewise)

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
            @test m0f.m0 ≈ m0d.m0.value
            @test m0d.m0 ≈ m0m.m0

            order = 0
            m0f = spectral_moment(order, m, r_psf, fasf, sdof, nodes=50)
            m0d = spectral_moment(order, m, r_psd, fasd, sdof, nodes=50)
            m0m = spectral_moment(order, m, r_psm, fasm, sdof, nodes=50)
            @test m0f.m0 ≈ m0d.m0.value
            @test m0d.m0 ≈ m0m.m0

            # @code_warntype spectral_moment(order, m, r_psf, fasf, sdof)
            # @code_warntype spectral_moment(order, m, r_psd, fasd, sdof)
            # @code_warntype spectral_moment(order, m, r_psm, fasm, sdof)

            # @time spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            # @time spectral_moments([0, 1, 2, 4], m, r_psd, fasd, sdof)
            # @time spectral_moments([0, 1, 2, 4], m, r_psm, fasm, sdof)

            # @code_warntype spectral_moments([0, 1, 2, 4], m, r_psf, fasf, sdof)
            # @code_warntype spectral_moments([0, 1, 2, 4], m, r_psd, fasd, sdof)
            # @code_warntype spectral_moments([0, 1, 2, 4], m, r_psm, fasm, sdof)

            smi = spectral_moments([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)
            smigk = StochasticGroundMotionSimulation.spectral_moments_gk([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)

            @test smi.m0 ≈ smigk.m0 rtol=1e-3
            @test smi.m1 ≈ smigk.m1 rtol=1e-3
            @test smi.m2 ≈ smigk.m2 rtol=1e-3
            @test smi.m3 ≈ smigk.m3 rtol=1e-3
            @test smi.m4 ≈ smigk.m4 rtol=1e-3

            smi = spectral_moments([0, 1, 2, 3, 4], m, r_psf, fasf, sdof, nodes=50)
            smigk = StochasticGroundMotionSimulation.spectral_moments_gk([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)

            @test smi.m0 ≈ smigk.m0 rtol = 1e-5
            @test smi.m1 ≈ smigk.m1 rtol = 1e-5
            @test smi.m2 ≈ smigk.m2 rtol = 1e-5
            @test smi.m3 ≈ smigk.m3 rtol = 1e-5
            @test smi.m4 ≈ smigk.m4 rtol = 1e-5


            sdof = Oscillator(1 / 3)
            smi = spectral_moments([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)
            # smai = spectral_moments_alt([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)
            smigk = StochasticGroundMotionSimulation.spectral_moments_gk([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)

            # @code_warntype spectral_moments([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)
            # @code_warntype spectral_moments_alt([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)

            # @benchmark spectral_moments([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)
            # @benchmark spectral_moments_alt([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)

            # @benchmark spectral_moments([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)
            # @profview @benchmark spectral_moments_alt([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)

            @test smi.m0 ≈ smigk.m0 rtol = 1e-3
            @test smi.m1 ≈ smigk.m1 rtol = 1e-3
            @test smi.m2 ≈ smigk.m2 rtol = 1e-3
            @test smi.m3 ≈ smigk.m3 rtol = 1e-3
            @test smi.m4 ≈ smigk.m4 rtol = 1e-3

            m = Dual(6.0)
            smdi = spectral_moments([0, 1, 2, 3, 4], m, r_psf, fasf, sdof)
            smfi = spectral_moments([0, 1, 2, 3, 4], m.value, r_psf, fasf, sdof)

            @test smdi.m0 ≈ smfi.m0 rtol = 1e-3
            @test smdi.m1 ≈ smfi.m1 rtol = 1e-3
            @test smdi.m2 ≈ smfi.m2 rtol = 1e-3
            @test smdi.m3 ≈ smfi.m3 rtol = 1e-3
            @test smdi.m4 ≈ smfi.m4 rtol = 1e-3

            m = Dual(6.0)
            smdi = spectral_moment(0, m, r_psf, fasf, sdof)
            smfi = spectral_moment(0, m.value, r_psf, fasf, sdof)
            @test smdi.m0 ≈ smfi.m0

            rvt = RandomVibrationParameters()
            Ti = [0.01, 0.1, 1.0]
            m = Dual(6.0)
            Sadi = rvt_response_spectrum(Ti, m, 10.0, fasf, rvt)
            Safi = rvt_response_spectrum(Ti, m.value, 10.0, fasf, rvt)
            for i = 1:length(Ti)
                @test Sadi[i].value ≈ Safi[i]
            end

            rvt_response_spectrum!(Sadi, Ti, m, 10.0, fasf, rvt)
            rvt_response_spectrum!(Safi, Ti, m.value, 10.0, fasf, rvt)
            for i = 1:length(Ti)
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
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0, 1]), :Piecewise)

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

            # @code_warntype StochasticGroundMotionSimulation.peak_factor_dk80(m, r_psf, fasf, sdof)

            @test pfps ≈ pfgk rtol = 1e-6
            @test pfpsn ≈ pfgk rtol = 1e-6

            @time pfps = StochasticGroundMotionSimulation.peak_factor_cl56(m, r_psf, fasf, sdof)
            @time pfpsn = StochasticGroundMotionSimulation.peak_factor_cl56(m, r_psf, fasf, sdof, nodes=40)
            @time pfgk = StochasticGroundMotionSimulation.peak_factor_cl56_gk(m, r_psf, fasf, sdof)

            @test pfps ≈ pfgk rtol = 1e-5
            @test pfpsn ≈ pfgk rtol = 1e-5

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
            Dexf = excitation_duration(m, r_ps, fasf, rvt)
            Dexd = excitation_duration(Dual(m), r_ps, fasf, rvt)
            m0f = spectral_moment(0, m, r_ps, fasf, sdof)
            m0d = spectral_moment(0, Dual(m), r_ps, fasf, sdof)
            rvt = RandomVibrationParameters(:PS)
            pf = peak_factor(Dexf, m0f, rvt)
            @test isnan(pf)
            pf = peak_factor(Dexd, m0d, rvt)
            @test isnan(pf)
            pf = peak_factor(Dexf, m0d, rvt)
            @test isnan(pf)
            rvt = RandomVibrationParameters(:CL56)
            pf = peak_factor(Dexf, m0f, rvt)
            @test isnan(pf)
            # @test pf ≈ peak_factor(m, r_ps, fasf, sdof, RandomVibrationParameters(:CL56))

            m0f = spectral_moments([0, 1, 2, 3, 4], m, r_ps, fasf, sdof)
            m0d = spectral_moments([0, 1, 2, 3, 4], Dual(m), r_ps, fasf, sdof)
            pff = peak_factor(Dexf, m0f, rvt)
            pfd = peak_factor(Dexd, m0d, rvt)
            pfm = peak_factor(Dexf, m0d, rvt)
            @test pff ≈ pfd.value
            @test pff ≈ pfm.value

            pfi = StochasticGroundMotionSimulation.peak_factor_integrand_cl56(0.0, 10.0, 10.0)
            @test pfi ≈ 1.0
            pfi = StochasticGroundMotionSimulation.peak_factor_integrand_cl56(Inf, 10.0, 10.0)
            @test pfi ≈ 0.0

            pfi = StochasticGroundMotionSimulation.peak_factor_integrand_cl56(0.0, 10.0, 10.0)
            @test pfi ≈ 1.0
            pfi = StochasticGroundMotionSimulation.peak_factor_integrand_cl56(Inf, 10.0, 10.0)
            @test pfi ≈ 0.0

            pf0 = StochasticGroundMotionSimulation.peak_factor_cl56(10.0, 10.0)
            pf1 = StochasticGroundMotionSimulation.peak_factor_cl56(10.0, 10.0, nodes=50)
            @test pf0 ≈ pf1 rtol=1e-7

            pf2 = StochasticGroundMotionSimulation.peak_factor_cl56(10.0, m0f)
            pf3 = StochasticGroundMotionSimulation.peak_factor_cl56(10.0, m0f, nodes=50)
            @test pf2 ≈ pf3 rtol=1e-6

            pf4 = StochasticGroundMotionSimulation.peak_factor_dk80(10.0, m0f)
            pf5 = StochasticGroundMotionSimulation.peak_factor_dk80(10.0, m0f, nodes=50)
            @test pf4 ≈ pf5 rtol = 1e-7

            rvt = RandomVibrationParameters(:DK80)
            @test rvt.dur_rms == :BT15

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
            geom = GeometricSpreadingParameters(Rrefi, [γ1f], [γ2d], BitVector([0, 1]), :Piecewise)

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

            Ti = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0]

            rvt = RandomVibrationParameters()

            Saif = rvt_response_spectrum(Ti, m, r_psf, fasf, rvt)
            Said = rvt_response_spectrum(Ti, m, r_psd, fasd, rvt)
            Saim = rvt_response_spectrum(Ti, m, r_psm, fasm, rvt)
            for i = 1:length(Ti)
                @test Saif[i] ≈ Said[i].value
            end
            @test all(isapprox.(Said, Saim))


            # test Beresnev source spectrum
            srcb1p0 = SourceParameters(100.0, 1.0)
            srcb1p5 = SourceParameters(100.0, 1.5)

            fasb1p0 = FourierParameters(srcb1p0, pathf, sitef)
            fasb1p5 = FourierParameters(srcb1p5, pathf, sitef)

            T = 0.05
            m = 6.0
            r = 10.0
            r_ps = equivalent_point_source_distance(r, m, fasb1p0)

            Sab1p0 = rvt_response_spectral_ordinate(T, m, r_ps, fasb1p0, rvt)
            Sab1p5 = rvt_response_spectral_ordinate(T, m, r_ps, fasb1p5, rvt)

            @test Sab1p0 > Sab1p5
            @test Sab1p5 < Sab1p0

            rvt = RandomVibrationParameters(:CL56)
            Sab1p0 = rvt_response_spectral_ordinate(T, m, r_ps, fasb1p0, rvt)
            Sab1p5 = rvt_response_spectral_ordinate(T, m, r_ps, fasb1p5, rvt)

            @test Sab1p0 > Sab1p5
            @test Sab1p5 < Sab1p0


            # @code_warntype rvt_response_spectrum(Ti, m, r_psf, fasf, rvt)
            # @code_warntype rvt_response_spectrum(Ti, m, r_psd, fasd, rvt)
            # @code_warntype rvt_response_spectrum(Ti, m, r_psm, fasm, rvt)

        end
    end
end
