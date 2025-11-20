using Test
using StochasticGroundMotionSimulation
using ForwardDiff

@testset "Custom Models" begin
    @testset "Functional FAS Model" begin
        # Define a simple exponential decay FAS model
        fas_func = (f, M, R, params) -> params[1] * exp(-params[2] * f * R)
        params = [1.0, 0.01]
        fas_model = FunctionalFASModel(fas_func, params)

        # Call the model directly as a function
        @test fas_model(1.0, 6.0, 10.0) ≈ 1.0 * exp(-0.01 * 1.0 * 10.0)

        # Test compute_fas function
        @test compute_fas(fas_model, 1.0, 6.0, 10.0) ≈ 1.0 * exp(-0.01 * 1.0 * 10.0)

        # Test ForwardDiff compatibility
        grad = ForwardDiff.gradient(p -> compute_fas(FunctionalFASModel(fas_func, p),
                1.0, 6.0, 10.0), [1.0, 0.01])
        @test length(grad) == 2
        @test !any(isnan.(grad))
    end

    @testset "Functional Duration Model" begin
        # Simple duration function
        dur_func = (M, R, params) -> params[1] + params[2] * M
        dur_params = [1.0, 0.5]
        dur_model = FunctionalDurationModel(dur_func, dur_params)

        # Test basic computation
        @test compute_duration(dur_model, 6.0, 10.0) ≈ 1.0 + 0.5 * 6.0

        # Test ForwardDiff compatibility
        grad = ForwardDiff.gradient(p -> compute_duration(FunctionalDurationModel(dur_func, p),
                6.0, 10.0), [1.0, 0.5])
        @test length(grad) == 2
        @test !any(isnan.(grad))
    end

    @testset "Custom Duration Model with log10" begin
        # Duration function with log10
        dur_func = (m, r, p) -> p[1] + p[2] * m + p[3] * log10(r)
        model = FunctionalDurationModel(dur_func, [1.0, 0.5, 0.3])

        @test compute_duration(model, 6.0, 10.0) ≈ 1.0 + 0.5 * 6.0 + 0.3 * 1.0
    end

    @testset "Integrated RVT Calculation with Custom Models" begin
        # Define a more realistic FAS model (simplified)
        fas_func = (f, M, R, params) -> begin
            # params: [constant, corner_freq_coef, Q0, kappa]
            C, fc_coef, Q0, kappa = params

            # Source (Brune spectrum)
            fc = fc_coef * 10^((M - 6.0) / 3.0)  # simplified corner frequency
            source = 1.0 / (1.0 + (f / fc)^2)

            # Path attenuation
            Q = Q0 * f^0.5
            path = exp(-π * f * R / (Q * 3.5))

            # Site attenuation
            site = exp(-π * kappa * f)

            return C * source * path * site
        end

        fas_model = FunctionalFASModel(fas_func, [1e-3, 5.0, 200.0, 0.04])

        # Define duration model
        dur_func = (M, R, params) -> params[1] + params[2] * M + params[3] * log10(R)
        dur_model = FunctionalDurationModel(dur_func, [1.0, 0.5, 0.3])

        # Compute response spectral ordinate using custom RVT function
        T = 1.0  # period in seconds
        mag = 6.0
        dist = 10.0

        Sa = rvt_response_spectral_ordinate_custom(T, mag, dist, fas_model, dur_model)

        @test Sa > 0
        @test !isnan(Sa)
        @test !isinf(Sa)

        println("  Sa(T=$T s) = $Sa g")
    end

    @testset "Backward Compatibility with Traditional Interface" begin
        # Test that old-style calls still work
        src = SourceParameters(100.0)
        geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
        ane = AnelasticAttenuationParameters(200.0, 0.4)
        sat = NearSourceSaturationParameters(:BT15)
        path = PathParameters(geo, sat, ane)
        site = SiteParameters(0.039)

        fas = FourierParameters(src, path, site)
        rvt = RandomVibrationParameters(:DK80)

        Sa_old = rvt_response_spectral_ordinate(1.0, 6.0, 10.0, fas, rvt)
        @test Sa_old > 0
        @test !isnan(Sa_old)

        println("  Traditional interface: Sa = $Sa_old g")
    end

    @testset "FourierParametersWrapper with Traditional Parameters" begin
        # Create traditional parameter objects
        src = SourceParameters(100.0)
        geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
        ane = AnelasticAttenuationParameters(200.0, 0.4)
        sat = NearSourceSaturationParameters(:BT15)
        path = PathParameters(geo, sat, ane)
        site = SiteParameters(0.039)

        # Create FourierParameters
        fourier_params = FourierParameters(src, path, site)

        # Wrap it to use with the AbstractFASModel interface
        wrapper = FourierParametersWrapper(fourier_params)

        # Test compute_fas
        fas_value = compute_fas(wrapper, 1.0, 6.0, 10.0)
        @test fas_value > 0
        @test !isnan(fas_value)

        println("  FAS via wrapper = $fas_value")
    end

    @testset "Validation Functions" begin
        # Test FAS model validation
        fas_func = (f, M, R, params) -> params[1] * exp(-params[2] * f)
        fas_model = FunctionalFASModel(fas_func, [1.0, 0.01])

        @test validate_fas_model(fas_model, mag=6.0, dist=10.0)

        # Test duration model validation
        dur_func = (M, R, params) -> params[1] + params[2] * M
        dur_model = FunctionalDurationModel(dur_func, [1.0, 0.5])

        @test validate_duration_model(dur_model, mag=6.0, dist=10.0)
    end
end
