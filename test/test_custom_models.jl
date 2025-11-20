using Test
using StochasticGroundMotionSimulation
using StochasticGroundMotionSimulation.CustomModels
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


    @testset "Custom Models - Comprehensive Coverage" begin

        # Test parameters
        test_freq = 1.0
        test_mag = 6.5
        test_dist = 20.0
        test_period = 1.0
        test_damping = 0.05

        @testset "Abstract Types" begin
            @test AbstractFASModel isa Type
            @test AbstractDurationModel isa Type
            @test isabstracttype(AbstractFASModel)
            @test isabstracttype(AbstractDurationModel)
        end

        @testset "FunctionalFASModel" begin
            # Basic functionality
            simple_fas = (f, m, r, p) -> p[1] * f^(-p[2]) * exp(-p[3] * f)
            params = [100.0, 1.0, 0.04]
            model = FunctionalFASModel(simple_fas, params)

            @test model isa AbstractFASModel
            @test model.func === simple_fas
            @test model.params == params

            # compute_fas
            result = compute_fas(model, test_freq, test_mag, test_dist)
            @test result isa Real
            @test result > 0
            @test isfinite(result)

            # Test with different input types
            result_f32 = compute_fas(model, Float32(test_freq), Float32(test_mag), Float32(test_dist))
            @test result_f32 isa Float32

            # Test vectorized frequencies
            freqs = [0.1, 1.0, 10.0]
            results = [compute_fas(model, f, test_mag, test_dist) for f in freqs]
            @test length(results) == 3
            @test all(r -> r > 0 && isfinite(r), results)

            # ForwardDiff compatibility
            @test validate_fas_model(model)

            # Test gradient computation
            test_m = 6.5
            gradient = ForwardDiff.derivative(m -> compute_fas(model, test_freq, m, test_dist), test_m)
            @test isfinite(gradient)

            # Test with Dual numbers directly
            dual_mag = ForwardDiff.Dual(test_mag, 1.0)
            dual_result = compute_fas(model, test_freq, dual_mag, test_dist)
            @test dual_result isa ForwardDiff.Dual
            @test isfinite(ForwardDiff.value(dual_result))
            @test isfinite(ForwardDiff.partials(dual_result)[1])
        end

        @testset "FunctionalDurationModel" begin
            # Basic functionality
            simple_duration = (m, r, p) -> p[1] + p[2] * (m - 6.0) + p[3] * r
            params = [5.0, 0.5, 0.1]
            model = FunctionalDurationModel(simple_duration, params)

            @test model isa AbstractDurationModel
            @test model.func === simple_duration
            @test model.params == params

            # compute_duration
            result = compute_duration(model, test_mag, test_dist)
            @test result isa Real
            @test result > 0
            @test isfinite(result)

            # Test with different input types
            result_f32 = compute_duration(model, Float32(test_mag), Float32(test_dist))
            @test result_f32 isa Float32

            # Test parameter sensitivity
            mags = [5.0, 6.0, 7.0]
            durations = [compute_duration(model, m, test_dist) for m in mags]
            @test issorted(durations)  # Duration should increase with magnitude

            # ForwardDiff compatibility
            @test validate_duration_model(model)

            # Test gradient computation
            gradient = ForwardDiff.derivative(m -> compute_duration(model, m, test_dist), test_mag)
            @test isfinite(gradient)
            @test gradient > 0  # Duration should increase with magnitude

            # Test with Dual numbers directly
            dual_mag = ForwardDiff.Dual(test_mag, 1.0)
            dual_result = compute_duration(model, dual_mag, test_dist)
            @test dual_result isa ForwardDiff.Dual
            @test isfinite(ForwardDiff.value(dual_result))
            @test isfinite(ForwardDiff.partials(dual_result)[1])
        end

        @testset "CustomFASModel" begin
            # Create a custom FAS model type
            struct TestCustomFAS <: AbstractFASModel
                amplitude::Float64
                decay::Float64
            end

            function CustomModels.compute_fas(model::TestCustomFAS, freq, mag, dist)
                model.amplitude * freq^(-1.0) * exp(-model.decay * freq) * mag / dist
            end

            model = TestCustomFAS(100.0, 0.04)

            @test model isa AbstractFASModel

            # Test compute_fas
            result = compute_fas(model, test_freq, test_mag, test_dist)
            @test result isa Real
            @test result > 0
            @test isfinite(result)

            # Test type stability with Float32
            result_f32 = compute_fas(model, Float32(test_freq), Float32(test_mag), Float32(test_dist))
            @test result_f32 isa Float32

            # ForwardDiff compatibility
            @test validate_fas_model(model)

            # Test gradient
            gradient = ForwardDiff.derivative(m -> compute_fas(model, test_freq, m, test_dist), test_mag)
            @test isfinite(gradient)
        end

        @testset "CustomDurationModel" begin
            # Create a custom duration model type
            struct TestCustomDuration <: AbstractDurationModel
                base_duration::Float64
                mag_scaling::Float64
                dist_scaling::Float64
            end

            function CustomModels.compute_duration(model::TestCustomDuration, mag, dist)
                model.base_duration + model.mag_scaling * (mag - 6.0) + model.dist_scaling * dist
            end

            model = TestCustomDuration(5.0, 0.5, 0.1)

            @test model isa AbstractDurationModel

            # Test compute_duration
            result = compute_duration(model, test_mag, test_dist)
            @test result isa Real
            @test result > 0
            @test isfinite(result)

            # Test type stability
            result_f32 = compute_duration(model, Float32(test_mag), Float32(test_dist))
            @test result_f32 isa Float32

            # ForwardDiff compatibility
            @test validate_duration_model(model)

            # Test gradient
            gradient = ForwardDiff.derivative(m -> compute_duration(model, m, test_dist), test_mag)
            @test isfinite(gradient)
        end

        @testset "HybridFASModel" begin
            # Create component models
            model1_func = (f, m, r, p) -> p[1] * f^(-1.0)
            model1 = FunctionalFASModel(model1_func, [100.0])

            model2_func = (f, m, r, p) -> exp(-p[1] * f)
            model2 = FunctionalFASModel(model2_func, [0.04])

            # Test with function combining results
            combine_func = (results, f, m, r, p) -> results[1] * results[2]
            hybrid = HybridFASModel([model1, model2], combine_func, Float64[])

            @test hybrid isa AbstractFASModel
            @test length(hybrid.models) == 2

            # Test compute_fas
            result = compute_fas(hybrid, test_freq, test_mag, test_dist)
            @test result isa Real
            @test result > 0
            @test isfinite(result)

            # Verify it's actually combining the models
            result1 = compute_fas(model1, test_freq, test_mag, test_dist)
            result2 = compute_fas(model2, test_freq, test_mag, test_dist)
            @test result ≈ result1 * result2

            # Test with additional parameters in combine function
            combine_with_params = (results, f, m, r, p) -> p[1] * results[1] * results[2]
            hybrid_params = HybridFASModel([model1, model2], combine_with_params, [2.0])
            result_params = compute_fas(hybrid_params, test_freq, test_mag, test_dist)
            @test result_params ≈ 2.0 * result

            # ForwardDiff compatibility
            @test validate_fas_model(hybrid)
            @test validate_fas_model(hybrid_params)

            # Test gradient computation
            gradient = ForwardDiff.derivative(m -> compute_fas(hybrid, test_freq, m, test_dist), test_mag)
            @test isfinite(gradient)

            # Test with Dual numbers
            dual_mag = ForwardDiff.Dual(test_mag, 1.0)
            dual_result = compute_fas(hybrid, test_freq, dual_mag, test_dist)
            @test dual_result isa ForwardDiff.Dual
            @test isfinite(ForwardDiff.value(dual_result))
            @test isfinite(ForwardDiff.partials(dual_result)[1])

            # Test with three models
            model3_func = (f, m, r, p) -> m / r
            model3 = FunctionalFASModel(model3_func, Float64[])
            combine_three = (results, f, m, r, p) -> results[1] * results[2] * results[3]
            hybrid3 = HybridFASModel([model1, model2, model3], combine_three, Float64[])
            result3 = compute_fas(hybrid3, test_freq, test_mag, test_dist)
            @test isfinite(result3)
        end

        @testset "FourierParametersWrapper" begin
            # Create a simple FAS model
            fas_func = (f, m, r, p) -> p[1] * f^(-1.0)
            fas_model = FunctionalFASModel(fas_func, [100.0])

            dur_func = (m, r, p) -> p[1] + p[2] * m
            dur_model = FunctionalDurationModel(dur_func, [5.0, 0.5])

            # Create wrapper
            wrapper = FourierParametersWrapper(fas_model, dur_model)

            @test wrapper.fas_model === fas_model
            @test wrapper.duration_model === dur_model

            # Test compute_fas delegation
            fas_result = compute_fas(wrapper.fas_model, test_freq, test_mag, test_dist)
            @test isfinite(fas_result)

            # Test compute_duration delegation
            dur_result = compute_duration(wrapper.duration_model, test_mag, test_dist)
            @test isfinite(dur_result)
        end

        @testset "rvt_response_spectral_ordinate_custom" begin
            # Create models
            fas_func = (f, m, r, p) -> p[1] * f^(-p[2]) * exp(-p[3] * f) * 10^(p[4] * m) / r
            fas_model = FunctionalFASModel(fas_func, [100.0, 1.0, 0.04, 0.3])

            dur_func = (m, r, p) -> p[1] + p[2] * (m - 6.0) + p[3] * log10(r)
            dur_model = FunctionalDurationModel(dur_func, [5.0, 0.5, 1.0])

            # Test basic RVT computation
            Sa = rvt_response_spectral_ordinate_custom(test_period, test_mag, test_dist,
                fas_model, dur_model, test_damping)
            @test Sa isa Real
            @test Sa > 0
            @test isfinite(Sa)

            # Test with different periods
            periods = [0.1, 0.5, 1.0, 2.0]
            Sas = [rvt_response_spectral_ordinate_custom(T, test_mag, test_dist,
                fas_model, dur_model, test_damping)
                   for T in periods]
            @test length(Sas) == 4
            @test all(sa -> sa > 0 && isfinite(sa), Sas)

            # Test ForwardDiff compatibility
            gradient_mag = ForwardDiff.derivative(
                m -> rvt_response_spectral_ordinate_custom(test_period, m, test_dist,
                    fas_model, dur_model, test_damping),
                test_mag
            )
            @test isfinite(gradient_mag)

            gradient_dist = ForwardDiff.derivative(
                r -> rvt_response_spectral_ordinate_custom(test_period, test_mag, r,
                    fas_model, dur_model, test_damping),
                test_dist
            )
            @test isfinite(gradient_dist)

            # Test with custom duration model type
            struct SimpleCustomDuration <: AbstractDurationModel
                base::Float64
            end
            CustomModels.compute_duration(m::SimpleCustomDuration, mag, dist) = m.base + 0.5 * mag

            custom_dur = SimpleCustomDuration(5.0)
            Sa_custom = rvt_response_spectral_ordinate_custom(test_period, test_mag, test_dist,
                fas_model, custom_dur, test_damping)
            @test isfinite(Sa_custom)

            # Test with HybridFASModel
            model1 = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0), [100.0])
            model2 = FunctionalFASModel((f, m, r, p) -> exp(-p[1] * f), [0.04])
            hybrid = HybridFASModel([model1, model2],
                (results, f, m, r, p) -> results[1] * results[2] * 10^(p[1] * m) / r,
                [0.3])

            Sa_hybrid = rvt_response_spectral_ordinate_custom(test_period, test_mag, test_dist,
                hybrid, dur_model, test_damping)
            @test isfinite(Sa_hybrid)

            # Test gradient with hybrid model
            gradient_hybrid = ForwardDiff.derivative(
                m -> rvt_response_spectral_ordinate_custom(test_period, m, test_dist,
                    hybrid, dur_model, test_damping),
                test_mag
            )
            @test isfinite(gradient_hybrid)
        end

        @testset "Validation Functions" begin
            # Valid FAS model
            valid_fas = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0), [100.0])
            @test validate_fas_model(valid_fas) == true

            # Valid duration model
            valid_dur = FunctionalDurationModel((m, r, p) -> p[1] + p[2] * m, [5.0, 0.5])
            @test validate_duration_model(valid_dur) == true

            # Test that validation catches type instability issues
            # (if your implementation has this check)

            # Test validation with custom types
            struct ValidCustomFAS <: AbstractFASModel end
            CustomModels.compute_fas(::ValidCustomFAS, freq, mag, dist) = 100.0 * freq^(-1.0) * mag / dist

            valid_custom = ValidCustomFAS()
            @test validate_fas_model(valid_custom) == true

            struct ValidCustomDuration <: AbstractDurationModel end
            CustomModels.compute_duration(::ValidCustomDuration, mag, dist) = 5.0 + 0.5 * mag

            valid_custom_dur = ValidCustomDuration()
            @test validate_duration_model(valid_custom_dur) == true
        end

        @testset "Edge Cases and Error Handling" begin
            # Test with extreme parameter values
            fas_model = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0), [1e-6])
            result = compute_fas(fas_model, test_freq, test_mag, test_dist)
            @test isfinite(result)

            # Test with very high frequencies
            high_freq_result = compute_fas(fas_model, 100.0, test_mag, test_dist)
            @test isfinite(high_freq_result)

            # Test with very low frequencies
            low_freq_result = compute_fas(fas_model, 0.01, test_mag, test_dist)
            @test isfinite(low_freq_result)

            # Test with empty parameter arrays
            fas_no_params = FunctionalFASModel((f, m, r, p) -> f * m / r, Float64[])
            @test isfinite(compute_fas(fas_no_params, test_freq, test_mag, test_dist))

            dur_no_params = FunctionalDurationModel((m, r, p) -> m + r, Float64[])
            @test isfinite(compute_duration(dur_no_params, test_mag, test_dist))

            # Test parameter mutation doesn't affect model
            mutable_params = [100.0, 1.0]
            fas_mutable = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-p[2]), mutable_params)
            result1 = compute_fas(fas_mutable, test_freq, test_mag, test_dist)
            mutable_params[1] = 200.0  # Modify original array
            result2 = compute_fas(fas_mutable, test_freq, test_mag, test_dist)
            @test result1 == result2  # Should be equal if params are copied properly
        end

        @testset "Type Stability" begin
            # Test that Float32 inputs produce Float32 outputs
            fas_model = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0), [100.0])
            dur_model = FunctionalDurationModel((m, r, p) -> p[1] + p[2] * m, [5.0, 0.5])

            f32_freq = Float32(test_freq)
            f32_mag = Float32(test_mag)
            f32_dist = Float32(test_dist)

            fas_f32 = compute_fas(fas_model, f32_freq, f32_mag, f32_dist)
            @test fas_f32 isa Float32

            dur_f32 = compute_duration(dur_model, f32_mag, f32_dist)
            @test dur_f32 isa Float32

            # Note: RVT computation might not preserve Float32 due to internal calculations
            # but it should at least not error
            Sa_f32 = rvt_response_spectral_ordinate_custom(Float32(test_period), f32_mag, f32_dist,
                fas_model, dur_model, Float32(test_damping))
            @test isfinite(Sa_f32)
        end

        @testset "Multiple Gradient Computation" begin
            # Test computing gradients with respect to multiple variables
            fas_model = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0) * 10^(p[2] * m) / r,
                [100.0, 0.3])
            dur_model = FunctionalDurationModel((m, r, p) -> p[1] + p[2] * m + p[3] * log10(r),
                [5.0, 0.5, 1.0])

            # Gradient with respect to magnitude
            grad_m = ForwardDiff.derivative(
                m -> rvt_response_spectral_ordinate_custom(test_period, m, test_dist,
                    fas_model, dur_model, test_damping),
                test_mag
            )
            @test isfinite(grad_m)

            # Gradient with respect to distance
            grad_r = ForwardDiff.derivative(
                r -> rvt_response_spectral_ordinate_custom(test_period, test_mag, r,
                    fas_model, dur_model, test_damping),
                test_dist
            )
            @test isfinite(grad_r)

            # Hessian computation (if needed for optimization)
            hess_m = ForwardDiff.derivative(
                m -> ForwardDiff.derivative(
                    m2 -> rvt_response_spectral_ordinate_custom(test_period, m2, test_dist,
                        fas_model, dur_model, test_damping),
                    m
                ),
                test_mag
            )
            @test isfinite(hess_m)
        end
    end

end
