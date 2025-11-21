using Test
using StochasticGroundMotionSimulation
using StochasticGroundMotionSimulation.CustomModels
using ForwardDiff

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
        @test typeof(result_f32) <: AbstractFloat

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

    @testset "FunctionalFASModel functor call" begin
        fas_func = (f, m, r, p) -> p[1] * f^(-p[2])
        fas_model = FunctionalFASModel(fas_func, [100.0, 1.0])

        # Direct compute_fas
        fas_val = compute_fas(fas_model, test_freq, test_mag, test_dist)

        # Functor call syntax: model(f, m, r)
        fas_call = fas_model(test_freq, test_mag, test_dist)

        @test fas_call ≈ fas_val
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
        @test typeof(result_f32) <: AbstractFloat

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

    @testset "FunctionalDurationModel functor call" begin
        dur_func = (m, r, p) -> p[1] + p[2] * m + p[3] * r
        dur_model = FunctionalDurationModel(dur_func, [5.0, 0.5, 0.1])

        dur_val = compute_duration(dur_model, test_mag, test_dist)
        dur_call = dur_model(test_mag, test_dist)

        @test dur_call ≈ dur_val

        dur_call_int = dur_model(Int(5), Int(10))
        @test dur_call_int isa Real
    end

    @testset "Functional models parameter type validation" begin
        # Non-Real element type → should throw for FAS
        fas_func = (f, m, r, p) -> p[1] * f
        @test_throws ArgumentError FunctionalFASModel(fas_func, ComplexF64[1+1im])

        # Another non-Real element type → should throw for Duration
        dur_func = (m, r, p) -> p[1]
        @test_throws ArgumentError FunctionalDurationModel(dur_func, ["not", "real"])
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
        @test typeof(result_f32) <: AbstractFloat

        # ForwardDiff compatibility
        @test validate_fas_model(model)

        # Test gradient
        gradient = ForwardDiff.derivative(m -> compute_fas(model, test_freq, m, test_dist), test_mag)
        @test isfinite(gradient)
    end

    @testset "Built-in CustomFASModel (PJScustomModels.jl)" begin
        # Valid model
        custom_fas = CustomFASModel(50.0, 200.0, 0.3, 0.04, 760.0)

        @test custom_fas isa AbstractFASModel

        # r <= 50 branch
        fas_near = compute_fas(custom_fas, test_freq, test_mag, 10.0)
        # r > 50 branch
        fas_far = compute_fas(custom_fas, test_freq, test_mag, 100.0)

        @test fas_near > fas_far > 0.0
        @test isfinite(fas_near)
        @test isfinite(fas_far)

        # Type promotion: Float32 inputs
        fas_f32 = compute_fas(
            custom_fas,
            Float32(test_freq),
            Float32(test_mag),
            Float32(test_dist),
        )
        @test typeof(fas_f32) <: AbstractFloat

        # ForwardDiff compatibility and gradient
        @test validate_fas_model(custom_fas)

        grad_m = ForwardDiff.derivative(
            m -> compute_fas(custom_fas, test_freq, m, test_dist),
            test_mag,
        )
        grad_r = ForwardDiff.derivative(
            r -> compute_fas(custom_fas, test_freq, test_mag, r),
            test_dist,
        )

        @test isfinite(grad_m)
        @test isfinite(grad_r)

        # Mixed-type constructor uses `promote`
        mixed_fas = CustomFASModel(50, 200.0, 0.3f0, 0.04, 760)
        @test mixed_fas isa CustomFASModel

        # Constructor validation branches
        @test_throws ArgumentError CustomFASModel(-1.0, 200.0, 0.3, 0.04, 760.0)   # stress_drop
        @test_throws ArgumentError CustomFASModel(50.0, 0.0, 0.3, 0.04, 760.0)     # Q0
        @test_throws ArgumentError CustomFASModel(50.0, 200.0, 0.3, -0.01, 760.0)  # kappa
        @test_throws ArgumentError CustomFASModel(50.0, 200.0, 0.3, 0.04, 0.0)     # vs30
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
        @test typeof(result_f32) <: AbstractFloat

        # ForwardDiff compatibility
        @test validate_duration_model(model)

        # Test gradient
        gradient = ForwardDiff.derivative(m -> compute_duration(model, m, test_dist), test_mag)
        @test isfinite(gradient)
    end

    @testset "Built-in CustomDurationModel (PJScustomModels.jl)" begin
        custom_dur = CustomDurationModel(0.1, 0.05, 0.1, 2.0)

        @test custom_dur isa AbstractDurationModel

        dur = compute_duration(custom_dur, test_mag, test_dist)
        @test dur > 0
        @test isfinite(dur)

        # Check dependence on distance (path term)
        dur_near = compute_duration(custom_dur, test_mag, test_dist / 2)
        dur_far = compute_duration(custom_dur, test_mag, test_dist * 2)
        @test dur_far > dur_near

        # Float32 inputs
        dur_f32 = compute_duration(
            custom_dur,
            Float32(test_mag),
            Float32(test_dist),
        )
        @test typeof(dur_f32) <: AbstractFloat

        # Validation and gradients
        @test validate_duration_model(custom_dur)

        grad_m = ForwardDiff.derivative(
            m -> compute_duration(custom_dur, m, test_dist),
            test_mag,
        )
        grad_r = ForwardDiff.derivative(
            r -> compute_duration(custom_dur, test_mag, r),
            test_dist,
        )

        @test isfinite(grad_m)
        @test isfinite(grad_r)

        # Mixed-type constructor uses `promote`
        mixed_dur = CustomDurationModel(1, 0.5f0, 0.1, 2)
        @test mixed_dur isa CustomDurationModel
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

    @testset "Delegating wrappers to existing spectrum / duration" begin
        SGS = StochasticGroundMotionSimulation

        # --- Build a simple FourierParameters instance ---

        # Source: just stress drop – uses default radiation pattern, etc.
        src = SGS.SourceParameters(100.0)

        # Geometric spreading: piecewise with a single segment (very simple)
        geo = SGS.GeometricSpreadingParameters([1.0, Inf], [1.0])

        # Saturation: no saturation
        sat = SGS.NearSourceSaturationParameters(:BT15)

        # Anelastic attenuation: simple Q0, η 
        anel = SGS.AnelasticAttenuationParameters(
            300.0,  # Q0 (constrained or variable, depending on constructor)
            0.0,     # η
            3.5 # c\_Q    
        )

        path = SGS.PathParameters(geo, sat, anel)

        site = SGS.SiteParameters(0.04, 0.0, 0.0, SGS.SiteAmpUnit())

        fas_params = SGS.FourierParameters(src, path, site)

        # Random vibration parameters (default DK80 / BT15 / BT15 / ACR)
        rvt_params = SGS.RandomVibrationParameters()

        # --- FourierParametersWrapper: compute_fas delegation ---

        fp_wrapper = FourierParametersWrapper(fas_params)
        fp_val = compute_fas(fp_wrapper, test_freq, test_mag, test_dist)
        @test isfinite(fp_val)

        # --- ExistingDurationWrapper: compute_duration delegation ---

        dur_wrapper = ExistingDurationWrapper(rvt_params, fas_params)
        dur_val = compute_duration(dur_wrapper, test_mag, test_dist)
        @test isfinite(dur_val)
        @test dur_val > 0
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

    @testset "rvt_response_spectral_ordinate_custom - zero FAS fallback branch" begin
        # FAS model that returns exactly zero at all frequencies
        zero_fas = FunctionalFASModel(
            (f, m, r, p) -> zero(eltype(p)),
            Float64[],   # params unused
        )

        # Simple constant duration
        const_dur = FunctionalDurationModel(
            (m, r, p) -> 10.0,
            Float64[],
        )

        Sa_zero = rvt_response_spectral_ordinate_custom(
            test_period, test_mag, test_dist,
            zero_fas, const_dur, test_damping;
            freq_range=(0.5, 2.0),  # non-default range
            nfreq=5,                # non-default number of frequencies
        )

        # With zero FAS everywhere, arms = 0, so spectral ordinate should be 0
        @test Sa_zero == 0.0
    end
   
    @testset "HybridFASModel constructor validation" begin
        combine_func = (results, f, m, r, p) -> results[1]

        # Empty model list should throw
        @test_throws ArgumentError HybridFASModel(
            AbstractFASModel[],  # empty
            combine_func,
            Float64[],
        )
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

    # @testset "Type Stability" begin
    #     # Test that Float32 inputs produce Float32 outputs
    #     fas_model = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0), [100.0])
    #     dur_model = FunctionalDurationModel((m, r, p) -> p[1] + p[2] * m, [5.0, 0.5])

    #     f32_freq = Float32(test_freq)
    #     f32_mag = Float32(test_mag)
    #     f32_dist = Float32(test_dist)

    #     fas_f32 = compute_fas(fas_model, f32_freq, f32_mag, f32_dist)
    #     @test typeof(fas_f32) <: AbstractFloat

    #     dur_f32 = compute_duration(dur_model, f32_mag, f32_dist)
    #     @test typeof(dur_f32) <: AbstractFloat

    #     # Note: RVT computation might not preserve Float32 due to internal calculations
    #     # but it should at least not error
    #     Sa_f32 = rvt_response_spectral_ordinate_custom(Float32(test_period), f32_mag, f32_dist,
    #         fas_model, dur_model, Float32(test_damping))
    #     @test isfinite(Sa_f32)
    # end

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

println("\n✓ All custom models tests completed successfully!")
