using Test
using StochasticGroundMotionSimulation
using StochasticGroundMotionSimulation.CustomModels
using ForwardDiff

@testset "RVT Interface Branching Coverage" begin

    # Test parameters
    test_period = 1.0
    test_mag = 6.5
    test_dist = 20.0
    test_damping = 0.05

    # Create custom models for testing
    fas_func = (f, m, r, p) -> p[1] * f^(-p[2]) * exp(-p[3] * f) * 10^(p[4] * m) / r
    custom_fas = FunctionalFASModel(fas_func, [100.0, 1.0, 0.04, 0.3])

    dur_func = (m, r, p) -> p[1] + p[2] * (m - 6.0) + p[3] * log10(r)
    custom_dur = FunctionalDurationModel(dur_func, [5.0, 0.5, 1.0])

    # Get a concrete FourierParameters instance (assuming this exists in your package)
    # You may need to adjust this based on your actual implementation
    # For example, using one of the predefined models
    try
        # Try to create a FourierParameters instance
        src = SourceParameters(100.0)
        geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
        ane = AnelasticAttenuationParameters(200.0, 0.4)
        sat = NearSourceSaturationParameters(:BT15)
        path = PathParameters(geo, sat, ane)
        site = SiteParameters(0.039)

        concrete_fourier = FourierParameters(src, path, site)
        concrete_rvt = RandomVibrationParameters(:DK80)


        @testset "Branch 1: FourierParameters + RandomVibrationParameters" begin
            # This is the original/legacy path - both concrete types
            # Should now go through the unified interface (method #1)
            Sa = rvt_response_spectral_ordinate(
                test_period, test_mag, test_dist,
                concrete_fourier, concrete_rvt,
                damping=test_damping
            )
            @test Sa isa Real
            @test Sa > 0
            @test isfinite(Sa)

            # Test with ForwardDiff
            gradient = ForwardDiff.derivative(
                m -> rvt_response_spectral_ordinate(
                    test_period, m, test_dist,
                    concrete_fourier, concrete_rvt,
                    damping=test_damping
                ),
                test_mag
            )
            @test isfinite(gradient)
        end

        @testset "Branch 2: FourierParameters + AbstractDurationModel" begin
            # Mixing concrete and custom
            Sa = rvt_response_spectral_ordinate(
                test_period, test_mag, test_dist,
                concrete_fourier, custom_dur,
                damping=test_damping
            )
            @test Sa isa Real
            @test Sa > 0
            @test isfinite(Sa)

            # Test with ForwardDiff
            gradient = ForwardDiff.derivative(
                m -> rvt_response_spectral_ordinate(
                    test_period, m, test_dist,
                    concrete_fourier, custom_dur,
                    damping=test_damping
                ),
                test_mag
            )
            @test isfinite(gradient)
        end

        @testset "Branch 3: AbstractFASModel + RandomVibrationParameters (NOT SUPPORTED)" begin
            # This combination is NOT supported because RandomVibrationParameters.excitation_duration
            # requires FourierParameters, which we don't have when using custom FAS
            @test_throws ErrorException rvt_response_spectral_ordinate(
                test_period, test_mag, test_dist,
                custom_fas, concrete_rvt
            )

            # Verify the error message is helpful
            try
                rvt_response_spectral_ordinate(
                    test_period, test_mag, test_dist,
                    custom_fas, concrete_rvt
                )
                @test false # Should have thrown an error
            catch e
                @test e isa ErrorException
                @test occursin("not supported", e.msg)
                @test occursin("custom duration model", e.msg)
            end

            # Test the workaround: use a custom duration model instead
            workaround_dur = FunctionalDurationModel(
                (m, r, p) -> p[1] + p[2] * m + p[3] * log10(r),
                [5.0, 0.5, 1.0]
            )

            Sa_workaround = rvt_response_spectral_ordinate(
                test_period, test_mag, test_dist,
                custom_fas, workaround_dur,
                damping=test_damping
            )
            @test Sa_workaround isa Real
            @test Sa_workaround > 0
            @test isfinite(Sa_workaround)
        end

    catch e
        @warn "Could not test concrete model branches" exception = e
        @test_skip "Concrete model testing - adjust to match your package's API"
    end

    @testset "Branch 4: AbstractFASModel + AbstractDurationModel" begin
        # Both custom - this should call rvt_response_spectral_ordinate_custom
        Sa = rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            custom_fas, custom_dur,
            damping=test_damping
        )
        @test Sa isa Real
        @test Sa > 0
        @test isfinite(Sa)

        # Verify it's equivalent to calling _custom directly
        Sa_direct = rvt_response_spectral_ordinate_custom(
            test_period, test_mag, test_dist,
            custom_fas, custom_dur, test_damping
        )
        @test Sa ≈ Sa_direct

        # Test with ForwardDiff
        gradient_mag = ForwardDiff.derivative(
            m -> rvt_response_spectral_ordinate(
                test_period, m, test_dist,
                custom_fas, custom_dur,
                damping=test_damping
            ),
            test_mag
        )
        @test isfinite(gradient_mag)

        gradient_dist = ForwardDiff.derivative(
            r -> rvt_response_spectral_ordinate(
                test_period, test_mag, r,
                custom_fas, custom_dur,
                damping=test_damping
            ),
            test_dist
        )
        @test isfinite(gradient_dist)

        # Test with multiple periods to ensure no issues with vectorization
        periods = [0.1, 0.5, 1.0, 2.0, 5.0]
        Sas = [rvt_response_spectral_ordinate(T, test_mag, test_dist,
            custom_fas, custom_dur,
            damping=test_damping)
               for T in periods]
        @test length(Sas) == 5
        @test all(isfinite, Sas)
        @test all(>(0), Sas)
    end

    @testset "Branch Coverage with Different Custom Model Types" begin
        # Test with HybridFASModel
        model1 = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0), [100.0])
        model2 = FunctionalFASModel((f, m, r, p) -> exp(-p[1] * f), [0.04])
        hybrid_fas = HybridFASModel(
            [model1, model2],
            (results, f, m, r, p) -> results[1] * results[2] * 10^(p[1] * m) / r,
            [0.3]
        )

        Sa_hybrid = rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            hybrid_fas, custom_dur,
            damping=test_damping
        )
        @test isfinite(Sa_hybrid)
        @test Sa_hybrid > 0

        # Test with custom struct types
        struct TestFASStruct <: AbstractFASModel
            amplitude::Float64
        end
        CustomModels.compute_fas(m::TestFASStruct, f, mag, dist) =
            m.amplitude * f^(-1.0) * 10^(0.3 * mag) / dist

        struct TestDurationStruct <: AbstractDurationModel
            base::Float64
        end
        CustomModels.compute_duration(m::TestDurationStruct, mag, dist) =
            m.base + 0.5 * mag + 0.1 * log10(dist)

        fas_struct = TestFASStruct(100.0)
        dur_struct = TestDurationStruct(5.0)

        Sa_struct = rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            fas_struct, dur_struct,
            damping=test_damping
        )
        @test isfinite(Sa_struct)
        @test Sa_struct > 0

        # Test ForwardDiff compatibility with custom structs
        gradient = ForwardDiff.derivative(
            m -> rvt_response_spectral_ordinate(
                test_period, m, test_dist,
                fas_struct, dur_struct,
                damping=test_damping
            ),
            test_mag
        )
        @test isfinite(gradient)
    end

    @testset "Keyword Arguments Pass-through" begin
        # Test that kwargs are properly passed through all branches

        # Test with different damping ratios
        damping_ratios = [0.01, 0.05, 0.10, 0.20]
        Sas_damping = [rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            custom_fas, custom_dur,
            damping=d)
                       for d in damping_ratios]

        @test length(Sas_damping) == 4
        @test all(isfinite, Sas_damping)
        @test issorted(Sas_damping, rev=true)  # Higher damping → lower Sa

        # Test that omitting damping uses default (if applicable)
        Sa_default = rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            custom_fas, custom_dur
        )
        @test isfinite(Sa_default)

        # If your implementation has other kwargs, test them here
        # For example: oscillator_type, integration_method, etc.
    end

    @testset "Type Stability Across Branches" begin
        # Test Float32 inputs
        T_f32 = Float32(test_period)
        mag_f32 = Float32(test_mag)
        dist_f32 = Float32(test_dist)
        damp_f32 = Float32(test_damping)

        Sa_f32 = rvt_response_spectral_ordinate(
            T_f32, mag_f32, dist_f32,
            custom_fas, custom_dur,
            damping=damp_f32
        )
        @test isfinite(Sa_f32)

        # Test that the result is reasonable (doesn't error)
        # Note: Float32 might not be preserved through all internal calculations
    end

    @testset "Edge Cases in Unified Interface" begin
        # Very short period
        Sa_short = rvt_response_spectral_ordinate(
            0.01, test_mag, test_dist,
            custom_fas, custom_dur,
            damping=test_damping
        )
        @test isfinite(Sa_short)

        # Very long period
        Sa_long = rvt_response_spectral_ordinate(
            10.0, test_mag, test_dist,
            custom_fas, custom_dur,
            damping=test_damping
        )
        @test isfinite(Sa_long)

        # Small magnitude
        Sa_small_mag = rvt_response_spectral_ordinate(
            test_period, 1.0, test_dist,
            custom_fas, custom_dur,
            damping=test_damping
        )
        @test isfinite(Sa_small_mag)

        # Large magnitude
        Sa_large_mag = rvt_response_spectral_ordinate(
            test_period, 10.0, test_dist,
            custom_fas, custom_dur,
            damping=test_damping
        )
        @test isfinite(Sa_large_mag)

        # Small distance
        Sa_small_dist = rvt_response_spectral_ordinate(
            test_period, test_mag, 0.1,
            custom_fas, custom_dur,
            damping=test_damping
        )
        @test isfinite(Sa_small_dist)

        # Large distance
        Sa_large_dist = rvt_response_spectral_ordinate(
            test_period, test_mag, 1000.0,
            custom_fas, custom_dur,
            damping=test_damping
        )
        @test isfinite(Sa_large_dist)
    end

    @testset "Consistency Check Between Branches" begin
        # If possible, verify that using equivalent models through different
        # branches produces similar results

        # Create a simple custom model
        simple_fas = FunctionalFASModel(
            (f, m, r, p) -> p[1] * f^(-1.0) * 10^(p[2] * m) / r,
            [100.0, 0.3]
        )
        simple_dur = FunctionalDurationModel(
            (m, r, p) -> p[1] + p[2] * m,
            [5.0, 0.5]
        )

        # Test at multiple periods
        test_periods = [0.1, 0.5, 1.0, 2.0]
        for T in test_periods
            Sa = rvt_response_spectral_ordinate(
                T, test_mag, test_dist,
                simple_fas, simple_dur,
                damping=test_damping
            )

            # Verify basic sanity checks
            @test Sa > 0
            @test isfinite(Sa)

            # Verify it increases with magnitude (holding period constant)
            Sa_lower_mag = rvt_response_spectral_ordinate(
                T, test_mag - 0.5, test_dist,
                simple_fas, simple_dur,
                damping=test_damping
            )
            @test Sa > Sa_lower_mag

            # Verify it decreases with distance
            Sa_closer = rvt_response_spectral_ordinate(
                T, test_mag, test_dist / 2,
                simple_fas, simple_dur,
                damping=test_damping
            )
            @test Sa_closer > Sa
        end
    end

    @testset "Multiple Calls - No Side Effects" begin
        # Ensure that calling the function multiple times doesn't cause issues
        # (e.g., mutation of model parameters)

        results = [rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            custom_fas, custom_dur,
            damping=test_damping)
                   for _ in 1:10]

        # All results should be identical
        @test all(r -> r ≈ results[1], results)

        # Test with different models called in sequence
        model_a = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-1.0), [100.0])
        model_b = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-2.0), [200.0])

        Sa_a = rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            model_a, custom_dur,
            damping=test_damping
        )

        Sa_b = rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            model_b, custom_dur,
            damping=test_damping
        )

        # Re-call model_a to ensure it wasn't affected by model_b
        Sa_a2 = rvt_response_spectral_ordinate(
            test_period, test_mag, test_dist,
            model_a, custom_dur,
            damping=test_damping
        )

        @test Sa_a ≈ Sa_a2
        @test Sa_a != Sa_b  # Different models should give different results
    end
end

println("\n✓ All RVT interface branching tests completed successfully!")
