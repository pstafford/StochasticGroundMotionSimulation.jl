module CustomModels

using ForwardDiff

# =====================================
# Abstract Types (imported from parent module)
# =====================================

# Import abstract types defined in parent module
using ..StochasticGroundMotionSimulation: AbstractFASModel, AbstractDurationModel


# =====================================
# Function Interfaces
# =====================================

"""
    compute_fas(model::AbstractFASModel, freq, mag, dist)

Compute the Fourier Amplitude Spectrum at a given frequency.

# Arguments
- `model`: FAS model instance containing model parameters
- `freq`: Frequency in Hz (can be Real or ForwardDiff.Dual)
- `mag`: Earthquake magnitude (can be Real or ForwardDiff.Dual)
- `dist`: Distance in km (can be Real or ForwardDiff.Dual)

# Returns
- FAS value at the given frequency

# Notes
- This function must be differentiable for ForwardDiff.jl compatibility
- Avoid branching logic that depends on the values of freq, mag, or dist
- Use smooth approximations for discontinuous functions
"""
function compute_fas end

"""
    compute_duration(model::AbstractDurationModel, mag, dist)

Compute the duration metric.

# Arguments
- `model`: Duration model instance containing model parameters
- `mag`: Earthquake magnitude (can be Real or ForwardDiff.Dual)
- `dist`: Distance in km (can be Real or ForwardDiff.Dual)

# Returns
- Duration value in seconds

# Notes
- This function must be differentiable for ForwardDiff.jl compatibility
"""
function compute_duration end

# =====================================
# Functional Interface (for quick prototyping)
# =====================================

"""
    FunctionalFASModel{F}

A wrapper that allows users to pass a function as a FAS model.
Useful for quick prototyping without defining a new type.

# Example
```julia
# Define a simple FAS function
my_fas = (f, m, r, params) -> begin
    stress_drop, Q0, kappa = params
    corner_freq = 4.9e6 * 1e-6 * (stress_drop/1e6)^(1/3) * 10^(-m/3)
    source = 1.0 / (1.0 + (f/corner_freq)^2)
    path = exp(-π * f * r / (Q0 * 30.0))  # simplified
    site = exp(-π * kappa * f)
    return source * path * site
end

# Create model with parameters
model = FunctionalFASModel(my_fas, [100.0, 200.0, 0.04])

# Use in computation
fas_value = compute_fas(model, 1.0, 6.0, 10.0)
```
"""
struct FunctionalFASModel{F,P} <: AbstractFASModel
    func::F
    params::P

    function FunctionalFASModel(func::F, params::P) where {F,P}
        # Validate that params can be used with ForwardDiff
        if !(eltype(params) <: Real)
            throw(ArgumentError("Parameters must be Real-valued for ForwardDiff compatibility"))
        end
        new{F,P}(func, copy(params))
    end
end

function compute_fas(model::FunctionalFASModel, freq, mag, dist)
    # Promote all inputs to common type for ForwardDiff
    # Handle both regular numbers and ForwardDiff.Dual types
    T = promote_type(typeof(freq), typeof(mag), typeof(dist), eltype(model.params))
    return model.func(
        convert(T, freq),
        convert(T, mag),
        convert(T, dist),
        convert(Vector{T}, model.params)
    )
end

# Make FunctionalFASModel callable directly
function (model::FunctionalFASModel)(freq::Real, mag::Real, dist::Real)
    return compute_fas(model, freq, mag, dist)
end


"""
    FunctionalDurationModel{F}

A wrapper that allows users to pass a function as a duration model.

# Example
```julia
# Define a simple duration function
my_duration = (m, r, params) -> begin
    a, b, c = params
    return a + b * m + c * log10(r)
end

# Create model with parameters
model = FunctionalDurationModel(my_duration, [1.0, 0.5, 0.3])

# Use in computation
duration = compute_duration(model, 6.0, 10.0)
```
"""
struct FunctionalDurationModel{F,P} <: AbstractDurationModel
    func::F
    params::P

    function FunctionalDurationModel(func::F, params::P) where {F,P}
        if !(eltype(params) <: Real)
            throw(ArgumentError("Parameters must be Real-valued for ForwardDiff compatibility"))
        end
        new{F,P}(func, copy(params))
    end
end

function compute_duration(model::FunctionalDurationModel, mag, dist)
    # Promote all inputs to common type for ForwardDiff
    # Handle both regular numbers and ForwardDiff.Dual types
    T = promote_type(typeof(mag), typeof(dist), eltype(model.params))
    return model.func(
        convert(T, mag),
        convert(T, dist),
        convert(Vector{T}, model.params)
    )
end

# Make FunctionalDurationModel callable directly
function (model::FunctionalDurationModel)(mag::Real, dist::Real)
    return compute_duration(model, mag, dist)
end


# =====================================
# Custom Type Examples
# =====================================

"""
    CustomFASModel

Example of a custom FAS model with structured parameters.
This approach is preferred for more complex models.
"""
struct CustomFASModel{T<:Real} <: AbstractFASModel
    stress_drop::T
    Q0::T
    eta::T
    kappa::T
    vs30::T

    # Optional: Add validation in constructor
    function CustomFASModel(stress_drop::T, Q0::T, eta::T, kappa::T, vs30::T) where T<:Real
        stress_drop > 0 || throw(ArgumentError("stress_drop must be positive"))
        Q0 > 0 || throw(ArgumentError("Q0 must be positive"))
        kappa >= 0 || throw(ArgumentError("kappa must be non-negative"))
        vs30 > 0 || throw(ArgumentError("vs30 must be positive"))
        new{T}(stress_drop, Q0, eta, kappa, vs30)
    end
end

# Constructor for mixed types (promotes to common type)
CustomFASModel(stress_drop, Q0, eta, kappa, vs30) =
    CustomFASModel(promote(stress_drop, Q0, eta, kappa, vs30)...)

function compute_fas(model::CustomFASModel, freq, mag, dist)
    # Promote all to common type for ForwardDiff compatibility
    T = promote_type(typeof(freq), typeof(mag), typeof(dist), typeof(model.stress_drop))
    f = convert(T, freq)
    m = convert(T, mag)
    r = convert(T, dist)

    # Source spectrum (Brune model)
    M0 = 10^(1.5 * m + 16.05)  # Moment in dyne-cm
    corner_freq = 4.9e6 * model.vs30 * (model.stress_drop / (M0^(1 / 3)))
    source_spectrum = (2π * f)^2 / (1 + (f / corner_freq)^2)

    # Path effects
    Q = model.Q0 * f^model.eta
    path_atten = exp(-π * f * r / (Q * 3.5))  # assuming vs = 3.5 km/s

    # Geometric spreading (simplified)
    if r <= 50
        geom_spread = 1 / r
    else
        geom_spread = 1 / (50 * sqrt(r / 50))
    end

    # Site effects
    site_atten = exp(-π * model.kappa * f)

    # Combine all effects
    return source_spectrum * path_atten * geom_spread * site_atten
end

"""
    CustomDurationModel

Example of a custom duration model.
"""
struct CustomDurationModel{T<:Real} <: AbstractDurationModel
    a::T
    b::T
    c::T
    d::T
end

CustomDurationModel(a, b, c, d) = CustomDurationModel(promote(a, b, c, d)...)

function compute_duration(model::CustomDurationModel, mag, dist)
    T = promote_type(typeof(mag), typeof(dist), typeof(model.a))
    m = convert(T, mag)
    r = convert(T, dist)

    # Example duration model
    source_duration = model.a * (10^(model.b * m))
    path_duration = model.c * r

    # Use smooth maximum to avoid discontinuities
    # smooth_max(x, y) ≈ max(x, y) but differentiable
    α = 10.0  # smoothing parameter
    total_duration = (1 / α) * log(exp(α * source_duration) + exp(α * path_duration))

    return total_duration + model.d
end

# =====================================
# Integration with RVT Calculations
# =====================================

"""
    rvt_response_spectral_ordinate_custom(
        T::Real,
        mag::Real,
        dist::Real,
        fas_model::AbstractFASModel,
        duration_model::AbstractDurationModel,
        damping::Real=0.05,
        glxi=nothing,
        glwi=nothing
    )

Compute response spectral ordinate using custom FAS and duration models.

# Arguments
- `T`: Oscillator period (seconds)
- `mag`: Earthquake magnitude
- `dist`: Distance (km)
- `fas_model`: Custom FAS model
- `duration_model`: Custom duration model
- `damping`: Damping ratio (default: 0.05)
- `glxi`: Gauss-Legendre quadrature points (optional, for compatibility)
- `glwi`: Gauss-Legendre quadrature weights (optional, for compatibility)

# Keyword Arguments (for backward compatibility)
- `freq_range`: Frequency range for integration (Hz)
- `nfreq`: Number of frequency points for integration

# Returns
- Response spectral acceleration (g)

# Notes
- The `glxi` and `glwi` parameters are included for API compatibility with the
  traditional RVT interface but are not currently used in this implementation,
  which uses log-spaced frequency points instead.
"""
function rvt_response_spectral_ordinate_custom(
    T::Real,
    mag::Real,
    dist::Real,
    fas_model::AbstractFASModel,
    duration_model::AbstractDurationModel,
    damping::Real=0.05,
    glxi=nothing,
    glwi=nothing;
    freq_range=(0.01, 100.0),
    nfreq=200
)
    # Natural frequency of oscillator
    fn = 1 / T

    # Get duration
    duration = compute_duration(duration_model, mag, dist)

    # Frequency array for integration (log-spaced)
    freq_min, freq_max = freq_range
    freqs = exp10.(range(log10(freq_min), log10(freq_max), length=nfreq))

    # Compute transfer function and FAS at each frequency
    m0 = zero(promote_type(typeof(mag), typeof(dist)))
    m1 = zero(m0)
    m2 = zero(m0)

    for f in freqs
        # Transfer function for SDOF oscillator
        freq_ratio = f / fn
        H_squared = 1 / ((1 - freq_ratio^2)^2 + (2 * damping * freq_ratio)^2)

        # FAS at this frequency
        fas = compute_fas(fas_model, f, mag, dist)
        fas_squared = fas^2

        # Moments of the process
        m0 += fas_squared * H_squared
        m1 += fas_squared * H_squared * f
        m2 += fas_squared * H_squared * f^2
    end

    # Normalize by frequency spacing (trapezoidal rule approximation)
    df = (freq_max - freq_min) / (nfreq - 1)
    m0 *= df
    m1 *= df * 2π
    m2 *= df * (2π)^2

    # RMS acceleration
    arms = sqrt(m0 / duration)

    # Peak factor (Vanmarcke approximation)
    if m0 > 0 && m2 > 0
        # Compute bandwidth parameter with bounds checking to avoid numerical issues
        delta_squared = 1 - m1^2 / (m0 * m2)

        # Clamp to valid range [0, 1] to handle numerical precision issues
        delta_squared = max(0.0, min(1.0, delta_squared))
        delta = sqrt(delta_squared)

        # Avoid division by zero when delta is very small
        if delta > 1e-10
            fz = sqrt(m2 / m0) / (2π)
            Ne = duration * fz / delta

            # Ne must be > 1 for log(Ne) to be positive
            if Ne > 1.0
                eta = sqrt(2 * log(Ne))
                peak_factor = eta + 0.5772 / eta
            else
                peak_factor = 2.5  # fallback value
            end
        else
            peak_factor = 2.5  # fallback for narrow-band processes
        end
    else
        peak_factor = 2.5  # fallback value
    end

    # Peak acceleration in g
    peak_accel = peak_factor * arms / 980.665  # convert cm/s² to g

    return peak_accel
end

# =====================================
# Hybrid Models (mix custom and built-in)
# =====================================

"""
    HybridFASModel

Combines multiple FAS models with a custom combining function.
Useful for mixing different model components or creating ensemble models.

# Fields
- `models`: Vector of AbstractFASModel instances to combine
- `combine_func`: Function that combines results: (results, freq, mag, dist, params) -> combined_result
- `params`: Additional parameters for the combining function

# Example
```julia
# Create component models
source_model = FunctionalFASModel((f, m, r, p) -> p[1] * f^(-2), [1e20])
path_model = FunctionalFASModel((f, m, r, p) -> exp(-π*f*r/(p[1]*3.5)), [200.0])

# Combine by multiplication
combine = (results, f, m, r, p) -> results[1] * results[2]
hybrid = HybridFASModel([source_model, path_model], combine, Float64[])

# Use like any other FAS model
fas = compute_fas(hybrid, 1.0, 6.0, 10.0)
```
"""
struct HybridFASModel{M,F,P} <: AbstractFASModel
    models::M  # Vector of models
    combine_func::F
    params::P

    function HybridFASModel(models::Vector{<:AbstractFASModel}, combine_func::F, params::P) where {F,P}
        length(models) >= 1 || throw(ArgumentError("HybridFASModel requires at least one component model"))
        new{typeof(models),F,P}(models, combine_func, params)
    end
end

function compute_fas(model::HybridFASModel, freq, mag, dist)
    # Compute FAS for each component model
    results = [compute_fas(m, freq, mag, dist) for m in model.models]

    # Promote types for ForwardDiff compatibility
    T = promote_type(typeof(freq), typeof(mag), typeof(dist), eltype(model.params))

    # Combine results using the provided function
    return model.combine_func(
        results,
        convert(T, freq),
        convert(T, mag),
        convert(T, dist),
        convert(Vector{T}, model.params)
    )
end

# =====================================
# Bridge to Traditional Parameters
# =====================================

# Import FourierParameters and related functions from parent module
using ..StochasticGroundMotionSimulation: FourierParameters, fourier_spectral_ordinate, RandomVibrationParameters, excitation_duration

"""
    FourierParametersWrapper <: AbstractFASModel

Wrapper to make existing FourierParameters compatible with AbstractFASModel interface.
This allows using traditional parameter-based models with the custom model interface.

# Example
```julia
# Create traditional parameters
src = SourceParameters(100.0)
geo = GeometricSpreadingParameters([1.0, 50.0, Inf], [1.0, 0.5])
ane = AnelasticAttenuationParameters(200.0, 0.4)
sat = NearSourceSaturationParameters(:BT15)
path = PathParameters(geo, sat, ane)
site = SiteParameters(0.039)
fourier_params = FourierParameters(src, path, site)

# Wrap for use with custom interface
wrapper = FourierParametersWrapper(fourier_params)

# Now can use compute_fas
fas = compute_fas(wrapper, freq, mag, dist)
```
"""
struct FourierParametersWrapper <: AbstractFASModel
    params::FourierParameters
end

function compute_fas(model::FourierParametersWrapper, freq, mag, dist)
    # Delegate to the existing fourier_spectral_ordinate function
    return fourier_spectral_ordinate(freq, mag, dist, model.params)
end


# simple container type (internal)
struct ModelPair
    fas_model::AbstractFASModel
    duration_model::AbstractDurationModel
end

# 2-arg constructor for tests and convenience
FourierParametersWrapper(
    fas_model::AbstractFASModel,
    duration_model::AbstractDurationModel,
) = ModelPair(fas_model, duration_model)



"""
    ExistingDurationWrapper <: AbstractDurationModel

Wrapper to make existing RandomVibrationParameters compatible with AbstractDurationModel interface.

Since `excitation_duration` requires both FourierParameters and RandomVibrationParameters,
this wrapper stores both to enable proper delegation to the existing duration calculation.

# Example
```julia
# Create traditional parameters
fourier_params = FourierParameters(...)
rvt_params = RandomVibrationParameters(:BT15)

# Wrap for use with custom interface
# Note: Must provide both parameters
wrapper = ExistingDurationWrapper(rvt_params, fourier_params)

# Now can use compute_duration
duration = compute_duration(wrapper, mag, dist)
```

# Important
This wrapper is only used when BOTH FourierParameters and RandomVibrationParameters
are provided (the fully concrete case). If you're using a custom FAS model, you must
also use a custom duration model.
"""
struct ExistingDurationWrapper <: AbstractDurationModel
    rvt_params::RandomVibrationParameters
    fourier_params::FourierParameters
end

function compute_duration(model::ExistingDurationWrapper, mag, dist)
    # Delegate to the existing excitation_duration function
    # This requires both FourierParameters and RandomVibrationParameters
    return excitation_duration(mag, dist, model.fourier_params, model.rvt_params)
end



# =====================================
# Utilities for Model Development
# =====================================

"""
    validate_fas_model(model::AbstractFASModel; freq_range=(0.01, 100), mag=6.0, dist=10.0)

Validate that a FAS model is properly implemented and ForwardDiff-compatible.
"""
function validate_fas_model(model::AbstractFASModel;
    freq_range=(0.01, 100),
    mag=6.0,
    dist=10.0)
    println("Validating FAS model...")

    # Test basic computation
    try
        fas = compute_fas(model, 1.0, mag, dist)
        println("✓ Basic computation: FAS = $fas")
    catch e
        println("✗ Basic computation failed: $e")
        return false
    end

    # Test ForwardDiff compatibility
    try
        # Test gradient with respect to frequency
        grad_f = ForwardDiff.derivative(f -> compute_fas(model, f, mag, dist), 1.0)
        println("✓ ForwardDiff w.r.t. frequency: ∂FAS/∂f = $grad_f")

        # Test gradient with respect to magnitude
        grad_m = ForwardDiff.derivative(m -> compute_fas(model, 1.0, m, dist), mag)
        println("✓ ForwardDiff w.r.t. magnitude: ∂FAS/∂m = $grad_m")

        # Test gradient with respect to distance
        grad_r = ForwardDiff.derivative(r -> compute_fas(model, 1.0, mag, r), dist)
        println("✓ ForwardDiff w.r.t. distance: ∂FAS/∂r = $grad_r")
    catch e
        println("✗ ForwardDiff compatibility failed: $e")
        return false
    end

    # Test numerical stability
    freqs = exp10.(range(log10(freq_range[1]), log10(freq_range[2]), length=50))
    fas_values = [compute_fas(model, f, mag, dist) for f in freqs]

    if any(isnan.(fas_values)) || any(isinf.(fas_values))
        println("⚠ Warning: NaN or Inf values detected in FAS computation")
    else
        println("✓ Numerical stability check passed")
    end

    println("Validation complete!")
    return true
end

"""
    validate_duration_model(model::AbstractDurationModel; mag=6.0, dist=10.0)

Validate that a duration model is properly implemented.
"""
function validate_duration_model(model::AbstractDurationModel; mag=6.0, dist=10.0)
    println("Validating duration model...")

    # Test basic computation
    try
        dur = compute_duration(model, mag, dist)
        println("✓ Basic computation: Duration = $dur seconds")

        if dur <= 0
            println("⚠ Warning: Duration is non-positive")
        end
    catch e
        println("✗ Basic computation failed: $e")
        return false
    end

    # Test ForwardDiff compatibility
    try
        grad_m = ForwardDiff.derivative(m -> compute_duration(model, m, dist), mag)
        println("✓ ForwardDiff w.r.t. magnitude: ∂Duration/∂m = $grad_m")

        grad_r = ForwardDiff.derivative(r -> compute_duration(model, mag, r), dist)
        println("✓ ForwardDiff w.r.t. distance: ∂Duration/∂r = $grad_r")
    catch e
        println("✗ ForwardDiff compatibility failed: $e")
        return false
    end

    println("Validation complete!")
    return true
end

# =====================================
# Export public interface
# =====================================

export compute_fas, compute_duration
export FunctionalFASModel, FunctionalDurationModel
export CustomFASModel, CustomDurationModel
export HybridFASModel
export FourierParametersWrapper, ExistingDurationWrapper
export rvt_response_spectral_ordinate_custom
export validate_fas_model, validate_duration_model

end
