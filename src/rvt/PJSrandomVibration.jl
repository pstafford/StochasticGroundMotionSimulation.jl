const xn31i, wn31i = gausslegendre(31)


@doc raw"""
	zeros_extrema_frequencies(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real}

Defines the frequencies of extrema and zero-crossings using moments ``m_0``, ``m_2`` and ``m_4``. Returns a tuple of ``(f_z,f_e)``.

The frequency of zero crossings is:
```math
f_{z} = \frac{1}{2\pi}\sqrt{\frac{m_2}{m_0}}
```
The frequency of extrema is:
```math
f_{e} = \frac{1}{2\pi}\sqrt{\frac{m_4}{m_2}}
```

See also: [`zeros_extrema_numbers`](@ref)
"""
function zeros_extrema_frequencies(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real}
    mi = spectral_moments([0, 2, 4], m, r_ps, fas, sdof, glxi=glxi, glwi=glwi)
    m_0 = mi.m0
    m_2 = mi.m2
    m_4 = mi.m4
    fzero = sqrt(m_2 / m_0) / (2π)
    fextrema = sqrt(m_4 / m_2) / (2π)
    return fzero, fextrema
end


@doc raw"""
	zeros_extrema_numbers(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Real,T<:Real}

Defines the numbers of extrema and zero-crossings using moments ``m_0``, ``m_2`` and ``m_4``. Returns a tuple of ``(2f_z D_{ex}, 2f_e D_{ex})``.

The frequency of zero crossings is:
```math
n_{z} = \frac{D_{ex}}{\pi}\sqrt{\frac{m_2}{m_0}}
```
The frequency of extrema is:
```math
n_{e} = \frac{D_{ex}}{\pi}\sqrt{\frac{m_4}{m_2}}
```
In both cases, the excitaion duration, ``D_{ex}`` is computed from the `excitation_duration` function.

See also: [`zeros_extrema_frequencies`](@ref), [`excitation_duration`](@ref)
"""
function zeros_extrema_numbers(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real}
    fz, fe = zeros_extrema_frequencies(m, r_ps, fas, sdof, glxi=glxi, glwi=glwi)
    Dex = excitation_duration(m, r_ps, fas, rvt)
    return 2fz * Dex, 2fe * Dex
end



@doc raw"""
	rvt_response_spectral_ordinate(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real}

Response spectral ordinate (units of ``g``) for the specified scenario.

The spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectral_ordinate(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real}
    # get the duration metrics
    Drms, Dex, ~ = rms_duration(m, r_ps, fas, sdof, rvt)
    # compute the necessary spectral moments
    if rvt.pf_method == :DK80
        order = [0, 1, 2]
    elseif rvt.pf_method == :CL56
        order = [0, 2, 4]
    end
    mi = spectral_moments(order, m, r_ps, fas, sdof, glxi=glxi, glwi=glwi)

    # get the rms response
    y_rms = sqrt(mi.m0 / Drms)
    # get the peak factor
    pf = peak_factor(Dex, mi, rvt, glxi=glxi, glwi=glwi)
    # response spectral value (in g)
    Sa = pf * y_rms / 9.80665
    return Sa
end


# @doc raw"""
# 	rvt_response_spectral_ordinate(period::U, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real,U<:Float64}

# Response spectral ordinate (units of ``g``) for the specified scenario.

# The spectral ordinate is computed using the expression:
# ```math
# S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
# ```
# where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`.

# See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum`](@ref), [`rvt_response_spectrum!`](@ref)
# """
# function rvt_response_spectral_ordinate(period::U, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i, damping=0.05) where {S<:Real,T<:Real,U<:Float64}
#     # create a sdof instance
#     sdof = Oscillator(1.0 / period, damping)
#     return rvt_response_spectral_ordinate(m, r_ps, fas, sdof, rvt, glxi=glxi, glwi=glwi)
# end



@doc raw"""
	rvt_response_spectrum(period::Vector{U}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real,U<:Float64}

Response spectrum (units of ``g``) for the vector of periods `period` and the specified scenario.

Each spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`. The various terms are all functions of the oscillator period.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectrum(period::Vector{U}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real,U<:Float64}
    if S <: Dual
        Sa = zeros(S, length(period))
    else
        V = get_parametric_type(fas)
        Sa = zeros(V, length(period))
    end
    for i in 1:length(period)
        @inbounds Sa[i] = rvt_response_spectral_ordinate(period[i], m, r_ps, fas, rvt, glxi=glxi, glwi=glwi)
    end
    return Sa
end


@doc raw"""
	rvt_response_spectrum!(Sa::Vector{U}, period::Vector{V}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters) where {S<:Real,T<:Real,U<:Real,V<:Float64}

In-place response spectrum (units of ``g``) for the vector of periods `period` and the specified scenario.

Each spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`. The various terms are all functions of the oscillator period.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectrum!(Sa::Vector{U}, period::Vector{V}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real,U<:Real,V<:Float64}
    for i in 1:length(period)
        @inbounds Sa[i] = rvt_response_spectral_ordinate(period[i], m, r_ps, fas, rvt, glxi=glxi, glwi=glwi)
    end
end


"""
    rvt_response_spectral_ordinate(T, mag, dist, fas_model, duration_model; kwargs...)

Unified interface for computing response spectral ordinates.

Can accept either:
- Traditional FourierParameters and RandomVibrationParameters
- Custom AbstractFASModel and AbstractDurationModel instances

# Supported Combinations
1. FourierParameters + RandomVibrationParameters (both concrete)
2. FourierParameters + AbstractDurationModel (concrete FAS, custom duration)
3. AbstractFASModel + AbstractDurationModel (both custom)

# Unsupported Combination
- AbstractFASModel + RandomVibrationParameters is NOT supported because
  RandomVibrationParameters.excitation_duration requires FourierParameters.
  Use a custom duration model instead.

# Implementation
Three method definitions ensure proper dispatch hierarchy and avoid
ambiguity with the earlier method at line 104.
"""

# 1. Both concrete types - this overrides the method at line 104
function rvt_response_spectral_ordinate(
    T::U,
    mag::S,
    dist::V,
    fas_model::FourierParameters,
    duration_model::RandomVibrationParameters;
    damping::Real=0.05,
    glxi::Vector{Float64}=xn31i,
    glwi::Vector{Float64}=wn31i,
    kwargs...
) where {S<:Real,V<:Real,U<:Float64}
    # Both are concrete - pass both to wrapper for proper duration calculation
    return rvt_response_spectral_ordinate_custom(
        T, mag, dist,
        FourierParametersWrapper(fas_model),
        ExistingDurationWrapper(duration_model, fas_model),  # Pass both!
        damping, glxi, glwi
    )
end

# 2. Concrete FAS with custom duration
function rvt_response_spectral_ordinate(
    T::Real,
    mag::Real,
    dist::Real,
    fas_model::FourierParameters,
    duration_model::AbstractDurationModel;
    damping::Real=0.05,
    glxi::Vector{Float64}=xn31i,
    glwi::Vector{Float64}=wn31i,
    kwargs...
)
    return rvt_response_spectral_ordinate_custom(
        T, mag, dist,
        FourierParametersWrapper(fas_model),
        duration_model,  # Already abstract, no wrapping needed
        damping, glxi, glwi
    )
end

# 3. Custom FAS with concrete duration - NOT SUPPORTED
function rvt_response_spectral_ordinate(
    T::Real,
    mag::Real,
    dist::Real,
    fas_model::AbstractFASModel,
    duration_model::RandomVibrationParameters;
    kwargs...
)
    error("""
    Mixing custom FAS models (AbstractFASModel) with concrete duration parameters 
    (RandomVibrationParameters) is not supported.

    Reason: RandomVibrationParameters.excitation_duration requires FourierParameters,
    which are not available when using a custom FAS model.

    Solutions:
    1. Use a custom duration model instead: FunctionalDurationModel((m,r,p) -> ...)
    2. Use concrete FAS parameters: FourierParameters(...)
    3. Use the old interface if you need both concrete types

    Example with custom duration:
        dur_model = FunctionalDurationModel((m, r, p) -> p[1] + p[2]*m + p[3]*log10(r), [5.0, 0.5, 1.0])
        Sa = rvt_response_spectral_ordinate(T, mag, dist, fas_model, dur_model)
    """)
end

# 4. Both custom types
function rvt_response_spectral_ordinate(
    T::Real,
    mag::Real,
    dist::Real,
    fas_model::AbstractFASModel,
    duration_model::AbstractDurationModel;
    damping::Real=0.05,
    glxi::Vector{Float64}=xn31i,
    glwi::Vector{Float64}=wn31i,
    kwargs...
)
    return rvt_response_spectral_ordinate_custom(
        T, mag, dist,
        fas_model,
        duration_model,
        damping, glxi, glwi
    )
end
