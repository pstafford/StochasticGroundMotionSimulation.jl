
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
	fzero = sqrt( m_2 / m_0 ) / (2π)
	fextrema = sqrt( m_4 / m_2 ) / (2π)
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
    return 2fz*Dex, 2fe*Dex
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


@doc raw"""
	rvt_response_spectral_ordinate(period::U, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real,U<:Float64}

Response spectral ordinate (units of ``g``) for the specified scenario.

The spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectral_ordinate(period::U, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters; glxi::Vector{Float64}=xn31i, glwi::Vector{Float64}=wn31i) where {S<:Real,T<:Real,U<:Float64}
  	# create a sdof instance
  	sdof = Oscillator(1.0/period)
	return rvt_response_spectral_ordinate(m, r_ps, fas, sdof, rvt, glxi=glxi, glwi=glwi)
end



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
		Sa = zeros(S,length(period))
	else
		V = get_parametric_type(fas)
		Sa = zeros(V,length(period))
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
