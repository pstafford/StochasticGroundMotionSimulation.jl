
"""
	spectral_moment(order::Int, m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=31, control_freqs::Vector{Float64}=[1e-3, 1e-1, 1.0, 10.0, 100.0, 300.0] ) where {S<:Real,T<:Real}

Compute spectral moment of a specified order.

Evaluates the expression:
```math
	m_k = 2\int_{0}^{\infty} \left(2\pi f\right)^k |H(f;f_n,\zeta_n)|^2 |A(f)|^2 df
```
where ``k`` is the order of the moment.

Integration is performed using Gauss-Legendre integration using `nodes` nodes and weights.
The integration domain is partitioned over the `control_freqs` as well as two inserted frequencies at `f_n/1.5` and `f_n*1.5` in order to ensure good approximation of the integral around the `sdof` resonant frequency.

See also: [`spectral_moments`](@ref)
"""
function spectral_moment(order::Int, m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=31, control_freqs::Vector{Float64}=[1e-3, 1e-1, 1.0, 10.0, 100.0, 300.0] ) where {S<:Real,T<:Real}
	# pre-allocate for squared fourier amplitude spectrum
	if S <: Dual
		Af2 = Vector{S}(undef, nodes)
		Hf2 = Vector{S}(undef, nodes)
		int_sum = zero(S)
	else
		U = get_parametric_type(fas)
		Af2 = Vector{U}(undef, nodes)
		Hf2 = Vector{U}(undef, nodes)
		int_sum = zero(U)
	end
	# partition the integration domain to make sure the integral captures the key change in the integrand
	f_n = sdof.f_n
	# note that we start slightly above zero to avoid a numerical issue with the frequency dependent Q(f) gradients
	fidLO = findall(control_freqs .< f_n/1.5)
	fidHI = findall(control_freqs .> f_n*1.5)
	# perform the integration with a logarithmic transformation
	lnflims = log.([ control_freqs[fidLO]; f_n/1.5; f_n*1.5; control_freqs[fidHI] ])

	# compute the Gauss Legendre nodes and weights
	xi, wi = gausslegendre(nodes)

	for i in 2:length(lnflims)
		@inbounds dfac = (lnflims[i]-lnflims[i-1])/2
		@inbounds pfac = (lnflims[i]+lnflims[i-1])/2
		lnfi = @. dfac * xi + pfac
		fi = exp.(lnfi)
		squared_transfer!(Hf2, fi, sdof)
		squared_fourier_spectrum!(Af2, fi, m, r_ps, fas)
		# note that the integrand here is scaled by exp(lnfi)=fi for the logarithmic transformation of the integrand
		Yf2 = @. (2π * fi)^order * Hf2 * Af2 * fi
		# weighted combination of amplitudes with Gauss-Legendre weights
		int_sum += dfac * dot( wi, Yf2 )
	end
    return 2.0 * int_sum
end


"""
	spectral_moments(order::Vector{Int}, m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=31, control_freqs::Vector{Float64}=[1e-3, 1e-1, 1.0, 10.0, 100.0, 300.0] ) where {S<:Real,T<:Real}

Compute a vector of spectral moments for the specified `order`.

Evaluates the expression:
```math
	m_k = 2\int_{0}^{\infty} \left(2\pi f\right)^k |H(f;f_n,\zeta_n)|^2 |A(f)|^2 df
```
for each order, where ``k`` is the order of the moment.

Integration is performed using Gauss-Legendre integration using `nodes` nodes and weights.
The integration domain is partitioned over the `control_freqs` as well as two inserted frequencies at `f_n/1.5` and `f_n*1.5` in order to ensure good approximation of the integral around the `sdof` resonant frequency.

See also: [`spectral_moment`](@ref), [`spectral_moments_gk`](@ref)
"""
function spectral_moments(order::Vector{Int}, m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=31, control_freqs::Vector{Float64}=[1e-3, 1e-1, 1.0, 10.0, 100.0, 300.0] ) where {S<:Real,T<:Real}
	# pre-allocate for squared fourier amplitude spectrum
	if S <: Dual
		Af2 = Vector{S}(undef, nodes)
		Hf2 = Vector{S}(undef, nodes)
		int_sumi = zeros(S, length(order))
	else
		U = get_parametric_type(fas)
		Af2 = Vector{U}(undef, nodes)
		Hf2 = Vector{U}(undef, nodes)
		int_sumi = zeros(U, length(order))
	end
	# partition the integration domain to make sure the integral captures the key change in the integrand
	f_n = sdof.f_n
 	# note that we start slightly above zero to avoid a numerical issue with the frequency dependent Q(f) gradients
	fidLO = findall(control_freqs .< f_n/1.5)
	fidHI = findall(control_freqs .> f_n*1.5)
	# perform the integration with a logarithmic transformation
	lnflims = log.([ control_freqs[fidLO]; f_n/1.5; f_n*1.5; control_freqs[fidHI] ])

	# compute the Gauss Legendre nodes and weights
	xi, wi = gausslegendre(nodes)

	# make sure the orders are listed as increasing for the following loop approach
	sort!(order)
	dorder = diff(order)

	for i in 2:length(lnflims)
		@inbounds dfac = (lnflims[i]-lnflims[i-1])/2
		@inbounds pfac = (lnflims[i]+lnflims[i-1])/2
		lnfi = @. dfac * xi + pfac
		fi = exp.(lnfi)
		squared_transfer!(Hf2, fi, sdof)
		squared_fourier_spectrum!(Af2, fi, m, r_ps, fas)
		# compute default zeroth order integrand amplitude
		# note that the integrand here is scaled by exp(lnfi)=fi for the logarithmic transformation of the integrand
		Yf2 = @. Hf2 * Af2 * fi

		for (idx, o) in enumerate(order)
			if idx == 1
				Yf2 .*= (2π * fi).^o
			else
				@inbounds Yf2 .*= (2π * fi).^(dorder[idx-1])
			end
			# weighted combination of amplitudes with Gauss-Legendre weights
			@inbounds int_sumi[idx] += dfac * dot( wi, Yf2 )
		end
	end
    return 2.0 * int_sumi
end



"""
	spectral_moments_gk(order::Vector{Int}, m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Real,T<:Real}

Compute a vector of spectral moments for the specified `order` using Gauss-Kronrod integration from the `QuadGK.jl` package.

Evaluates the expression:
```math
	m_k = 2\int_{0}^{\infty} \left(2\pi f\right)^k |H(f;f_n,\zeta_n)|^2 |A(f)|^2 df
```
for each order, where ``k`` is the order of the moment.

Integration is performed using adaptive Gauss-Kronrod integration with the domain split over two intervals from ``[0,f_n]`` and ``[f_n,\infty]`` to ensure that the resonant peak is not missed.

Note that due to the default tolerances, the moments computed by this method are more accurate than those from `spectral_moments` using the Gauss-Legendre approximation. However, this method is also significantly slower, and cannot be used within an automatic differentiation environment.

See also: [`spectral_moment`](@ref), [`spectral_moments`](@ref)
"""
function spectral_moments_gk(order::Vector{Int}, m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Real,T<:Real}
	int_sumi = zeros(length(order))
	for (i, o) in enumerate(order)
		moment_integrand(f) = squared_transfer(f, sdof) * squared_fourier_spectral_ordinate(f, m, r_ps, fas) * (2π*f)^o
		int_sumi[i] = quadgk(moment_integrand, 0.0, sdof.f_n, Inf)[1]
	end
    return 2.0 * int_sumi
end


"""
	zeros_extrema_frequencies(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Real,T<:Real}

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
function zeros_extrema_frequencies(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Real,T<:Real}
	mi = spectral_moments([0, 2, 4], m, r_ps, fas, sdof)
	m_0 = mi[1]
	m_2 = mi[2]
	m_4 = mi[3]
	fzero = sqrt( m_2 / m_0 ) / (2π)
	fextrema = sqrt( m_4 / m_2 ) / (2π)
	return fzero, fextrema
end


"""
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
function zeros_extrema_numbers(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Real,T<:Real}
    fz, fe = zeros_extrema_frequencies(m, r_ps, fas, sdof)
    Dex = excitation_duration(m, r_ps, fas, rvt)
    return 2fz*Dex, 2fe*Dex
end



"""
	rvt_response_spectral_ordinate(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Real,T<:Real}

Response spectral ordinate (units of ``g``) for the specified scenario.

The spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectral_ordinate(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Real,T<:Real}
	# get the duration metrics
	Drms, Dex, Dratio = rms_duration(m, r_ps, fas, sdof, rvt)
	# get the rms response
	m_0 = spectral_moment(0, m, r_ps, fas, sdof)
	y_rms = sqrt( m_0 / Drms )
	# get the peak factor
	pf = peak_factor(m, r_ps, Dex, m_0, fas, sdof, rvt)
	# response spectral value (in g)
	Sa = pf * y_rms / 9.80665
	return Sa
end


"""
	rvt_response_spectral_ordinate(period::U, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters) where {S<:Real,T<:Real,U<:Float64}

Response spectral ordinate (units of ``g``) for the specified scenario.

The spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectral_ordinate(period::U, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters) where {S<:Real,T<:Real,U<:Float64}
  	# create a sdof instance
  	sdof = Oscillator(1.0/period)
	return rvt_response_spectral_ordinate(m, r_ps, fas, sdof, rvt)
end



"""
	rvt_response_spectrum(period::Vector{U}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters) where {S<:Real,T<:Real,U<:Float64}

Response spectrum (units of ``g``) for the vector of periods `period` and the specified scenario.

Each spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`. The various terms are all functions of the oscillator period.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectrum(period::Vector{U}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters) where {S<:Real,T<:Real,U<:Float64}
	if S <: Dual
		Sa = zeros(S,length(period))
	else
		V = get_parametric_type(fas)
		Sa = zeros(V,length(period))
	end
	for i in 1:length(period)
		@inbounds Sa[i] = rvt_response_spectral_ordinate(period[i], m, r_ps, fas, rvt)
  	end
  	return Sa
end


"""
	rvt_response_spectrum!(Sa::Vector{U}, period::Vector{V}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters) where {S<:Real,T<:Real,U<:Real,V<:Float64}

In-place response spectrum (units of ``g``) for the vector of periods `period` and the specified scenario.

Each spectral ordinate is computed using the expression:
```math
S_a = \psi \sqrt{\frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment from `spectral_moment`, and ``D_{rms}`` is the RMS duration computed from `rms_duration`. The various terms are all functions of the oscillator period.

See also: [`rvt_response_spectral_ordinate`](@ref), [`rvt_response_spectrum!`](@ref)
"""
function rvt_response_spectrum!(Sa::Vector{U}, period::Vector{V}, m::S, r_ps::T, fas::FourierParameters, rvt::RandomVibrationParameters) where {S<:Real,T<:Real,U<:Real,V<:Float64}
	for i in 1:length(period)
		@inbounds Sa[i] = rvt_response_spectral_ordinate(period[i], m, r_ps, fas, rvt)
  	end
end
