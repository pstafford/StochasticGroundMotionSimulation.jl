"""
	@kwdef struct SpectralMoments{T<:Real}

Struct holding spectral moments `m0` to `m4`
"""
@kwdef struct SpectralMoments{T<:Real}
    m0::T = NaN
    m1::T = NaN
    m2::T = NaN
    m3::T = NaN
    m4::T = NaN
end


"""
	create_spectral_moments(order::Vector{Int}, value::Vector{T}) where {T<:Real}

Create a `SpectralMoments` instance from vectors of order integers and moment values.
Allows for encapsulation of named moments within the returned instance.
"""
function create_spectral_moments(order::Vector{Int}, value::Vector{T}) where {T<:Real}
	dict = Dict{Symbol, T}()
	for (i, o) in enumerate(order)
		dict[Symbol("m$o")] = value[i]
	end
	return SpectralMoments{T}(; dict...)
end


@doc raw"""
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

	for i in 2:lastindex(lnflims)
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
    # return 2.0 * int_sum
	return create_spectral_moments([order], [2.0 * int_sum])
end


@doc raw"""
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

	for i in 2:lastindex(lnflims)
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
    # return 2.0 * int_sumi
    return create_spectral_moments(order, 2.0 * int_sumi)
end



@doc raw"""
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
    # return 2.0 * int_sumi
    return create_spectral_moments(order, 2.0 * int_sumi)
end