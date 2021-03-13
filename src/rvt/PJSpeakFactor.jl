


# function to compute the peak over rms ratio, a.k.a. peak factor
function peak_factor_cl56(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; intervals::Int=101) where T<:Float64
	# get the numbers of zero crossing and extrema
	n_z, n_e = zeros_extrema_numbers(m, r, fas, sdof)
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e
	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e
	# compute the peak factor
	xx = collect(range(0.0, stop=10.0, length=intervals))
	yy = integrand.(xx)
	int_sum = trapezoidal_rule(xx, yy)
	# int_sum = simpsons_rule(xx, yy)

	pf = sqrt(2.0) * int_sum
	return pf
end

function peak_factor_cl56(n_z::T, n_e::T; intervals::Int=101) where T<:Float64
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e
	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e
	# compute the peak factor
	xx = collect(range(0.0, stop=10.0, length=intervals))
	yy = integrand.(xx)
	# int_sum = simpsons_rule(xx, yy)
	int_sum = trapezoidal_rule(xx, yy)

	pf = sqrt(2.0) * int_sum
	return pf
end


function peak_factor_integrand_cl56(z, n_z, n_e)
  ξ = n_z / n_e
  return 1.0 - (1.0 - ξ*exp(-z^2))^n_e
end

function vanmarcke_cdf(x, n_z, δeff)
	if x < eps()
		return 0.0
	else
		xsq = x^2
		eδe = exp( -sqrt(π/2) * δeff * x )
		Fx = ( 1.0 - exp(-xsq/2) ) * exp( -n_z * (1.0 - eδe)/(exp(xsq/2) - 1.0) )
		return Fx
	end
end

function vanmarcke_ccdf(x, n_z, δeff)
	return 1.0 - vanmarcke_cdf(x, n_z, δeff)
end

function peak_factor_dk80(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; intervals::Int=101) where T<:Float64
	# compute first three spectral moments
	mi = spectral_moments([0, 1, 2], m, r, fas, sdof)
	m0 = mi[1]
	m1 = mi[2]
	m2 = mi[3]

	# bandwidth, and effective bandwidth
	δ = sqrt( 1.0 - (m1^2 / (m0*m2)) )
	δeff = δ^1.2

	# excitation duration
	Dex = boore_thompson_2014(m, r, fas)
	# number of zero crossings
	n_z = Dex * sqrt( m2/m0 ) / π

	xx = collect(range(0.0, stop=10.0, length=intervals))
	yy = vanmarcke_ccdf.(xx, n_z, δeff)
	# return simpsons_rule(xx, yy)

	return trapezoidal_rule(xx, yy)
end


function peak_factor_dk80_gk(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}) where T<:Float64
	# compute first three spectral moments
	mi = spectral_moments([0, 1, 2], m, r, fas, sdof)
	m0 = mi[1]
	m1 = mi[2]
	m2 = mi[3]

	# bandwidth, and effective bandwidth
	δ = sqrt( 1.0 - (m1^2 / (m0*m2)) )
	δeff = δ^1.2

	# excitation duration
	Dex = boore_thompson_2014(m, r, fas)
	# number of zero crossings
	n_z = Dex * sqrt( m2/m0 ) / π

	integrand(x) = vanmarcke_ccdf(x, n_z, δeff)

	return quadgk(integrand, 0.0, Inf)[1]
end





"""
	peak_factor(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; pf_method::Symbol=:DK80) where T<:Float64

Peak factor u_max / u_rms with a switch of `pf_method` to determine the approach adopted.
`pf_method` can currently be one of:
	- `:CL56` for Cartright Longuet-Higgins (1956)
	- `:DK80` for Der Kiureghian (1980), building on Vanmarcke (1975)

Defaults to `:DK80`.
"""
function peak_factor(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; pf_method::Symbol=:DK80) where T<:Float64
	if pf_method == :CL56
		return peak_factor_cl56(m, r, fas, sdof)
	else
		return peak_factor_dk80(m, r, fas, sdof)
	end
end
