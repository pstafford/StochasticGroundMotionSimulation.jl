


# function to compute the peak over rms ratio, a.k.a. peak factor
function peak_factor_cl56(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=35) where {S<:Float64,T<:Real}
	# get the numbers of zero crossing and extrema
	rvt = RandomVibrationParameters(:CL56)
	n_z, n_e = zeros_extrema_numbers(m, r_ps, fas, sdof, rvt)
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e

	# compute the Gauss Legendre nodes and weights
	xi, wi = gausslegendre(nodes)

	z_min = 0.0
	z_max = 8.0
	dfac = (z_max-z_min)/2
	pfac = (z_max+z_min)/2

	# scale the nodes to the appropriate range
	zi = @. dfac * xi + pfac

	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e

	int_sum = dfac * dot( wi, integrand.(zi) )

	pf = sqrt(2.0) * int_sum
	return pf
end


function peak_factor_cl56(m::S, r_ps::T, Dex::U, m0::V, fas::FourierParameters, sdof::Oscillator; nodes::Int=35) where {S<:Float64,T<:Real,U<:Real,V<:Real}
	# get all necessary spectral moments
	mi = spectral_moments([2, 4], m, r_ps, fas, sdof)
	m2 = mi[1]
	m4 = mi[2]
	# get the peak factor
	n_z = Dex * sqrt( m2 / m0 ) / π
	n_e = Dex * sqrt( m4 / m2 ) / π
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e

	# compute the Gauss Legendre nodes and weights
	xi, wi = gausslegendre(nodes)

	z_min = 0.0
	z_max = 8.0
	dfac = (z_max-z_min)/2
	pfac = (z_max+z_min)/2

	# scale the nodes to the appropriate range
	zi = @. dfac * xi + pfac

	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e

	int_sum = dfac * dot( wi, integrand.(zi) )

	pf = sqrt(2.0) * int_sum
	return pf
end



function peak_factor_cl56(n_z::T, n_e::T; nodes::Int=35) where T<:Real
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e

	# compute the Gauss Legendre nodes and weights
	xi, wi = gausslegendre(nodes)

	z_min = 0.0
	z_max = 8.0
	dfac = (z_max-z_min)/2
	pfac = (z_max+z_min)/2

	# scale the nodes to the appropriate range
	zi = @. dfac * xi + pfac

	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e

	int_sum = dfac * dot( wi, integrand.(zi) )

	pf = sqrt(2.0) * int_sum
	return pf
end


function peak_factor_cl56_gk(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Float64,T<:Real}
	# get the numbers of zero crossing and extrema
	rvt = RandomVibrationParameters(:CL56)
	n_z, n_e = zeros_extrema_numbers(m, r_ps, fas, sdof, rvt)
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e

	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e

	return sqrt(2.0) * quadgk(integrand, 0.0, Inf)[1]
end


function peak_factor_integrand_cl56(z, n_z, n_e)
  ξ = n_z / n_e
  return 1.0 - (1.0 - ξ*exp(-z^2))^n_e
end



function vanmarcke_cdf(x, n_z::T, δeff::T) where T<:Real
	if x < eps()
		return zero(T)
	else
		xsq = x^2
		eδe = exp( -sqrt(π/2) * δeff * x )
		Fx = ( 1.0 - exp(-xsq/2) ) * exp( -n_z * (1.0 - eδe)/(exp(xsq/2) - 1.0) )
		return Fx
	end
end

function vanmarcke_ccdf(x, n_z::T, δeff::T) where T<:Real
	return 1.0 - vanmarcke_cdf(x, n_z, δeff)
end

function peak_factor_dk80(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=30) where {S<:Float64,T<:Real}
	# compute first three spectral moments
	mi = spectral_moments([0, 1, 2], m, r_ps, fas, sdof)
	m0 = mi[1]
	m1 = mi[2]
	m2 = mi[3]

	# bandwidth, and effective bandwidth
	δ = sqrt( 1.0 - (m1^2 / (m0*m2)) )
	δeff = δ^1.2

	# excitation duration
	Dex = boore_thompson_2014(m, r_ps, fas)
	# number of zero crossings
	n_z = Dex * sqrt( m2/m0 ) / π

	# compute the Gauss Legendre nodes and weights
	xi, wi = gausslegendre(nodes)

	z_min = 0.0
	z_max = 6.0
	dfac = (z_max-z_min)/2
	pfac = (z_max+z_min)/2

	# scale the nodes to the appropriate range
	zi = @. dfac * xi + pfac
	# evaluate integrand
	yy = vanmarcke_ccdf.(zi, n_z, δeff)

	return int_sum = dfac * dot( wi, yy )
end


function peak_factor_dk80(m::S, r_ps::T, Dex::U, m0::V, fas::FourierParameters, sdof::Oscillator; nodes::Int=30) where {S<:Float64,T<:Real,U<:Real,V<:Real}
	# compute first three spectral moments
	mi = spectral_moments([1, 2], m, r_ps, fas, sdof)
	m1 = mi[1]
	m2 = mi[2]

	# bandwidth, and effective bandwidth
	δ = sqrt( 1.0 - (m1^2 / (m0*m2)) )
	δeff = δ^1.2

	# number of zero crossings
	n_z = Dex * sqrt( m2/m0 ) / π

	# compute the Gauss Legendre nodes and weights
	xi, wi = gausslegendre(nodes)

	z_min = 0.0
	z_max = 6.0
	dfac = (z_max-z_min)/2
	pfac = (z_max+z_min)/2

	# scale the nodes to the appropriate range
	zi = @. dfac * xi + pfac
	# evaluate integrand
	yy = vanmarcke_ccdf.(zi, n_z, δeff)

	return int_sum = dfac * dot( wi, yy )
end


function peak_factor_dk80_gk(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Float64,T<:Real}
	# compute first three spectral moments
	mi = spectral_moments([0, 1, 2], m, r_ps, fas, sdof)
	m0 = mi[1]
	m1 = mi[2]
	m2 = mi[3]

	# bandwidth, and effective bandwidth
	δ = sqrt( 1.0 - (m1^2 / (m0*m2)) )
	δeff = δ^1.2

	# excitation duration
	Dex = boore_thompson_2014(m, r_ps, fas)
	# number of zero crossings
	n_z = Dex * sqrt( m2/m0 ) / π

	integrand(x) = vanmarcke_ccdf(x, n_z, δeff)

	return quadgk(integrand, 0.0, Inf)[1]
end





"""
	peak_factor(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; pf_method::Symbol=:DK80) where {S<:Float64,T<:Real}

Peak factor u_max / u_rms with a switch of `pf_method` to determine the approach adopted.
`pf_method` can currently be one of:
	- `:CL56` for Cartright Longuet-Higgins (1956)
	- `:DK80` for Der Kiureghian (1980), building on Vanmarcke (1975)

Defaults to `:DK80`.
"""
function peak_factor(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real}
	if rvt.pf_method == :CL56
		return peak_factor_cl56(m, r_ps, fas, sdof)
	elseif rvt.pf_method == :DK80
		return peak_factor_dk80(m, r_ps, fas, sdof)
	else
		U = get_parametric_type(fas)
		return U(NaN)
	end
end


"""
	peak_factor(m::T, r::T, fas::FourierParameters, sdof::Oscillator; pf_method::Symbol=:DK80) where T<:Float64

Peak factor u_max / u_rms with a switch of `pf_method` to determine the approach adopted.
`pf_method` can currently be one of:
	- `:CL56` for Cartright Longuet-Higgins (1956)
	- `:DK80` for Der Kiureghian (1980), building on Vanmarcke (1975)

Defaults to `:DK80`.
"""
function peak_factor(m::S, r_ps::T, Dex::U, m0::V, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real,U<:Real,V<:Real}
	if rvt.pf_method == :CL56
		return peak_factor_cl56(m, r_ps, Dex, m0, fas, sdof)
	elseif rvt.pf_method == :DK80
		return peak_factor_dk80(m, r_ps, Dex, m0, fas, sdof)
	else
		if U <: Float64
			return V(NaN)
		else
			return U(NaN)
		end
	end
end
