
@doc raw"""
	peak_factor_cl56(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=35) where {S<:Real,T<:Real}

Peak factor computed using the Cartwright and Longuet-Higgins (1956) formulation.

Evaluates the integral to obtain the peak factor ``\psi``:
```math
\psi = \sqrt{2} \int_0^\infty 1 - \left( 1 - \xi \exp\left( -z^2 \right)\right)^{n_e} dz
```
where ``n_e`` is the number of extrema, ``\xi`` is the ratio ``n_z/n_e`` with ``n_z`` being the number of zero crossings.

The integral is evaluated using Gauss-Legendre integration -- and is suitable for use within an automatic differentiation environment.

See also: [`peak_factor`](@ref), [`peak_factor_dk80`](@ref)
"""
function peak_factor_cl56(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=35) where {S<:Real,T<:Real}
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

	return sqrt(2.0) * int_sum
end


@doc raw"""
	peak_factor_cl56(Dex::U, mi::SpectralMoments; nodes::Int=35) where {U<:Real}

Peak factor computed using the Cartwright and Longuet-Higgins (1956) formulation, using pre-computed `Dex` and `m0` values.

Evaluates the integral to obtain the peak factor ``\psi``:
```math
\psi = \sqrt{2} \int_0^\infty 1 - \left( 1 - \xi \exp\left( -z^2 \right)\right)^{n_e} dz
```
where ``n_e`` is the number of extrema, ``\xi`` is the ratio ``n_z/n_e`` with ``n_z`` being the number of zero crossings.

The integral is evaluated using Gauss-Legendre integration -- and is suitable for use within an automatic differentiation environment.

See also: [`peak_factor`](@ref), [`peak_factor_dk80`](@ref)
"""
function peak_factor_cl56(Dex::U, mi::SpectralMoments; nodes::Int=35) where {U<:Real}
	# get all necessary spectral moments
	m0 = mi.m0
	m2 = mi.m2
	m4 = mi.m4
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

	return sqrt(2.0) * int_sum
end


@doc raw"""
	peak_factor_cl56(n_z::T, n_e::T; nodes::Int=35) where T<:Real

Peak factor computed using the Cartwright and Longuet-Higgins (1956) formulation, using pre-computed `n_z` and `n_e` values.

Evaluates the integral to obtain the peak factor ``\psi``:
```math
\psi = \sqrt{2} \int_0^\infty 1 - \left( 1 - \xi \exp\left( -z^2 \right)\right)^{n_e} dz
```
where ``n_e`` is the number of extrema, ``\xi`` is the ratio ``n_z/n_e`` with ``n_z`` being the number of zero crossings.

The integral is evaluated using Gauss-Legendre integration -- and is suitable for use within an automatic differentiation environment.

See also: [`peak_factor`](@ref), [`peak_factor_dk80`](@ref)
"""
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


@doc raw"""
	peak_factor_cl56_gk(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Real,T<:Real}

Peak factor computed using the Cartwright and Longuet-Higgins (1956) formulation, using adaptive Gauss-Kronrod integration.

Evaluates the integral to obtain the peak factor ``\psi``:
```math
\psi = \sqrt{2} \int_0^\infty 1 - \left( 1 - \xi \exp\left( -z^2 \right)\right)^{n_e} dz
```
where ``n_e`` is the number of extrema, ``\xi`` is the ratio ``n_z/n_e`` with ``n_z`` being the number of zero crossings.

The integral is evaluated using Gauss-Kronrod integration via the QuadGK.jl package -- and is not suitable for use within an automatic differentiation environment.

See also: [`peak_factor`](@ref), [`peak_factor_dk80`](@ref)
"""
function peak_factor_cl56_gk(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Real,T<:Real}
	# get the numbers of zero crossing and extrema
	rvt = RandomVibrationParameters(:CL56)
	n_z, n_e = zeros_extrema_numbers(m, r_ps, fas, sdof, rvt)
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e

	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e

	return sqrt(2.0) * quadgk(integrand, 0.0, Inf)[1]
end


"""
	peak_factor_integrand_cl56(z, n_z, n_e)

Integrand of the Cartwright Longuet-Higgins (1956) peak factor formulation.

See also: [`peak_factor_cl56`](@ref)
"""
function peak_factor_integrand_cl56(z, n_z, n_e)
  ξ = n_z / n_e
  return 1.0 - (1.0 - ξ*exp(-z^2))^n_e
end


"""
	vanmarcke_cdf(x, n_z::T, δeff::T) where T<:Real

CDF of the Vanmarcke (1975) peak factor distribution -- used within the Der Kiureghian (1980) peak factor expression.
"""
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


"""
	vanmarcke_ccdf(x, n_z::T, δeff::T) where T<:Real

CCDF of the Vanmarcke (1975) peak factor distribution -- used within the Der Kiureghian (1980) peak factor expression.
"""
function vanmarcke_ccdf(x, n_z::T, δeff::T) where T<:Real
	return 1.0 - vanmarcke_cdf(x, n_z, δeff)
end


@doc raw"""
	peak_factor_dk80(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=30) where {S<:Real,T<:Real}

Peak factor computed using the Der Kiureghian (1980)/Vanmarcke (1975) formulation.

Evaluates the integral to obtain the peak factor ``\psi``:
```math
\psi = \int_0^\infty 1 - F_Z(z)dz
```
where ``F_Z(z)`` is the CDF defined by:
```math
F_Z(z) = \left[ 1-\exp\left(-\frac{z^2}{2}\right) \right] \exp\left[ -n_z \frac{1 - \exp\left( -\sqrt{\frac{\pi}{2}}\delta_{eff} z \right)}{\exp\left(\frac{z^2}{2}\right) - 1} \right]
```

The effective bandwidth ``\delta_{eff}=\delta^{1.2}`` where the bandwidth is:
```math
\delta = \sqrt{ 1 - \frac{m_1^2}{m_0 m_2} }
```

The integral is evaluated using Gauss-Legendre integration -- and is suitable for use within an automatic differentiation environment.

See also: [`peak_factor`](@ref), [`peak_factor_cl56`](@ref)
"""
function peak_factor_dk80(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; nodes::Int=30) where {S<:Real,T<:Real}
	# compute first three spectral moments
	mi = spectral_moments([0, 1, 2], m, r_ps, fas, sdof)
	m0 = mi.m0
	m1 = mi.m1
	m2 = mi.m2

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

	return dfac * dot( wi, yy )
end


@doc raw"""
	peak_factor_dk80(Dex::U, mi::SpectralMoments; nodes::Int=30) where {U<:Real}

Peak factor computed using the Der Kiureghian (1980)/Vanmarcke (1975) formulation, using precomputed `Dex` and `m0`.

Evaluates the integral to obtain the peak factor ``\psi``:
```math
\psi = \int_0^\infty 1 - F_Z(z)dz
```
where ``F_Z(z)`` is the CDF defined by:
```math
F_Z(z) = \left[ 1-\exp\left(-\frac{z^2}{2}\right) \right] \exp\left[ -n_z \frac{1 - \exp\left( -\sqrt{\frac{\pi}{2}}\delta_{eff} z \right)}{\exp\left(\frac{z^2}{2}\right) - 1} \right]
```

The effective bandwidth ``\delta_{eff}=\delta^{1.2}`` where the bandwidth is:
```math
\delta = \sqrt{ 1 - \frac{m_1^2}{m_0 m_2} }
```

The integral is evaluated using Gauss-Legendre integration -- and is suitable for use within an automatic differentiation environment.

See also: [`peak_factor`](@ref), [`peak_factor_cl56`](@ref)
"""
function peak_factor_dk80(Dex::U, mi::SpectralMoments; nodes::Int=30) where {U<:Real}
	# compute first three spectral moments
	m0 = mi.m0
	m1 = mi.m1
	m2 = mi.m2

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

	return dfac * dot( wi, yy )
end


@doc raw"""
	peak_factor_dk80_gk(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Float64,T<:Real}

Peak factor computed using the Der Kiureghian (1980)/Vanmarcke (1975) formulation, using precomputed `Dex` and `m0`.

Evaluates the integral to obtain the peak factor ``\psi``:
```math
\psi = \int_0^\infty 1 - F_Z(z)dz
```
where ``F_Z(z)`` is the CDF defined by:
```math
F_Z(z) = \left[ 1-\exp\left(-\frac{z^2}{2}\right) \right] \exp\left[ -n_z \frac{1 - \exp\left( -\sqrt{\frac{\pi}{2}}\delta_{eff} z \right)}{\exp\left(\frac{z^2}{2}\right) - 1} \right]
```

The effective bandwidth ``\delta_{eff}=\delta^{1.2}`` where the bandwidth is:
```math
\delta = \sqrt{ 1 - \frac{m_1^2}{m_0 m_2} }
```

The integral is evaluated using adaptive Gauss-Kronrod integration from the QuadGK.jl package -- and is therefore not suitable for use within an automatic differentiation environment.

See also: [`peak_factor`](@ref), [`peak_factor_dk80`](@ref), [`peak_factor_cl56`](@ref)
"""
function peak_factor_dk80_gk(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator) where {S<:Float64,T<:Real}
	# compute first three spectral moments
	mi = spectral_moments([0, 1, 2], m, r_ps, fas, sdof)
	m0 = mi.m0
	m1 = mi.m1
	m2 = mi.m2

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





@doc raw"""
	peak_factor(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator; rvt::RandomVibrationParameters) where {S<:Real,T<:Real}

Peak factor ``u_{max} / u_{rms}`` with a switch of `pf_method` to determine the approach adopted.
`rvt.pf_method` can currently be one of:
	- `:CL56` for Cartright Longuet-Higgins (1956)
	- `:DK80` for Der Kiureghian (1980), building on Vanmarcke (1975)

Defaults to `:DK80`.
"""
function peak_factor(m::S, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Real,T<:Real}
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
	peak_factor(Dex::U, m0::V, rvt::RandomVibrationParameters) where {U<:Real}

Peak factor u_max / u_rms with a switch of `pf_method` to determine the approach adopted.
`pf_method` can currently be one of:
	- `:CL56` for Cartright Longuet-Higgins (1956)
	- `:DK80` for Der Kiureghian (1980), building on Vanmarcke (1975)

Defaults to `:DK80`.
"""
function peak_factor(Dex::U, mi::SpectralMoments, rvt::RandomVibrationParameters) where {U<:Real}
	if rvt.pf_method == :CL56
		return peak_factor_cl56(Dex, mi)
	elseif rvt.pf_method == :DK80
		return peak_factor_dk80(Dex, mi)
	else
		return U(NaN)
	end
end
