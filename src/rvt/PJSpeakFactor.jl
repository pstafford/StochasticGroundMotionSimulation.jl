


# function to compute the peak over rms ratio, a.k.a. peak factor
function peak_factor_cl56(m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14) where T<:Real
	# get the numbers of zero crossing and extrema
	n_z, n_e = zeros_extrema_numbers( m, r, fas, sdof; fc_fun=fc_fun, amp_model=amp_model )
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e
	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e
	# compute the peak factor
	xx = collect(range(0.0, stop=10.0, length=101))
	yy = integrand.(xx)
	int_sum = simpsons_rule(xx, yy)

	pf = sqrt(2.0) * int_sum
	return pf
end

function peak_factor_cl56(n_z::T, n_e::T) where T<:Real
	# get the ratio of zero crossings to extrema
	ξ = n_z / n_e
	# define the integrand
	integrand(z) = 1.0 - (1.0 - ξ*exp(-z^2))^n_e
	# compute the peak factor
	xx = collect(range(0.0, stop=10.0, length=101))
	yy = integrand.(xx)
	int_sum = simpsons_rule(xx, yy)

	pf = sqrt(2.0) * int_sum
	return pf
end


function peak_factor_integrand_cl56(z, n_z, n_e)
  ξ = n_z / n_e
  return 1.0 - (1.0 - ξ*exp(-z^2))^n_e
end

function peak_factor_dk80(m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14) where T<:Real
	# get the numbers of zero crossing and extrema
	n_z, n_e = zeros_extrema_numbers(m, r, fas, sdof; fc_fun=fc_fun, amp_model=amp_model)

end
