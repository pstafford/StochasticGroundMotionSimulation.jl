



function spectral_moment(order::Int, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14, intervals::Int=101 ) where T<:Real
	# pre-allocate for squared fourier amplitude spectrum
	Af2 = zeros(Real,intervals)
	Hf2 = zeros(Real,intervals)
	# partition the integration domain to make sure the integral captures the key change in the integrand
	# the key frequencies are likely to be
	# - f_c (corner frequency)
	# - f_osc (oscillator frequency and two values either side)
	# - f_Q (combined kappa_0 & kappa_r fall-off frequency)
	f_n = sdof.f_n
	# note that we start slightly above zero to avoid a numerical issue with the frequency dependent Q(f) gradients
	fi = [ 1e-10, 0.1, 1.0, 10.0, f_n/1.1, f_n, f_n*1.1, 100.0, 300.0 ]
	sort!(fi)

	int_sum = 0.0
	for i in 1:length(fi)-1
		if fi[i] != fi[i+1]
			xx = collect(range(fi[i], stop=fi[i+1], length=intervals))
			squared_transfer!(Hf2, xx, sdof)
			squared_fourier_spectrum!(Af2, xx, m, r, fas; fc_fun=fc_fun, amp_model=amp_model)
			yy = (2π * xx).^order .* Hf2 .* Af2
			int_sum += simpsons_rule(xx,yy)
		end
	end
    return 2int_sum
end


function spectral_moments(order::Vector{Int}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14, intervals::Int=101 ) where T<:Real
	# pre-allocate for squared fourier amplitude spectrum
	Af2 = zeros(Real, intervals)
	Hf2 = zeros(Real, intervals)
	# partition the integration domain to make sure the integral captures the key change in the integrand
	# the key frequencies are likely to be
	# - f_c (corner frequency)
	# - f_osc (oscillator frequency and two values either side)
	# - f_Q (combined kappa_0 & kappa_r fall-off frequency)
	f_n = sdof.f_n
 	# note that we start slightly above zero to avoid a numerical issue with the frequency dependent Q(f) gradients
	fi = [ 1e-10, 0.1, 1.0, 10.0, f_n/1.1, f_n, f_n*1.1, 100.0, 300.0 ]
	sort!(fi)

	# make sure the orders are listed as increasing for the following loop approach
	sort!(order)
	dorder = diff(order)
	int_sumi = zeros(length(order))

	for i in 1:length(fi)-1
		if fi[i] != fi[i+1] # this can arise because of the mix of static and dynamic frequencies in the fi definition
			xx = collect(range(fi[i], stop=fi[i+1], length=intervals))
			squared_transfer!(Hf2, xx, sdof)
			squared_fourier_spectrum!(Af2, xx, m, r, fas; fc_fun=fc_fun, amp_model=amp_model)
			# compute default zeroth order integrand amplitude
			yy = Hf2 .* Af2

			for (idx, o) in enumerate(order)
				if idx == 1
					yy .*= (2π*xx).^o
				else
					yy .*= (2π*xx).^(dorder[idx-1])
				end
				@inbounds int_sumi[idx] += simpsons_rule(xx,yy)
			end
		end
	end
    return 2.0 * int_sumi
end



# function to compute spectral moments (m_0, m_2 & m_4)
function spectral_moments( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14, intervals::Int=101 ) where T<:Real
	# pre-allocate for squared fourier amplitude spectrum
	Af2 = zeros(Real, intervals)
	Hf2 = zeros(Real, intervals)
	# partition the integration domain to make sure the integral captures the key change in the integrand
	# the key frequencies are likely to be
	# - f_c (corner frequency)
	# - f_osc (oscillator frequency and two values either side)
	# - f_Q (combined kappa_0 & kappa_r fall-off frequency)
	f_n = sdof.f_n
 	# note that we start slightly above zero to avoid a numerical issue with the frequency dependent Q(f) gradients
	fi = [ 1e-10, 0.1, 1.0, 10.0, f_n/1.1, f_n, f_n*1.1, 100.0, 300.0 ]
	sort!(fi)

	int_sum0 = 0.0
	int_sum2 = 0.0
	int_sum4 = 0.0
	for i in 1:length(fi)-1
		if fi[i] != fi[i+1]
			xx = collect(range(fi[i], stop=fi[i+1], length=intervals))
			squared_transfer!(Hf2, xx, sdof)
			squared_fourier_spectrum!(Af2, xx, m, r, fas; fc_fun=fc_fun, amp_model=amp_model )
			yy = Hf2 .* Af2
			int_sum0 += simpsons_rule(xx,yy)
			yy .*= (2π*xx).^2
			int_sum2 += simpsons_rule(xx,yy)
			yy .*= (2π*xx).^2
			int_sum4 += simpsons_rule(xx,yy)
		end
	end
    return 2int_sum0, 2int_sum2, 2int_sum4
end


function spectral_moments_cy( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14, intervals::Int=101 ) where T<:Real
	# pre-allocate for squared fourier amplitude spectrum
	Af2 = zeros(Real, intervals)
	Hf2 = zeros(Real, intervals)
	# partition the integration domain to make sure the integral captures the key change in the integrand
	f_n = sdof.f_n
 	# note that we start slightly above zero to avoid a numerical issue with the frequency dependent Q(f) gradients
	fi = [ 1e-10, 0.1, 1.0, 10.0, f_n/1.1, f_n, f_n*1.1, 100.0, 300.0 ]
	sort!(fi)

	int_sum0 = 0.0
	int_sum2 = 0.0
	int_sum4 = 0.0
	for i in 1:length(fi)-1
		if fi[i] != fi[i+1]
			xx = collect(range(fi[i], stop=fi[i+1], length=intervals))
			squared_transfer!(Hf2, xx, sdof)
			squared_fourier_spectrum_cy!(Af2, xx, m, r, fas; fc_fun=fc_fun, amp_model=amp_model )
			yy = Hf2 .* Af2
			int_sum0 += simpsons_rule(xx,yy)
			yy .*= (2π*xx).^2
			int_sum2 += simpsons_rule(xx,yy)
			yy .*= (2π*xx).^2
			int_sum4 += simpsons_rule(xx,yy)
		end
	end
    return 2int_sum0, 2int_sum2, 2int_sum4
end


# function to define functions to compute the extrema and zero-crossing frequencies
function zeros_extrema_frequencies( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	m_0, m_2, m_4 = spectral_moments( m, r, fas, sdof; fc_fun=fc_fun, amp_model=amp_model )

	fzero = sqrt( m_2 / m_0 ) / (2π)
	fextrema = sqrt( m_4 / m_2 ) / (2π)
	return fzero, fextrema
end


# function to define the numbers of zero crossings and extrema
function zeros_extrema_numbers( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
    (fz, fe) = zeros_extrema_frequencies( m, r, fas, sdof; fc_fun=fc_fun, amp_model=amp_model )
    Dex = boore_thompson_2014( m, r, fas; fc_fun=fc_fun )
    return (2fz*Dex, 2fe*Dex)
end





function rvt_response_spectral_ordinate( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
  # get all relevant spectral moments
  m_0, m_2, m_4 = spectral_moments( m, r, fas, sdof; fc_fun=fc_fun, amp_model=amp_model )
  # get the duration metrics
  (Drms, Dex, Dratio) = boore_thompson_2012( m, r, fas, sdof; region=region, fc_fun=fc_fun )
  # get the rms response
  y_rms = sqrt( m_0 / Drms )
  # get the peak factor
  n_z = Dex * sqrt( m_2 / m_0 ) / π
  n_e = Dex * sqrt( m_4 / m_2 ) / π
  pf = peak_factor( n_z, n_e )
  # response spectral value (in g)
  Sa = pf * y_rms / 9.80665
  return Sa
end

function rvt_response_spectral_ordinate_cy( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, sdof::Oscillator{T}; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
  # get all relevant spectral moments
  m_0, m_2, m_4 = spectral_moments_cy( m, r, fas, sdof; fc_fun=fc_fun, amp_model=amp_model )
  # get the duration metrics
  (Drms, Dex, Dratio) = boore_thompson_2012( m, r, fas, sdof; region=region, fc_fun=fc_fun )
  # get the rms response
  y_rms = sqrt( m_0 / Drms )
  # get the peak factor
  n_z = Dex * sqrt( m_2 / m_0 ) / π
  n_e = Dex * sqrt( m_4 / m_2 ) / π
  pf = peak_factor( n_z, n_e )
  # response spectral value (in g)
  Sa = pf * y_rms / 9.80665
  return Sa
end


function rvt_response_spectral_ordinate( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, period::T; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
  	# create a sdof instance
  	sdof = Oscillator(1.0/period)
	# call other function
	return rvt_response_spectral_ordinate(m, r, fas, sdof; region=region, fc_fun=fc_fun, amp_model=amp_model)
end

function rvt_response_spectral_ordinate_cy( m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}, period::T; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
  	# create a sdof instance
  	sdof = Oscillator(1.0/period)
	# call other function
	return rvt_response_spectral_ordinate_cy(m, r, fas, sdof; region=region, fc_fun=fc_fun, amp_model=amp_model)
end


function rvt_response_spectrum( period::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	Sa = zeros(typeof(fas.Δσ),length(period))
	for i = 1:length(period)
		Sa[i] = rvt_response_spectral_ordinate( m, r, fas, period[i]; region=region, fc_fun=fc_fun, amp_model=amp_model )
  	end
  	return Sa
end

function rvt_response_spectrum_cy( period::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	Sa = zeros(typeof(fas.Δσ),length(period))
	for i = 1:length(period)
		Sa[i] = rvt_response_spectral_ordinate_cy( m, r, fas, period[i]; region=region, fc_fun=fc_fun, amp_model=amp_model )
  	end
  	return Sa
end

function rvt_response_spectrum!( Sa::Vector, period::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	for i = 1:length(period)
		Sa[i] = rvt_response_spectral_ordinate( m, r, fas, period[i]; region=region, fc_fun=fc_fun, amp_model=amp_model )
  	end
end

function rvt_response_spectrum_cy!( Sa::Vector, period::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; region::Symbol=:WNA, fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	for i = 1:length(period)
		Sa[i] = rvt_response_spectral_ordinate_cy( m, r, fas, period[i]; region=region, fc_fun=fc_fun, amp_model=amp_model )
  	end
end
