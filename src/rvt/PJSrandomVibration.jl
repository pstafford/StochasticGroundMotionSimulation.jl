



function spectral_moment(order::Int, m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; intervals::Int=101, control_freqs::Vector{Float64}=[1e-3, 1e-1, 1.0, 10.0, 100.0, 300.0] ) where T<:Float64
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
	fidLO = findall(control_freqs .< f_n/1.5)
	fidHI = findall(control_freqs .> f_n*1.5)
	fi = [ control_freqs[fidLO]; f_n/1.5; f_n*1.5; control_freqs[fidHI] ]

	int_sum = 0.0
	for i in 1:length(fi)-1
		if fi[i] != fi[i+1]
			xx = collect(range(fi[i], stop=fi[i+1], length=intervals))
			squared_transfer!(Hf2, xx, sdof)
			squared_fourier_spectrum!(Af2, xx, m, r, fas)
			yy = @. (2π * xx)^order * Hf2 * Af2
			# int_sum += simpsons_rule(xx,yy)
			int_sum += trapezoidal_rule(xx,yy)
		end
	end
    return 2int_sum
end


function spectral_moments(order::Vector{Int}, m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; intervals::Int=101, control_freqs::Vector{Float64}=[1e-3, 1e-1, 1.0, 10.0, 100.0, 300.0] ) where T<:Float64
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
	# fi = [ 1e-10, 0.1, 1.0, 10.0, f_n/1.1, f_n, f_n*1.1, 100.0, 300.0 ]
	fidLO = findall(control_freqs .< f_n/1.5)
	fidHI = findall(control_freqs .> f_n*1.5)
	fi = [ control_freqs[fidLO]; f_n/1.5; f_n*1.5; control_freqs[fidHI] ]

	# make sure the orders are listed as increasing for the following loop approach
	sort!(order)
	dorder = diff(order)
	int_sumi = zeros(length(order))

	for i in 1:length(fi)-1
		if fi[i] != fi[i+1] # this can arise because of the mix of static and dynamic frequencies in the fi definition
			xx = collect(range(fi[i], stop=fi[i+1], length=intervals))
			squared_transfer!(Hf2, xx, sdof)
			squared_fourier_spectrum!(Af2, xx, m, r, fas)
			# compute default zeroth order integrand amplitude
			yy = @. Hf2 * Af2

			for (idx, o) in enumerate(order)
				if idx == 1
					yy .*= (2π*xx).^o
				else
					@inbounds yy .*= (2π*xx).^(dorder[idx-1])
				end
				# @inbounds int_sumi[idx] += simpsons_rule(xx,yy)
				@inbounds int_sumi[idx] += trapezoidal_rule(xx,yy)
			end
		end
	end
    return 2.0 * int_sumi
end

function spectral_moments_ln(order::Vector{Int}, m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; intervals::Int=101, control_freqs::Vector{Float64}=[1e-3, 1e-1, 1.0, 10.0, 100.0, 300.0] ) where T<:Float64
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
	fidLO = findall(control_freqs .< f_n/1.5)
	fidHI = findall(control_freqs .> f_n*1.5)
	fi = [ control_freqs[fidLO]; f_n/1.5; f_n*1.5; control_freqs[fidHI] ]
	lnfi = log.(fi)

	# make sure the orders are listed as increasing for the following loop approach
	sort!(order)
	dorder = diff(order)
	int_sumi = zeros(length(order))

	for i in 1:length(lnfi)-1
		if lnfi[i] != lnfi[i+1] # this can arise because of the mix of static and dynamic frequencies in the fi definition
			uu = collect(range(lnfi[i], stop=lnfi[i+1], length=intervals))
			xx = exp.(uu)
			squared_transfer!(Hf2, xx, sdof)
			squared_fourier_spectrum!(Af2, xx, m, r, fas)
			# compute default zeroth order integrand amplitude
			yy = @. Hf2 * Af2 * xx

			for (idx, o) in enumerate(order)
				if idx == 1
					yy .*= (2π*xx).^o
				else
					@inbounds yy .*= (2π*xx).^(dorder[idx-1])
				end
				# @inbounds int_sumi[idx] += simpsons_rule(uu,yy)
				@inbounds int_sumi[idx] += trapezoidal_rule(uu,yy)
			end
		end
	end
    return 2.0 * int_sumi
end

function spectral_moments_gk(order::Vector{Int}, m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; intervals::Int=101 ) where T<:Float64
	int_sumi = zeros(length(order))

	for (i, o) in enumerate(order)
		moment_integrand(f) = squared_transfer(f, sdof) * ( fourier_spectral_ordinate(f, m, r, fas)^2 ) * (2π*f)^o
		int_sumi[i] = quadgk(moment_integrand, 0, sdof.f_n, Inf)[1]
	end

    return 2.0 * int_sumi
end


# function to define functions to compute the extrema and zero-crossing frequencies
function zeros_extrema_frequencies(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}) where T<:Float64
	m_0, m_2, m_4 = spectral_moments(m, r, fas, sdof)
	fzero = sqrt( m_2 / m_0 ) / (2π)
	fextrema = sqrt( m_4 / m_2 ) / (2π)
	return fzero, fextrema
end


# function to define the numbers of zero crossings and extrema
function zeros_extrema_numbers(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}) where T<:Float64
    (fz, fe) = zeros_extrema_frequencies(m, r, fas, sdof)
    Dex = boore_thompson_2014(m, r, fas)
    return (2fz*Dex, 2fe*Dex)
end




function rvt_response_spectral_ordinate(m::T, r::T, fas::FourierParameters, sdof::Oscillator{T}; region::Symbol=:WNA) where T<:Float64
	# get all relevant spectral moments
	m_0, m_2, m_4 = spectral_moments(m, r, fas, sdof)
	# get the duration metrics
	(Drms, Dex, Dratio) = boore_thompson_2012(m, r, fas, sdof; region=region)
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


function rvt_response_spectral_ordinate(period::T, m::T, r::T, fas::FourierParameters; region::Symbol=:WNA) where T<:Float64
  	# create a sdof instance
  	sdof = Oscillator(1.0/period)
	# call other function
	return rvt_response_spectral_ordinate(m, r, fas, sdof; region=region)
end



function rvt_response_spectrum(period::Vector{T}, m::T, r::T, fas::FourierParameters; region::Symbol=:WNA) where T<:Float64
	Sa = zeros(typeof(fas.source.Δσ),length(period))
	for i = 1:length(period)
		@inbounds Sa[i] = rvt_response_spectral_ordinate(period[i], m, r, fas; region=region)
  	end
  	return Sa
end



function rvt_response_spectrum!(Sa::Vector, period::Vector{T}, m::T, r::T, fas::FourierParameters; region::Symbol=:WNA) where T<:Float64
	for i = 1:length(period)
		Sa[i] = rvt_response_spectral_ordinate(period[i], m, r, fas; region=region)
  	end
end
