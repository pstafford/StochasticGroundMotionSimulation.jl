
"""
	fourier_constant(src::SourceParameters)

Define the constant source term for the Fourier Amplitude Spectrum.
Constant set to permit distances to be passed in km, densities in t/m^3, and velocities in km/s
"""
function fourier_constant(src::SourceParameters)
	# note that the 1e-20 factor allows one to pass in distances in km, density in t/m^3, velocity in km/s
	Rref0 = 1.0
	return src.RΘϕ * src.V * src.F * 1e-20 / ( 4π * src.ρ * src.β^3 * Rref0 )
end

fourier_constant(fas::FourierParameters) = fourier_constant(fas.source)

"""
	fourier_source(f::T, m::T, src::SourceParameters) where T<:Float64

Source Fourier Amplitude Spectrum of displacement, without the constant term.
This simply includes the seismic moment and the source spectral shape.
"""
function fourier_source(f::T, m::T, src::SourceParameters) where T<:Float64
	Mo = magnitude_to_moment(m)
	(fa, fb, ε) = corner_frequency(m, src)
	if src.source_model == :Brune
		# use single corner, with fa=fc
		return Mo / (1.0 + (f/fa)^2)
	elseif src.source_model == :Atkinson_Silva_2000
		# use double corner model
		return Mo * ( (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2) )
	else
		# use single corner, with fa=fc
		return Mo / (1.0 + (f/fa)^2)
	end
end

fourier_source(f::T, m::T, fas::FourierParameters) where T = fourier_source(f, m, fas.source)

function fourier_source_shape(f::T, m::T, src::SourceParameters) where T<:Real
	(fa, fb, ε) = corner_frequency(m, fas)
	if src.source_model == :Brune
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	elseif src.source_model == :Atkinson_Silva_2000
		# use double corner model
		return (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2)
	else
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	end
end

fourier_source_shape(f::T, m::T, fas::FourierParameters) where T = fourier_source_shape(f, m, fas.source)

function fourier_source_shape(f::Float64, fa::Union{Float64,Dual{Float64}}, fb::Union{Float64,Dual{Float64}}, ε::Union{Float64,Dual{Float64}})
	if isnan(ε) == true
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	else
		# use double corner model
		return (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2)
	end
end



function fourier_path(f::T, r::T, path::PathParameters) where T<:Float64
	Zr = geometric_spreading(r, path)
	return Zr * exp( -π*f*r / ( path.Q0 * f^path.η * path.cQ ) )
end

fourier_path(f, r, fas::FourierParameters) = fourier_path(f, r, fas.path)


function fourier_attenuation(f::T, r::T, path::PathParameters, site::SiteParameters) where T<:Float64
	if f < eps()
		return 1.0
	else
		return exp( -π*f*r / ( path.Q0 * f^path.η * path.cQ ) - π*f*site.κ0 )
	end
end

fourier_attenuation(f::T, r::T, fas::FourierParameters) where T<:Float64 = fourier_attenuation(f, r, fas.path, fas.site)


function fourier_site(f::Float64, site::SiteParameters)
	# kappa filter
    Sf = exp( -π*f*site.κ0 )
    # site impedence function
    Sif = site_amplification(f, site)
	return Sf*Sif
end

fourier_site(f::Float64, fas::FourierParameters) = fourier_site(f, fas.site)


function fourier_spectral_ordinate(f::T, m::T, r::T, fas::FourierParameters) where T<:Float64
	# define all constant terms here
  	C = fourier_constant(fas)
	# source term
	Ef = fourier_source(f, m, fas)
	# use point-source distance metric
	h = near_source_saturation(m, fas)
	r_ps = sqrt(r^2 + h^2)
	# geometric spreading
	Gr = geometric_spreading(r_ps, fas)
	# combined attenuation
	Kf = fourier_attenuation(f, r_ps, fas)
	# site impedance
    Sf = site_amplification(f, fas)

	# overall displacement spectrum (in cms) is Ef*Pf*Sf
	# the division by 10^2 below is to convert from units of cm/s to m/s
	return (2π*f)^2 * C * Ef * Gr * Kf * Sf / 100.0
end



function fourier_spectrum(f::Vector{T}, m::T, r::T, fas::FourierParameters) where T<:Float64
	# define all frequency independent terms here
  	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = near_source_saturation(m, fas)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading(r_ps, fas)
	factor = 4π^2 * C * Mo * Gr / 100.0

	Af = ones(typeof(fas.source.Δσ), size(f))
  	for i = 1:length(f)
	  	# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i], fas)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af[i] = Ef * Kf * Sf * factor * f[i]^2
    end
    return Af
end



function fourier_spectrum!(Af::Vector{Union{T,Dual{T}}}, f::Vector{T}, m::T, r::T, fas::FourierParameters) where T<:Float64
	# define all frequency independent terms here
	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas)
	# use point-source distance metric
	h = near_source_saturation(m, fas)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading(r_ps, fas)
	factor = 4π^2 * C * Mo * Gr / 100.0

	for i = 1:length(f)
		# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i], fas)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af[i] = Ef * Kf * Sf * factor * f[i]^2
  	end
  	nothing
end



function squared_fourier_spectrum!(Af2::Vector{Union{T,Dual{T}}}, f::Vector{T}, m::T, r::T, fas::FourierParameters) where T<:Float64
	# define all frequency independent terms here
	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas)
	# use point-source distance metric
	h = near_source_saturation(m, fas)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading(r_ps, fas)
	factor = (4π^2 * C * Mo * Gr / 100.0)^2

	for i = 1:length(f)
		# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i], fas)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af2[i] = (Ef * Kf * Sf)^2 * factor * f[i]^4
  	end
  	nothing
end



function combined_kappa_frequency(r::Float64, fas::FourierParameters)
	# function gives the frequency at which the combined κ_r and κ_0 filters (squared versions) give a value of 0.5
  	if fas.path.η < 0.1
    	# a closed form solution exists
    	return log(2.0) / ( 2π * ( r / (fas.path.Q0 * fas.path.cQ) + fas.site.κ0 ) )
  	else
		g(f) = 0.5 - (fourier_attenuation(f, r, fas)^2)
		f_0 = find_zero(g, (0.01,100.0), Bisection(); xatol=1e-2)
		fk = min(max(f_0, 0.2), 1.0)
	    return fk
  	end
end
