
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
	fourier_source_shape(f::T, m::T, src::SourceParameters) where T<:Float64

Source shape of the Fourier Amplitude Spectrum of _displacement_, without the constant term or seismic moment.
This simply includes the source spectral shape.
"""
function fourier_source_shape(f::T, m::T, src::SourceParameters) where T<:Float64
	fa, fb, ε = corner_frequency(m, src)
	if src.model == :Brune
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	elseif src.model == :Atkinson_Silva_2000
		# use double corner model
		return (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2)
	else
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	end
end

fourier_source_shape(f::T, m::T, fas::FourierParameters) where T = fourier_source_shape(f, m, fas.source)


function fourier_source_shape(f::Float64, fa::T, fb::T, ε::T, model::Symbol) where T<:Real
	if model == :Brune
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	elseif model == :Atkinson_Silva_2000
		# use double corner model
		return (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2)
	else
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	end
end



"""
	fourier_source(f::T, m::T, src::SourceParameters) where T<:Float64

Source Fourier Amplitude Spectrum of displacement, without the constant term.
This simply includes the seismic moment and the source spectral shape.
"""
function fourier_source(f::T, m::T, src::SourceParameters) where T<:Float64
	Mo = magnitude_to_moment(m)
	return Mo * fourier_source_shape(f, m, src)
end

fourier_source(f::T, m::T, fas::FourierParameters) where T = fourier_source(f, m, fas.source)





function fourier_path(f::S, r_ps::T, geo::GeometricSpreadingParameters, ane::AnelasticAttenuationParameters) where {S<:Float64,T<:Real}
	Zr = geometric_spreading(r_ps, geo)
	Qr = anelastic_attenuation(f, r_ps, ane)
	return Zr * Qr
end

fourier_path(f, r, path::PathParameters) = fourier_path(f, r, path.geometric, path.anelastic)
fourier_path(f, r, fas::FourierParameters) = fourier_path(f, r, fas.path)


"""
	fourier_attenuation(f::S, r_ps::T, ane::AnelasticAttenuationParameters{U,V}, site::SiteParameters{W}) where {S<:Float64,T<:Real,U<:Real,V<:Real,W<:Real}

Combined full-path attenuation, including `Q(f)` effects and `κ0` filter for frequency `f`
Distance defined in terms of an equivalent point source distance `r_ps`
"""
function fourier_attenuation(f::S, r_ps::T, ane::AnelasticAttenuationParameters{U,V}, site::SiteParameters{W}) where {S<:Float64,T<:Real,U<:Real,V<:Real,W<:Real}
	if f < eps()
		return oneunit(T)*oneunit(U)*oneunit(V)*oneunit(W)
	else
		Qr = anelastic_attenuation(f, r_ps, ane)
		Kf = kappa_filter(f, site)
		return Qr * Kf
	end
end

fourier_attenuation(f::S, r_ps::T, path::PathParameters, site::SiteParameters) where {S<:Float64,T<:Real} = fourier_attenuation(f, r_ps, path.anelastic, site)
fourier_attenuation(f::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real} = fourier_attenuation(f, r_ps, fas.path, fas.site)


"""
	fourier_site(f::Float64, site::SiteParameters)

Combined site amplification and kappa filter for frequency `f`
"""
function fourier_site(f::Float64, site::SiteParameters)
	# kappa filter
    Kf = kappa_filter(f, site)
    # site impedance function
    Sf = site_amplification(f, site)
	return Sf * Kf
end

fourier_site(f::Float64, fas::FourierParameters) = fourier_site(f, fas.site)



"""
	fourier_spectral_ordinate(f::S, m::S, r_ps::T, src::SourceParameters, geo::GeometricSpreadingParameters, ane::AnelasticAttenuationParameters, site::SiteParameters) where {S<:Float64,T<:Real}

Fourier acceleration spectral ordinate (m/s) based upon an equivalent point source distance `r_ps`
- `f` is frequency (Hz)
- `m` is magnitude
- `r_ps` is the equivalent point source distance including saturation effects (km)
- `src` are the source parameters `SourceParameters`
- `geo` are the geometric spreading parameters `GeometricSpreadingParameters`
- `ane` are the anelastic attenuation parameters `AnelasticAttenuationParameters`
- `site` are the site parameters `SiteParameters`
"""
function fourier_spectral_ordinate(f::S, m::S, r_ps::T, src::SourceParameters, geo::GeometricSpreadingParameters, ane::AnelasticAttenuationParameters, site::SiteParameters) where {S<:Float64,T<:Real}
	# define all constant terms here
  	C = fourier_constant(src)
	# source term
	Ef = fourier_source(f, m, src)
	# geometric spreading
	Gr = geometric_spreading(r_ps, geo)
	# combined attenuation of both path (κr) and site (κ0)
	Kf = fourier_attenuation(f, r_ps, ane, site)
	# site impedance
    Sf = site_amplification(f, site)

	# overall displacement spectrum (in cms) is Ef*Pf*Sf
	# the division by 10^2 below is to convert from units of cm/s to m/s
	return (2π*f)^2 * C * Ef * Gr * Kf * Sf / 100.0
end

"""
	fourier_spectral_ordinate(f::S, m::S, r_ps::T, src::SourceParameters, path::PathParameters, site::SiteParameters) where {S<:Float64,T<:Real}

Fourier acceleration spectral ordinate (m/s) based upon an equivalent point source distance `r_ps`
- `f` is frequency (Hz)
- `m` is magnitude
- `r_ps` is the equivalent point source distance including saturation effects (km)
- `src` are the source parameters `SourceParameters`
- `path` are the path parameters `PathParameters`
- `site` are the site parameters `SiteParameters`
"""
fourier_spectral_ordinate(f::S, m::S, r_ps::T, src::SourceParameters, path::PathParameters, site::SiteParameters) where {S<:Float64,T<:Real} = fourier_spectral_ordinate(f, m, r_ps, src, path.geometric, path.anelastic, site)

"""
	fourier_spectral_ordinate(f::S, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real}

Fourier acceleration spectral ordinate (m/s) based upon an equivalent point source distance `r_ps`
- `f` is frequency (Hz)
- `m` is magnitude
- `r_ps` is the equivalent point source distance including saturation effects (km)
- `fas` are the Fourier spectral parameters `FourierParameters`
"""
fourier_spectral_ordinate(f::S, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real} = fourier_spectral_ordinate(f, m, r_ps, fas.source, fas.path, fas.site)


"""
	squared_fourier_spectral_ordinate(f::S, m::S, r_ps::T, src::SourceParameters, geo::GeometricSpreadingParameters, ane::AnelasticAttenuationParameters, site::SiteParameters) where {S<:Float64,T<:Real}

Squared Fourier acceleration spectral ordinate (m^2/s^2) based upon an equivalent point source distance `r_ps`
- `f` is frequency (Hz)
- `m` is magnitude
- `r_ps` is the equivalent point source distance including saturation effects (km)
- `src` are the source parameters `SourceParameters`
- `geo` are the geometric spreading parameters `GeometricSpreadingParameters`
- `ane` are the anelastic attenuation parameters `AnelasticAttenuationParameters`
- `site` are the site parameters `SiteParameters`
"""
function squared_fourier_spectral_ordinate(f::S, m::S, r_ps::T, src::SourceParameters, geo::GeometricSpreadingParameters, ane::AnelasticAttenuationParameters, site::SiteParameters) where {S<:Float64,T<:Real}
	return fourier_spectral_ordinate(f, m, r_ps, src, geo, ane, site)^2
end

squared_fourier_spectral_ordinate(f::S, m::S, r_ps::T, src::SourceParameters, path::PathParameters, site::SiteParameters) where {S<:Float64,T<:Real} = squared_fourier_spectral_ordinate(f, m, r_ps, src, path.geometric, path.anelastic, site)
squared_fourier_spectral_ordinate(f::S, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real} = squared_fourier_spectral_ordinate(f, m, r_ps, fas.source, fas.path, fas.site)



"""
	fourier_spectrum(f::Vector{S}, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real}

Fourier acceleration spectrum (m/s) based upon an equivalent point source distance `r_ps`
- `f` is `Vector` of frequencies (Hz)
- `m` is magnitude
- `r_ps` is the equivalent point source distance including saturation effects (km)
- `fas` are the Fourier spectral parameters `FourierParameters`
"""
function fourier_spectrum(f::Vector{S}, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real}
	numf = length(f)
	U = get_parametric_type(fas)
	if numf == 0
		return Vector{U}()
	else
		# define all frequency independent terms here
	  	C = fourier_constant(fas)
		# source amplitude and corner frequencies
		Mo = magnitude_to_moment(m)
		fa, fb, ε = corner_frequency(m, fas)
		# geometric spreading
		Gr = geometric_spreading(r_ps, fas)
		factor = 4π^2 * C * Mo * Gr / 100.0

		Af = Vector{U}(undef, numf)
	  	for i in 1:numf
			@inbounds fi = f[i]
		  	# source term
			Ef = fourier_source_shape(fi, fa, fb, ε, fas.source.model)
	  		# combined attenuation
	  		Kf = fourier_attenuation(fi, r_ps, fas)
	  		# site impedance
	      	Sf = site_amplification(fi, fas)
			# apply factor and convert to acceleration in appropriate units
	  		@inbounds Af[i] = Ef * Kf * Sf * factor * fi^2
	    end
	    return Af
	end
end

function squared_fourier_spectrum(f::Vector{S}, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real}
	Af = fourier_spectrum(f, m, r_ps, fas)
	return Af.^2
end


"""
	fourier_spectrum!(Af::Vector{U}, f::Vector{S}, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real,U<:Real}

Fourier acceleration spectrum (m/s) based upon an equivalent point source distance `r_ps`
- `Af` is the vector of fas amplitudes to be filled (m/s)
- `f` is `Vector` of frequencies (Hz)
- `m` is magnitude
- `r_ps` is the equivalent point source distance including saturation effects (km)
- `fas` are the Fourier spectral parameters `FourierParameters`
"""
function fourier_spectrum!(Af::Vector{U}, f::Vector{S}, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real,U<:Real}
	numf = length(f)
	numAf = length(Af)
	numf == numAf || error("length of `f` and `Af` must be equal")
	if numf == 0
		Af = Vector{U}()
	else
		# define all frequency independent terms here
	  	C = fourier_constant(fas)
		# source amplitude and corner frequencies
		Mo = magnitude_to_moment(m)
		fa, fb, ε = corner_frequency(m, fas)
		# geometric spreading
		Gr = geometric_spreading(r_ps, fas)
		factor = 4π^2 * C * Mo * Gr / 100.0

		for i in 1:numf
			@inbounds fi = f[i]
			# source term
			Ef = fourier_source_shape(fi, fa, fb, ε, fas.source.model)
	  		# combined attenuation
	  		Kf = fourier_attenuation(fi, r_ps, fas)
	  		# site impedance
	      	Sf = site_amplification(fi, fas)
			# apply factor and convert to acceleration in appropriate units
	  		@inbounds Af[i] = Ef * Kf * Sf * factor * fi^2
	  	end
	end
	return nothing
end



"""
	squared_fourier_spectrum!(Afsq::Vector{U}, f::Vector{S}, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real,U<:Real}

Fourier acceleration spectrum (m/s) based upon an equivalent point source distance `r_ps`
- `Afsq` is the vector of squared fas amplitudes to be filled (m^2/s^2)
- `f` is `Vector` of frequencies (Hz)
- `m` is magnitude
- `r_ps` is the equivalent point source distance including saturation effects (km)
- `fas` are the Fourier spectral parameters `FourierParameters`
"""
function squared_fourier_spectrum!(Afsq::Vector{U}, f::Vector{S}, m::S, r_ps::T, fas::FourierParameters) where {S<:Float64,T<:Real,U<:Real}
	numf = length(f)
	numAfsq = length(Afsq)
	numf == numAfsq || error("length of `f` and `Afsq` must be equal")
	if numf == 0
		Af = Vector{U}()
	else
		# define all frequency independent terms here
	  	C = fourier_constant(fas)
		# source amplitude and corner frequencies
		Mo = magnitude_to_moment(m)
		fa, fb, ε = corner_frequency(m, fas)
		# geometric spreading
		Gr = geometric_spreading(r_ps, fas)
		factor = (4π^2 * C * Mo * Gr / 100.0)^2

		for i in 1:numf
			@inbounds fi = f[i]
			# source term
			Ef = fourier_source_shape(fi, fa, fb, ε, fas.source.model)
	  		# combined attenuation
	  		Kf = fourier_attenuation(fi, r_ps, fas)
	  		# site impedance
	      	Sf = site_amplification(fi, fas)
			# apply factor and convert to acceleration in appropriate units
	  		@inbounds Afsq[i] = (Ef * Kf * Sf)^2 * factor * fi^4
	  	end
	end
	return nothing
end

"""
	combined_kappa_frequency(r_ps::T, ane::AnelasticAttenuationParameters, site::SiteParameters) where T<:Real

Frequency at which the combined κ_r and κ_0 filters (squared versions) give a value of 0.5
"""
function combined_kappa_frequency(r_ps::T, ane::AnelasticAttenuationParameters, site::SiteParameters) where T<:Real
  	if ane.η < 0.1
    	# a closed form solution exists (for effectively η=0)
    	return log(2.0) / ( 2π * ( r_ps / (ane.Q0 * ane.cQ) + site.κ0 ) )
  	else
		g(f) = 0.5 - (fourier_attenuation(f, r_ps, ane, site)^2)
		f_0 = find_zero(g, (0.01,100.0), Bisection(); xatol=1e-2)
		fk = min(max(f_0, 0.2), 1.0)
		U = get_parametric_type(ane)
		V = get_parametric_type(site)
		if T <: Float64
			if U <: Float64
				unit = oneunit(V)
			else
				unit = oneunit(U)
			end
		else
			unit = oneunit(T)
		end
	    return fk * unit
  	end
end

combined_kappa_frequency(r_ps::T, path::PathParameters, site::SiteParameters) where T<:Real = combined_kappa_frequency(r_ps, path.anelastic, site)
combined_kappa_frequency(r_ps::T, fas::FourierParameters) where T<:Real = combined_kappa_frequency(r_ps, fas.path, fas.site)
