
"""
    FASParams

Custom type for the parameters the Fourier amplitude spectrum. Note that type doesn't include the impedance effects. The parameters defined as `Union{Real, ForwardDiff.Dual{Real}}` are differentiable.

# Source parameters
    - `Δσ::Union{Real,ForwardDiff.Dual{Real}}` is the stress parameter in bars
    - `RΘϕ::Float64` is the radiation pattern
    - `V::Float64` is the partition factor (for splitting to horizontal components)
    - `F::Float64` is the free surface factor
    - `β::Float64` is the source velocity in units of km/s
    - `ρ::Float64` is the source density in units of t/m^3 or g/cm^3

# Path parameters
    - `Rrefi::Vector{Float64}` is a vector of reference distances for geometric spreading in units of km
    - `γi::Vector{Union{Real,ForwardDiff.Dual{Real}}}` is a vector of geometric spreading rates
    - `Q0::Union{Real,ForwardDiff.Dual{Real}}` is the Quality Constant
    - `η::Union{Real,ForwardDiff.Dual{Real}}` is the Quality Exponent
    - `cQ::Float64` is the velocity used to determine Q(f) (or the travel time) in units of km/s

# Site parameters
    - `κ0::Union{Real,ForwardDiff.Dual{Real}}` is the site kappa in units of s

# Examples
```julia-repl
    # define fas with
    FASParams( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, [ 1.0, 100.0, 140.0, 1000.0 ], [1.0, 0.0, 0.5, 0.5], 400.0, 0.0, 3.5, κ0 )
```
"""
struct FASParams
	# source parameters
	Δσ::Union{Real,ForwardDiff.Dual{Real}}             # stressParameter
	RΘϕ::Float64                                       # radiationPattern
	V::Float64                                         # partitionFactor
	F::Float64                                         # freeSurfaceFactor
	β::Float64                                         # sourceVelocity
  	ρ::Float64                                         # sourceDensity

	# path parameters
  	Rrefi::Vector{Float64}                             # geometricReferenceDistances
	γi::Vector{Union{Real,ForwardDiff.Dual{Real}}}     # geometricSpreadingRates
	Q0::Union{Real,ForwardDiff.Dual{Real}}             # qualityConstant
	η::Union{Real,ForwardDiff.Dual{Real}}              # qualityExponent
	cQ::Float64                                        # qualityVelocity

	# site parameters
	κ0::Union{Real,ForwardDiff.Dual{Real}}             # site kappa
end

# simplified constructor that takes just values of the Δσ and κ0
FASParams( Δσ, κ0 ) = FASParams( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, [ 1.0, 100.0, 140.0, 1000.0 ], [1.0, 0.0, 0.5, 0.5], 400.0, 0.0, 3.5, κ0 )
# simplified constructor that takes values of Δσ, Q0, η and κ0
FASParams( Δσ, Q0, η, κ0 ) = FASParams( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, [ 1.0, 100.0, 140.0, 1000.0 ], [1.0, 0.0, 0.5, 0.5], Q0, η, 3.5, κ0 )
# simplified constructor that takes values of Δσ, Q0 and κ0
FASParams( Δσ, Q0, κ0 ) = FASParams( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, [ 1.0, 100.0, 140.0, 1000.0 ], [1.0, 0.0, 0.5, 0.5], Q0, 0.0, 3.5, κ0 )
# simplified constructor that takes values of Δσ, Rrefi, γi, Q0 and κ0
FASParams( Δσ, Rrefi, γi, Q0, κ0 ) = FASParams( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, Q0, 0.0, 3.5, κ0 )
# simplified constructor that takes values of Δσ, Rrefi, γi, Q0, η and κ0
FASParams( Δσ, Rrefi, γi, Q0, η, κ0 ) = FASParams( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, Q0, η, 3.5, κ0 )


# create an alternative struct that allows for a more generic allocation of the geometric spreading rates
# for example, this allows the spreading rate for a particular segment to be constrained while others may
# adopt a ForwardDiff.Dual type
struct FASParamsGeo
	# source parameters
	Δσ::Union{Real,Dual{Real}}           # stressParameter
	RΘϕ::Float64                         # radiationPattern
	V::Float64                           # partitionFactor
	F::Float64                           # freeSurfaceFactor
	β::Float64                           # sourceVelocity
  	ρ::Float64                           # sourceDensity

	# path parameters
  	Rrefi::Vector{Float64}               # geometricReferenceDistances
	γi::Vector{Union{Real,Dual{Real}}}   # geometricSpreadingRates (to final segment)
	γf::Float64        	                 # fixed rate for last segment
	Q0::Union{Real,Dual{Real}}           # qualityConstant
	η::Union{Real,Dual{Real}}            # qualityExponent
	cQ::Float64                          # qualityVelocity

	# site parameters
	κ0::Union{Real,Dual{Real}}           # site kappa
end

# simplified constructor that takes values of Δσ, Rrefi, γi, Q0 and κ0
FASParamsGeo( Δσ, Rrefi, γi, γf, Q0, κ0 ) = FASParamsGeo( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, γf, Q0, 0.0, 3.5, κ0 )
# simplified constructor that takes values of Δσ, Rrefi, γi, Q0, η and κ0
FASParamsGeo( Δσ, Rrefi, γi, γf, Q0, η, κ0 ) = FASParamsGeo( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, γf, Q0, η, 3.5, κ0 )


struct FASParamsQr
	# source parameters
	Δσ::Union{Real,Dual{Real}}           # stressParameter
	RΘϕ::Float64                         # radiationPattern
	V::Float64                           # partitionFactor
	F::Float64                           # freeSurfaceFactor
	β::Float64                           # sourceVelocity
  	ρ::Float64                           # sourceDensity

	# path parameters
  	Rrefi::Vector{Float64}               # geometricReferenceDistances
	γi::Vector{Union{Real,Dual{Real}}}   # geometricSpreadingRates (to final segment)
	QRrefi::Vector{Float64}				 # qualityFactor reference distances
	Q0i::Vector{Union{Real,Dual{Real}}}  # qualityConstants
	ηi::Vector{Union{Real,Dual{Real}}}   # qualityExponent
	cQ::Float64                          # qualityVelocity

	# site parameters
	κ0::Union{Real,Dual{Real}}           # site kappa
end

# simplified constructor that takes values of Δσ, Rrefi, γi, Q0 and κ0
FASParamsQr( Δσ, Rrefi, γi, QRrefi, Q0i, κ0 ) = FASParamsQr( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, QRrefi, Q0i, zeros(length(Q0i)), 3.5, κ0 )
# simplified constructor that takes values of Δσ, Rrefi, γi, Q0, η and κ0
FASParamsQr( Δσ, Rrefi, γi, QRrefi, Q0i, ηi, κ0 ) = FASParamsQr( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, QRrefi, Q0i, ηi, 3.5, κ0 )


struct FASParamsGeoQr
	# source parameters
	Δσ::Union{Real,Dual{Real}}           # stressParameter
	RΘϕ::Float64                         # radiationPattern
	V::Float64                           # partitionFactor
	F::Float64                           # freeSurfaceFactor
	β::Float64                           # sourceVelocity
  	ρ::Float64                           # sourceDensity

	# path parameters
  	Rrefi::Vector{Float64}               # geometricReferenceDistances
	γi::Vector{Union{Real,Dual{Real}}}   # geometricSpreadingRates (to final segment)
	γf::Float64        	                 # fixed rate for last segment
	QRrefi::Vector{Float64}				 # qualityFactor reference distances
	Q0i::Vector{Union{Real,Dual{Real}}}  # qualityConstants
	ηi::Vector{Union{Real,Dual{Real}}}   # qualityExponent
	cQ::Float64                          # qualityVelocity

	# site parameters
	κ0::Union{Real,Dual{Real}}           # site kappa
end

# simplified constructor that takes values of Δσ, Rrefi, γi, Q0 and κ0
FASParamsGeoQr( Δσ, Rrefi, γi, γf, QRrefi, Q0i, κ0 ) = FASParamsGeoQr( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, γf, QRrefi, Q0i, zeros(length(Q0i)), 3.5, κ0 )
# simplified constructor that takes values of Δσ, Rrefi, γi, Q0, η and κ0
FASParamsGeoQr( Δσ, Rrefi, γi, γf, QRrefi, Q0i, ηi, κ0 ) = FASParamsGeoQr( Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, Rrefi, γi, γf, QRrefi, Q0i, ηi, 3.5, κ0 )




function fourier_constant(fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr})
	# note that the 1e-20 factor allows one to pass in distances in km, density in t/m^2, velocity in km/s
	return fas.RΘϕ * fas.V * fas.F * 1e-20 / ( 4π * fas.ρ * fas.β^3 * fas.Rrefi[1] )
end

# time_c(fas) = @time fourier_constant(fas)
# time_c(fas)
# time_c(gfas)

function fourier_source(f::T, m::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune) where T<:Real
	Mo = magnitude_to_moment(m)
	(fa, fb, ε) = corner_frequency(m, fas; fc_fun=fc_fun)
	if isnan(ε) == true
		# use single corner, with fa=fc
		return Mo / (1.0 + (f/fa)^2)
	else
		# use double corner model
		return Mo * ( (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2) )
	end
	# fc = corner_frequency(m, fas; fc_fun=fc_fun)
end

# time_Ef(f,m,fas) = @time fourier_source(f,m,fas)
# time_Ef(1.0,5.0,fas)
# time_Ef(1.0,5.0,gfas)

function fourier_source_shape(f::T, m::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune) where T<:Real
	# fc = corner_frequency(m, fas; fc_fun=fc_fun)
	# return 1.0 / (1.0 + (f/fc)^2)
	(fa, fb, ε) = corner_frequency(m, fas; fc_fun=fc_fun)
	if isnan(ε) == true
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	else
		# use double corner model
		return (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2)
	end
end

function fourier_source_shape(f::T, fa::Real, fb::Real, ε::Real) where T<:Real
	if isnan(ε) == true
		# use single corner, with fa=fc
		return 1.0 / (1.0 + (f/fa)^2)
	else
		# use double corner model
		return (1.0 - ε)/(1.0 + (f/fa)^2 ) + ε/(1.0 + (f/fb)^2)
	end
end



function fourier_path(f::T, r::T, fas::Union{FASParams,FASParamsGeo}) where T<:Real
	Zr = geometric_spreading( r, fas )
	return Zr * exp( -π*f*r / ( fas.Q0 * f^fas.η * fas.cQ ) )
end

function fourier_path(f::T, r::T, fas::Union{FASParamsQr,FASParamsGeoQr}) where T<:Real
	Zr = geometric_spreading( r, fas )
	Qr = 1.0
	if r <= fas.QRrefi[1]
		Qr *= exp( -π*f*r / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
	else
		Qr *= exp( -π*f*fas.QRrefi[1] / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
		for i in 2:length(fas.QRrefi)
			if r <= fas.QRrefi[i]
				Qr *= exp( -π*f*(r-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
			else
				Qr *= exp( -π*f*(fas.QRrefi[i]-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
			end
		end
	end
	return Zr * Qr
end

function fourier_path_cy(f::T, r::T, fas::Union{FASParams,FASParamsGeo}) where T<:Real
	Zr = geometric_spreading_cy( r, fas )
	return Zr * exp( -π*f*r / ( fas.Q0 * f^fas.η * fas.cQ ) )
end

function fourier_path_cy(f::T, r::T, fas::Union{FASParamsQr,FASParamsGeoQr}) where T<:Real
	Zr = geometric_spreading_cy( r, fas )
	Qr = 1.0
	if r <= fas.QRrefi[1]
		Qr *= exp( -π*f*r / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
	else
		Qr *= exp( -π*f*fas.QRrefi[1] / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
		for i in 2:length(fas.QRrefi)
			if r <= fas.QRrefi[i]
				Qr *= exp( -π*f*(r-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
			else
				Qr *= exp( -π*f*(fas.QRrefi[i]-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
			end
		end
	end
	return Zr * Qr
end


# time_Pf(f,r,fas) = @time fourier_path(f,r,fas)
# time_Pf(1.0,100.0,fas)
# time_Pf(1.0,100.0,gfas)

function fourier_attenuation(f::T, r::T, fas::Union{FASParams,FASParamsGeo}) where T<:Real
	if f < eps()
		return 1.0
	else
		return exp( -π*f*r / ( fas.Q0 * f^fas.η * fas.cQ ) - π*f*fas.κ0 )
	end
end

function fourier_attenuation(f::Real, r::T, fas::Union{FASParams,FASParamsGeo}) where T<:Real
	if f < eps()
		return 1.0
	else
		return exp( -π*f*r / ( fas.Q0 * f^fas.η * fas.cQ ) - π*f*fas.κ0 )
	end
end

function fourier_attenuation(f::T, r::T, fas::Union{FASParamsQr,FASParamsGeoQr}) where T<:Real
	if f < eps()
		return 1.0
	else
		Kr = exp( -π*f*fas.κ0 )
		Qr = 1.0
		if r <= fas.QRrefi[1]
			Qr *= exp( -π*f*r / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
		else
			Qr *= exp( -π*f*fas.QRrefi[1] / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
			for i in 2:length(fas.QRrefi)
				if r <= fas.QRrefi[i]
					Qr *= exp( -π*f*(r-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
				else
					Qr *= exp( -π*f*(fas.QRrefi[i]-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
				end
			end
		end
		return Kr * Qr
	end
end

function fourier_attenuation(f::Real, r::T, fas::Union{FASParamsQr,FASParamsGeoQr}) where T<:Real
	if f < eps()
		return 1.0
	else
		Kr = exp( -π*f*fas.κ0 )
		Qr = 1.0
		if r <= fas.QRrefi[1]
			Qr *= exp( -π*f*r / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
		else
			Qr *= exp( -π*f*fas.QRrefi[1] / ( fas.Q0i[1] * f^fas.ηi[1] * fas.cQ ) )
			for i in 2:length(fas.QRrefi)
				if r <= fas.QRrefi[i]
					Qr *= exp( -π*f*(r-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
				else
					Qr *= exp( -π*f*(fas.QRrefi[i]-fas.QRrefi[i-1]) / ( fas.Q0i[i] * f^fas.ηi[i] * fas.cQ ) )
				end
			end
		end
		return Kr * Qr
	end
end

# time_κ(f,r,fas) = @time fourier_attenuation(f,r,fas)
# time_κ(1.0,100.0,fas)
# time_κ(1.0,100.0,gfas)

function fourier_site(f::Real, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; amp_model::Symbol=:AlAtik2021_cy14)
	# kappa filter
    Sf = exp( -π*f*fas.κ0 )
    # site impedence function
    Sif = site_amplification(f; amp_model=amp_model )
	return Sf*Sif
end

# time_Sf(f, fas) = @time fourier_site(f, fas)
# time_Sf(1.0,fas)
# time_Sf(1.0,gfas)

function fourier_spectral_ordinate( f::T, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	# define all constant terms here
  	C = fourier_constant(fas)
	# source term
	Ef = fourier_source(f, m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = finite_fault_factor(m)
	r_ps = sqrt(r^2 + h^2)
	# geometric spreading
	Gr = geometric_spreading(r_ps, fas)
	# combined attenuation
	Kf = fourier_attenuation(f, r_ps, fas)
	# site impedance
    Sf = site_amplification(f; amp_model=amp_model)

	# overall displacement spectrum (in cms) is Ef*Pf*Sf
	# the division by 10^2 below is to convert from units of cm/s to m/s
	return (2π*f)^2 * C * Ef * Gr * Kf * Sf / 100.0
end

# time_fas(f,m,r,fas) = @time fourier_spectral_ordinate(f, m, r, fas)
# time_fas(1.0,5.0,60.0,fas)
# time_fas(1.0,5.0,60.0,gfas)


function fourier_spectrum( f::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	# define all frequency independent terms here
  	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = finite_fault_factor(m)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading(r_ps, fas)
	factor = 4π^2 * C * Mo * Gr / 100.0

	Af = ones(typeof(fas.Δσ), size(f))
  	for i = 1:length(f)
	  	# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i]; amp_model=amp_model)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af[i] = Ef * Kf * Sf * factor * f[i]^2
    end
    return Af
end



function fourier_spectrum!( Af::Vector, f::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	# define all frequency independent terms here
	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = finite_fault_factor(m)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading(r_ps, fas)
	factor = 4π^2 * C * Mo * Gr / 100.0

	for i = 1:length(f)
		# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i]; amp_model=amp_model)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af[i] = Ef * Kf * Sf * factor * f[i]^2
  	end
  	nothing
end

# @btime fourier_spectrum(ff,5.0,10.0,fas)
# ff = collect(range(0.01, stop=10.0, step=0.01))
# time_fs(f,m,r,fas) = @time fourier_spectrum(f,m,r,fas)
# time_fs(ff,5.0,10.0,fas)
# time_fs(ff,5.0,10.0,gfas)


function squared_fourier_spectrum!( Af2::Vector, f::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	# define all frequency independent terms here
	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = finite_fault_factor(m)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading(r_ps, fas)
	factor = (4π^2 * C * Mo * Gr / 100.0)^2

	for i = 1:length(f)
		# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i]; amp_model=amp_model)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af2[i] = (Ef * Kf * Sf)^2 * factor * f[i]^4
  	end
  	nothing
end

function squared_fourier_spectrum!( Af2::Vector, f::Vector, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	# define all frequency independent terms here
	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = finite_fault_factor(m)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading(r_ps, fas)
	factor = (4π^2 * C * Mo * Gr / 100.0)^2

	for i = 1:length(f)
		# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i]; amp_model=amp_model)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af2[i] = (Ef * Kf * Sf)^2 * factor * f[i]^4
  	end
  	nothing
end


function squared_fourier_spectrum_cy!( Af2::Vector, f::Vector{T}, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	# define all frequency independent terms here
	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = finite_fault_factor(m)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading_cy(r_ps, fas)
	factor = (4π^2 * C * Mo * Gr / 100.0)^2

	for i = 1:length(f)
		# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i]; amp_model=amp_model)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af2[i] = (Ef * Kf * Sf)^2 * factor * f[i]^4
  	end
  	nothing
end


function squared_fourier_spectrum_cy!( Af2::Vector, f::Vector, m::T, r::T, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune, amp_model::Symbol=:AlAtik2021_cy14 ) where T<:Real
	# define all frequency independent terms here
	C = fourier_constant(fas)
	Mo = magnitude_to_moment(m)
	fa, fb, ε = corner_frequency(m, fas; fc_fun=fc_fun)
	# use point-source distance metric
	h = finite_fault_factor(m)
	r_ps = sqrt(r^2 + h^2)
	Gr = geometric_spreading_cy(r_ps, fas)
	factor = (4π^2 * C * Mo * Gr / 100.0)^2

	for i = 1:length(f)
		# source term
		@inbounds Ef = fourier_source_shape(f[i], fa, fb, ε)
  		# combined attenuation
  		@inbounds Kf = fourier_attenuation(f[i], r_ps, fas)
  		# site impedance
      	@inbounds Sf = site_amplification(f[i]; amp_model=amp_model)
		# apply factor and convert to acceleration in appropriate units
  		@inbounds Af2[i] = (Ef * Kf * Sf)^2 * factor * f[i]^4
  	end
  	nothing
end



# function gives the frequency at which the combined κ_r and κ_0 filters (squared versions) give a value of 0.5
function combined_kappa_frequency( r::Real, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr} )
  if fas.η < 0.1
    # a closed form solution exists
    return log(2.0) / ( 2π * ( r / (fas.Q0 * fas.cQ) + fas.κ0 ) )
  else
	g(f) = 0.5 - (fourier_attenuation(f, r, fas)^2)
    # g(f) = 0.5 - exp(-2π*f*(r/(fas.Q0*f^fas.η*fas.cQ) + fas.κ0))
	f_0 = find_zero( g, (0.01,100.0), Bisection(); xatol=1e-2 )
	fk = min(max(f_0, 0.2), 1.0)
    return fk
  end
end
