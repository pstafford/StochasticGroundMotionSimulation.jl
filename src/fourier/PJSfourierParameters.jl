
"""
	SourceParameters

Custom type defining the source parameters of a Fourier spectrum

- `Δσ::Union{Float64,Dual{Float64}}` is the stress parameter in bars
- `RΘϕ::Float64` is the radiation pattern
- `V::Float64` is the partition factor (for splitting to horizontal components)
- `F::Float64` is the free surface factor
- `β::Float64` is the source velocity in units of km/s
- `ρ::Float64` is the source density in units of t/m^3 or g/cm^3
- `src_model::Symbol` identifies the type of source spectrum (:Brune, Atkinson_Silva_2000)
"""
struct SourceParameters
	# source parameters
	Δσ::Union{Float64,Dual{Float64}}	# stressParameter
	RΘϕ::Float64                            # radiationPattern
	V::Float64                              # partitionFactor
	F::Float64                              # freeSurfaceFactor
	β::Float64                              # sourceVelocity
  	ρ::Float64                              # sourceDensity
	src_model::Symbol 						# source spectrum model
end

SourceParameters(Δσ) = SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, :Brune)
SourceParameters(Δσ, βs, ρs) = SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, βs, ρs, :Brune)



"""
	PathParameters

Custom type defining the path parameters of a Fourier spectrum

- `Rrefi::Vector{Float64}` is a vector of reference distances for geometric spreading in units of km
- `γi::Vector{Union{Float64,Dual{Float64}}` is a vector of geometric spreading rates. Note that to initialize individual elements we need to pass `Union{Float64,Dual{Float64}}[ 1.0, Dual{Float64}(0.5) ]`, for example
- `Q0::Union{Float64,Dual{Float64}}` is the Quality Constant
- `η::Union{Float64,Dual{Float64}}` is the Quality Exponent
- `cQ::Float64` is the velocity used to determine Q(f) (or the travel time) in units of km/s
- `geo_model::Symbol` indicates the type of geometric spreading `:Piecewise` or `:CY14`
- `sat_model::Symbol` indicates the near source saturation model `:BT15`
"""
struct PathParameters
	# path parameters
  	Rrefi::Vector{Float64} 	# geometricReferenceDistances
	γi::Vector{Union{Float64,Dual{Float64,Float64}}}        # geometricSpreadingRates
	Q0::Union{Float64,Dual{Float64}} # qualityConstant
	η::Union{Float64,Dual{Float64}}  # qualityExponent
	cQ::Float64             # qualityVelocity
	geo_model::Symbol 		# geometric spreading model
	sat_model::Symbol		# saturation model
	heff::Union{Float64,Dual{Float64}} # near-source saturation
end

PathParameters(Rrefi, γi, Q0) = PathParameters(Rrefi, γi, Q0, 0.0, 3.5, :Piecewise, :BT15, NaN)
PathParameters(Rrefi, γi, Q0, η) = PathParameters(Rrefi, γi, Q0, η, 3.5, :Piecewise, :BT15, NaN)


"""
	SiteParameters

Custom type defining the site parameters of a Fourier spectrum

- `κ0::Union{Float64,Dual{Float64}}` is the site kappa in units of s
- `amp_model::Symbol` is a symbol identifying the impedance function
"""
struct SiteParameters
    # site parameters
	κ0::Union{Float64,Dual{Float64}}            # site kappa
	amp_model::Symbol	# site amplification model
end

SiteParameters(κ0) = SiteParameters(κ0, :AlAtik2021_cy14)



"""
    FourierParameters

Custom type for the parameters the Fourier amplitude spectrum. This type is comprised of source, path and site types.
"""
struct FourierParameters
	source::SourceParameters	# source parameters
	path::PathParameters		# path parameters
	site::SiteParameters 		# site parameters
end

# simplified constructor that takes just values of the Δσ and κ0
FourierParameters(Δσ, κ0) = FourierParameters(SourceParameters(Δσ), PathParameters([1.0, 50.0, Inf], [1.0, 0.5], 150.0, 0.4), SiteParameters(κ0))

# simplified constructor that takes values of Δσ, Q0, η and κ0
FourierParameters(Δσ, Q0, η, κ0) = FourierParameters(SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, :Brune), PathParameters([1.0, 50.0, Inf], [1.0, 0.5], Q0, η, 3.5, :Piecewise, :BT15), SiteParameters(κ0, :AlAtik2021_cy14))
# simplified constructor that takes values of Δσ, Q0 and κ0
FourierParameters(Δσ, Q0, κ0) = FourierParameters(SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, :Brune), PathParameters([1.0, 50.0, Inf], [1.0, 0.5], Q0, 0.0, 3.5, :Piecewise, :BT15), SiteParameters(κ0, :AlAtik2021_cy14))
# simplified constructor that takes values of Δσ, Rrefi, γi, Q0 and κ0
FourierParameters(Δσ, Rrefi, γi, Q0, κ0) = FourierParameters(SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, :Brune), PathParameters(Rrefi, γi, Q0, 0.0, 3.5, :Piecewise, :BT15), SiteParameters(κ0, :AlAtik2021_cy14))
# simplified constructor that takes values of Δσ, Rrefi, γi, Q0, η and κ0
FourierParameters(Δσ, Rrefi, γi, Q0, η, κ0) = FourierParameters(SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, :Brune), PathParameters(Rrefi, γi, Q0, η, 3.5, :Piecewise, :BT15), SiteParameters(κ0, :AlAtik2021_cy14))
