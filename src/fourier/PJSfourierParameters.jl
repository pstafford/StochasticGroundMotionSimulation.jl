
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
struct SourceParameters{S<:Float64, T<:Real}
	# source parameters
	Δσ::T	# stressParameter
	RΘϕ::S                            # radiationPattern
	V::S                              # partitionFactor
	F::S                              # freeSurfaceFactor
	β::S                              # sourceVelocity
  	ρ::S                              # sourceDensity
	model::Symbol  					  # source spectrum model
end

SourceParameters(Δσ::T) where {T<:Float64} = SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, :Brune)
SourceParameters(Δσ::T) where {T<:Real} = SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, :Brune)
SourceParameters(Δσ::T, model::Symbol) where {T<:Real} = SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, 3.5, 2.75, model)
SourceParameters{S, T}(Δσ::T, βs::S, ρs::S) where {S<:Float64,T<:Real} = SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, βs, ρs, :Brune)
SourceParameters{T, T}(Δσ::T, βs::T, ρs::T) where {T<:Float64} = SourceParameters(Δσ, 0.55, 1.0/sqrt(2.0), 2.0, βs, ρs, :Brune)

"""
	get_parametric_type(src::SourceParameters{S,T}) where {S,T} = T

Extract type of `T` from parametric `SourceParameters` struct
"""
get_parametric_type(src::SourceParameters{S,T}) where {S,T} = T


"""
	GeometricSpreadingParameters{S<:Real,T<:Real,U<:AbstractVector{Bool}}

Struct for geometric spreading parameters.
Holds fields:
- `Rrefi` are reference distances, these are `<:Real` but will generally be `Float64` values
- `γconi` are constant spreading rates, meaning that they will not be free for AD purposes
- `γvari` are variable spreading rates, meaning that they can be represented as `Dual` numbers for AD
- `γfree` is a vector of `Bool` instances, or a `BitVector` that indicates which segments are constant or variable. Variable spreading rates are given `1` or `true`
- `model` is a symbol defining the type of spreading model `:Piecewise` or `:CY14`
"""
struct GeometricSpreadingParameters{S<:Real,T<:Real,U<:AbstractVector{Bool}}
	Rrefi::Vector{S}
	γconi::Vector{S}
	γvari::Vector{T}
	γfree::U
	model::Symbol
end

GeometricSpreadingParameters(Rrefi, γconi) = GeometricSpreadingParameters(Rrefi, γconi, Vector{Float64}(), BitVector(undef,length(γconi)), :Piecewise )
GeometricSpreadingParameters(Rrefi, γconi, model) = GeometricSpreadingParameters(Rrefi, γconi, Vector{Float64}(), BitVector(undef,length(γconi)), model )
GeometricSpreadingParameters(Rrefi::Vector{T}, γconi::Vector{T}) where T = GeometricSpreadingParameters{T,T,BitVector}(Rrefi, γconi, Vector{T}(), BitVector(undef,length(γconi)), :Piecewise )
GeometricSpreadingParameters(Rrefi::Vector{S}, γvari::Vector{T}) where {S<:Float64,T<:Dual} = GeometricSpreadingParameters{S,T,BitVector}(Rrefi, Vector{S}(), γvari, BitVector(ones(length(γvari))), :Piecewise )

"""
	get_parametric_type(geo::GeometricSpreadingParameters{S,T,U}) where {S,T,U} = T

Extract type of `T` from parametric `GeometricSpreadingParameters` struct
"""
get_parametric_type(geo::GeometricSpreadingParameters{S,T,U}) where {S,T,U} = T


"""
	NearSourceSaturationParameters{S<:Real,T<:Real,U<:AbstractVector{Bool}}

Struct for near-source saturation parameters. Mimic structure of the `GeometricSpreadingParameters` struct.
Holds fields:
- `mRefi` reference magnitudes
- `hconi` constrained coefficients, not free for AD purposes
- `hvari` variable coefficeints, free for AD purposes
- `hfree` is a vector of `Bool` instances, or a `BitVector` indicating which parameters are constant or variable
- `model` is a symbol defining the type of saturation model:
	- `:BT15` is Boore & Thompson (2015)
	- `:YA15` is Yenier & Atkinson (2015)
	- `:None` returns zero saturation length
	- `:ConstantConstrained` is a fixed saturation length not subject to AD operations
	- `:ConstantVariable` is a fixed saturation length that is subject to AD operations (i.e., is a `<:Dual`)
"""
struct NearSourceSaturationParameters{S<:Real,T<:Real,U<:AbstractVector{Bool}}
	mRefi::Vector{S}
	hconi::Vector{S}
	hvari::Vector{T}
	hfree::U
	model::Symbol
end

NearSourceSaturationParameters(model::Symbol) = NearSourceSaturationParameters(Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), BitVector(), model)
NearSourceSaturationParameters{T,T,U}(mRefi::Vector{T}, hconi::Vector{T}, model::Symbol) where {T,U} = NearSourceSaturationParameters(mRefi, hconi, Vector{T}(), BitVector(undef, length(hconi)), model)
NearSourceSaturationParameters{S,T,U}(mRefi::Vector{S}, hvari::Vector{T}, model::Symbol) where {S,T,U} = NearSourceSaturationParameters(mRefi, Vector{S}(), hvari, BitVector(ones(length(hvari))), model)

NearSourceSaturationParameters{T,T,U}(mRefi::Vector{T}, hconi::Vector{T}) where {T,U} = NearSourceSaturationParameters(mRefi, hconi, Vector{T}(), BitVector(undef, length(hconi)), :FullyConstrained)
NearSourceSaturationParameters{S,T,U}(mRefi::Vector{S}, hvari::Vector{T}) where {S,T,U} = NearSourceSaturationParameters(mRefi, Vector{S}(), hvari, BitVector(ones(length(hvari))), :FullyVariable)

# specialisation for a constant saturation term
NearSourceSaturationParameters(hcon::Float64) = NearSourceSaturationParameters(Vector{Float64}(), [ hcon ], Vector{Float64}(), BitVector(undef,1), :ConstantConstrained)
NearSourceSaturationParameters(hvar::T) where {T<:Dual} = NearSourceSaturationParameters(Vector{Float64}(), Vector{Float64}(), Vector{T}([ hvar ]), BitVector([1]), :ConstantVariable)

"""
	get_parametric_type(sat::NearSourceSaturationParameters{S,T,U}) where {S,T,U} = T

Extract type of `T` from parametric `NearSourceSaturationParameters` struct
"""
get_parametric_type(sat::NearSourceSaturationParameters{S,T,U}) where {S,T,U} = T


"""
	AnelasticAttenuationParameters{S<:Real,T<:Real}

Struct for anelastic attenuation parameters.
Holds fields:
- `Q0` quality factor at 1 Hz
- `η` quality exponent ∈ [0,1)
- `cQ` velocity (km/s) along propagation path used to determine `Q(f)`
"""
struct AnelasticAttenuationParameters{S<:Real,T<:Real}
	Q0::S
	η::T
	cQ::Float64
end

AnelasticAttenuationParameters(Q0) = AnelasticAttenuationParameters(Q0, 0.0, 3.5)
AnelasticAttenuationParameters(Q0, η) = AnelasticAttenuationParameters(Q0, η, 3.5)
AnelasticAttenuationParameters(Q0::T, η::T) where T = AnelasticAttenuationParameters{T,T}(Q0, η, 3.5)
AnelasticAttenuationParameters(Q0::T, η::T, cQ::Float64) where T = AnelasticAttenuationParameters{T,T}(Q0, η, cQ)

"""
	get_parametric_type(anelastic::AnelasticAttenuationParameters{S,T}) where {S,T}

Extract type most elaborate type `S` or `T` from parametric `AnelasticAttenuationParameters` struct
"""
function get_parametric_type(anelastic::AnelasticAttenuationParameters{S,T}) where {S,T}
	if S <: Float64
		if T <: Float64
			return S
		else
			return T
		end
	else
		return S
	end
end



"""
	PathParameters

Custom type defining the path parameters of a Fourier spectrum.
Consists of three other custom structs
- `geometric` is a `GeometricSpreadingParameters` type
- `saturation` is a `NearSourceSaturationParameters` type
- `anelastic` is an `AnelasticAttenuationParameters` type
"""
struct PathParameters{G<:GeometricSpreadingParameters,S<:NearSourceSaturationParameters,A<:AnelasticAttenuationParameters}
	geometric::G
	saturation::S
	anelastic::A
end

PathParameters(geometric::GeometricSpreadingParameters, anelastic::AnelasticAttenuationParameters) = PathParameters(geometric, NearSourceSaturationParameters(:None), anelastic)

"""
	get_parametric_type(path::PathParameters)

Extract type most elaborate type from parametric `PathParameters` struct.
This requires dropping down to lower level structs within `path`.
"""
function get_parametric_type(path::PathParameters)
	T = get_parametric_type(path.geometric)
	U = get_parametric_type(path.saturation)
	V = get_parametric_type(path.anelastic)
	if T <: Float64
		if U <: Float64
			return V
		else
			return U
		end
	else
		return T
	end
end



"""
	SiteParameters

Custom type defining the site parameters of a Fourier spectrum

- `κ0::Union{Float64,Dual{Float64}}` is the site kappa in units of s
- `amp_model::Symbol` is a symbol identifying the impedance function
"""
struct SiteParameters{T<:Real}
    # site parameters
	κ0::T            	# site kappa
	model::Symbol		# site amplification model
end

SiteParameters(κ0::T) where {T} = SiteParameters(κ0, :AlAtik2021_cy14)

"""
	get_parametric_type(site::SiteParameters{T}) where {T} = T

Extract type of `T` from parametric `SiteParameters` struct
"""
get_parametric_type(site::SiteParameters{T}) where {T} = T


"""
    FourierParameters

Custom type for the parameters the Fourier amplitude spectrum. This type is comprised of source, path and site types.
"""
struct FourierParameters{S<:SourceParameters, T<:PathParameters, U<:SiteParameters}
	source::S	# source parameters
	path::T		# path parameters
	site::U		# site parameters
end

# initialiser focussing upon source and path response only. Uses zero kappa and unit amplification
FourierParameters(src::S, path::T) where {S<:SourceParameters, T<:PathParameters} = FourierParameters(src, path, SiteParameters(0.0, :Unit))

"""
	get_parametric_type(fas::FourierParameters)

Extract type most elaborate type from parametric `FourierParameters` struct.
This requires dropping down to lower level structs within `fas`.
"""
function get_parametric_type(fas::FourierParameters)
	T = get_parametric_type(fas.source)
	U = get_parametric_type(fas.path)
	V = get_parametric_type(fas.site)
	if T <: Float64
		if U <: Float64
			return V
		else
			return U
		end
	else
		return T
	end
end

# simplified constructor that takes just values of the Δσ and κ0
# FourierParameters(Δσ, κ0) = FourierParameters(SourceParameters(Δσ), PathParameters([1.0, 50.0, Inf], [1.0, 0.5], 150.0, 0.4), SiteParameters(κ0))
# # simplified constructor that takes values of Δσ, Q0 and κ0
# FourierParameters(Δσ, Q0, κ0) = FourierParameters(SourceParameters(Δσ), PathParameters([1.0, 50.0, Inf], [1.0, 0.5], Q0), SiteParameters(κ0))
# # simplified constructor that takes values of Δσ, Q0, η and κ0
# FourierParameters(Δσ, Q0, η, κ0) = FourierParameters(SourceParameters(Δσ), PathParameters([1.0, 50.0, Inf], [1.0, 0.5], Q0, η), SiteParameters(κ0))
# # simplified constructor that takes values of Δσ, Rrefi, γi, Q0 and κ0
# FourierParameters(Δσ, Rrefi, γi, Q0, κ0) = FourierParameters(SourceParameters(Δσ), PathParameters(Rrefi, γi, Q0), SiteParameters(κ0))
# # simplified constructor that takes values of Δσ, Rrefi, γi, Q0, η and κ0
# FourierParameters(Δσ, Rrefi, γi, Q0, η, κ0) = FourierParameters(SourceParameters(Δσ), PathParameters(Rrefi, γi, Q0, η), SiteParameters(κ0))
