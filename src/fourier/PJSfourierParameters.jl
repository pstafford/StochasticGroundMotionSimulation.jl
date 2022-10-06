
"""
	SourceParameters

Custom type defining the source parameters of a Fourier spectrum.

Constructed with signature `SourceParameters{S<:Float64, T<:Real}` with fields:
- `Δσ::T` is the stress parameter in bars
- `RΘϕ::S` is the radiation pattern
- `V::S` is the partition factor (for splitting to horizontal components)
- `F::S` is the free surface factor
- `β::S` is the source velocity in units of km/s
- `ρ::S` is the source density in units of t/m³ or g/cm³
- `model::Symbol` identifies the type of source spectrum (`:Brune`, `:Atkinson_Silva_2000`)
"""
struct SourceParameters{S<:Float64,T<:Real,U<:Real}
    # source parameters
    Δσ::T# stressParameter
    RΘϕ::S                            # radiationPattern
    V::S                              # partitionFactor
    F::S                              # freeSurfaceFactor
    β::S                              # sourceVelocity
    ρ::S                              # sourceDensity
    n::U                              # high-frequency fall-off rate for Beresnev (2019)
    model::Symbol    # source spectrum model
end

SourceParameters(Δσ::T) where {T<:Float64} = SourceParameters(Δσ, 0.55, 1.0 / sqrt(2.0), 2.0, 3.5, 2.75, 1.0, :Brune)
SourceParameters(Δσ::T) where {T<:Real} = SourceParameters(Δσ, 0.55, 1.0 / sqrt(2.0), 2.0, 3.5, 2.75, 1.0, :Brune)
SourceParameters(Δσ::T, model::Symbol) where {T<:Real} = SourceParameters(Δσ, 0.55, 1.0 / sqrt(2.0), 2.0, 3.5, 2.75, 1.0, model)
SourceParameters(Δσ::T, βs::S, ρs::S) where {S<:Float64,T<:Real} = SourceParameters(Δσ, 0.55, 1.0 / sqrt(2.0), 2.0, βs, ρs, 1.0, :Brune)
SourceParameters(Δσ::T, βs::T, ρs::T) where {T<:Float64} = SourceParameters(Δσ, 0.55, 1.0 / sqrt(2.0), 2.0, βs, ρs, 1.0, :Brune)
SourceParameters(Δσ::T, n::T) where {T<:Real} = SourceParameters(Δσ, 0.55, 1.0 / sqrt(2.0), 2.0, 3.5, 2.75, n, :Beresnev_2019)
SourceParameters(Δσ::T, n::U) where {T<:Real,U<:Real} = SourceParameters(Δσ, 0.55, 1.0 / sqrt(2.0), 2.0, 3.5, 2.75, n, :Beresnev_2019)


"""
	get_parametric_type(src::SourceParameters{S,T}) where {S,T} = T

Extract type of `T` from parametric `SourceParameters` struct
"""
get_parametric_type(src::SourceParameters{S,T}) where {S,T} = T


"""
	GeometricSpreadingParameters

Struct for geometric spreading parameters.
Holds fields:
- `Rrefi` are reference distances, these are `<:Real` but will generally be `Float64` values
- `γconi` are constant spreading rates, meaning that they will not be free for AD purposes
- `γvari` are variable spreading rates, meaning that they can be represented as `Dual` numbers for AD
- `γfree` is a vector of `Bool` instances, or a `BitVector` that indicates which segments are constant or variable. Variable spreading rates are given `1` or `true`
- `model` is a symbol defining the type of spreading model `:Piecewise`, `:CY14`, `:CY14mod`
"""
struct GeometricSpreadingParameters{S<:Real,T<:Real,U<:AbstractVector{Bool}}
    Rrefi::Vector{S}
    γconi::Vector{S}
    γvari::Vector{T}
    γfree::U
    model::Symbol
end

GeometricSpreadingParameters(Rrefi::Vector{T}, γconi::Vector{T}) where {T} = GeometricSpreadingParameters{T,T,BitVector}(Rrefi, γconi, Vector{T}(), BitVector(zeros(length(γconi))), :Piecewise)
GeometricSpreadingParameters(Rrefi::Vector{T}, γconi::Vector{T}, model::Symbol) where {T} = GeometricSpreadingParameters{T,T,BitVector}(Rrefi, γconi, Vector{T}(), BitVector(zeros(length(γconi))), model)
GeometricSpreadingParameters(Rrefi::Vector{S}, γvari::Vector{T}) where {S<:Float64,T<:Dual} = GeometricSpreadingParameters{S,T,BitVector}(Rrefi, Vector{S}(), γvari, BitVector(ones(length(γvari))), :Piecewise)
GeometricSpreadingParameters(Rrefi::Vector{S}, γvari::Vector{T}, model::Symbol) where {S<:Float64,T<:Dual} = GeometricSpreadingParameters{S,T,BitVector}(Rrefi, Vector{S}(), γvari, BitVector(ones(length(γvari))), model)

"""
	get_parametric_type(geo::GeometricSpreadingParameters{S,T,U}) where {S,T,U} = T

Extract type of `T` from parametric `GeometricSpreadingParameters` struct
"""
get_parametric_type(geo::GeometricSpreadingParameters{S,T,U}) where {S,T,U} = T


"""
	NearSourceSaturationParameters

Struct for near-source saturation parameters. Mimic structure of the `GeometricSpreadingParameters` struct.
Holds fields:
- `mRefi` reference magnitudes
- `hconi` constrained coefficients, not free for AD purposes
- `hvari` variable coefficients, free for AD purposes
- `hfree` is a vector of `Bool` instances, or a `BitVector` indicating which parameters are constant or variable
- `exponent` is the exponent used within equivalent point-source distance calculations: ``r_{ps} = \\left[r_{rup}^n + h(m)^n\\right]^{1/n}``
- `model` is a symbol defining the type of saturation model:
  - `:BT15` is Boore & Thompson (2015)
  - `:YA15` is Yenier & Atkinson (2015)
  - `:CY14` is average Chiou & Youngs (2014)
  - `:None` returns zero saturation length
  - `:ConstantConstrained` is a fixed saturation length not subject to AD operations
  - `:ConstantVariable` is a fixed saturation length that is subject to AD operations (i.e., is a `<:Dual`)

"""
struct NearSourceSaturationParameters{S<:Real,T<:Real,U<:AbstractVector{Bool}}
    mRefi::Vector{S}
    hconi::Vector{S}
    hvari::Vector{T}
    hfree::U
    exponent::Int
    model::Symbol
end

NearSourceSaturationParameters(model::Symbol) = NearSourceSaturationParameters(Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), BitVector(), 2, model)
NearSourceSaturationParameters(exponent::Int, model::Symbol) = NearSourceSaturationParameters(Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), BitVector(), exponent, model)
NearSourceSaturationParameters(mRefi::Vector{T}, hconi::Vector{T}, model::Symbol) where {T} = NearSourceSaturationParameters(mRefi, hconi, Vector{T}(), BitVector(undef, length(hconi)), 2, model)
NearSourceSaturationParameters(mRefi::Vector{S}, hvari::Vector{T}, model::Symbol) where {S,T} = NearSourceSaturationParameters(mRefi, Vector{S}(), hvari, BitVector(ones(length(hvari))), 2, model)

NearSourceSaturationParameters(mRefi::Vector{T}, hconi::Vector{T}) where {T} = NearSourceSaturationParameters(mRefi, hconi, Vector{T}(), BitVector(undef, length(hconi)), 2, :FullyConstrained)
NearSourceSaturationParameters(mRefi::Vector{S}, hvari::Vector{T}) where {S,T} = NearSourceSaturationParameters(mRefi, Vector{S}(), hvari, BitVector(ones(length(hvari))), 2, :FullyVariable)

# specialisation for a constant saturation term
NearSourceSaturationParameters(hcon::Float64) = NearSourceSaturationParameters(Vector{Float64}(), [hcon], Vector{Float64}(), BitVector(undef, 1), 2, :ConstantConstrained)
NearSourceSaturationParameters(hvar::T) where {T<:Dual} = NearSourceSaturationParameters(Vector{Float64}(), Vector{Float64}(), Vector{T}([hvar]), BitVector([1]), 2, :ConstantVariable)
# with exponents
NearSourceSaturationParameters(hcon::Float64, exponent::Int) = NearSourceSaturationParameters(Vector{Float64}(), [hcon], Vector{Float64}(), BitVector(undef, 1), exponent, :ConstantConstrained)
NearSourceSaturationParameters(hvar::T, exponent::Int) where {T<:Dual} = NearSourceSaturationParameters(Vector{Float64}(), Vector{Float64}(), Vector{T}([hvar]), BitVector([1]), exponent, :ConstantVariable)



"""
	get_parametric_type(sat::NearSourceSaturationParameters{S,T,U}) where {S,T,U} = T

Extract type of `T` from parametric `NearSourceSaturationParameters` struct
"""
get_parametric_type(sat::NearSourceSaturationParameters{S,T,U}) where {S,T,U} = T


"""
	AnelasticAttenuationParameters

Struct for anelastic attenuation parameters.
Holds fields:
- `Q0` quality factor at 1 Hz
- `η` quality exponent ∈ [0,1)
- `cQ` velocity (km/s) along propagation path used to determine `Q(f)`
- `rmetric` is a symbol `:Rrup` or `:Rps` to define which distance metric is used for anelastic attenuation
"""
struct AnelasticAttenuationParameters{S<:Real,T<:Real,U<:AbstractVector{Bool}}
    Rrefi::Vector{S}
    Q0coni::Vector{S}
    Q0vari::Vector{T}
    ηconi::Vector{S}
    ηvari::Vector{T}
    cQ::Vector{Float64}
    Qfree::U
    ηfree::U
    rmetric::Symbol
end

AnelasticAttenuationParameters(Q0::T) where {T<:Float64} = AnelasticAttenuationParameters([0.0, Inf], [Q0], Vector{Float64}(), [0.0], Vector{Float64}(), [3.5], BitVector(zeros(1)), BitVector(zeros(1)), :Rps)
AnelasticAttenuationParameters(Q0::T) where {T<:Real} = AnelasticAttenuationParameters([0.0, Inf], Vector{Float64}(), [Q0], [0.0], Vector{T}(), [3.5], BitVector(ones(1)), BitVector(zeros(1)), :Rps)
AnelasticAttenuationParameters(Q0::T, η::T) where {T<:Float64} = AnelasticAttenuationParameters([0.0, Inf], [Q0], Vector{T}(), [η], Vector{T}(), [3.5], BitVector(zeros(1)), BitVector(zeros(1)), :Rps)
AnelasticAttenuationParameters(Q0::T, η::T, rmetric::Symbol) where {T<:Float64} = AnelasticAttenuationParameters([0.0, Inf], [Q0], Vector{T}(), [η], Vector{T}(), [3.5], BitVector(zeros(1)), BitVector(zeros(1)), rmetric)
AnelasticAttenuationParameters(Q0::T, η::T) where {T<:Real} = AnelasticAttenuationParameters([0.0, Inf], Vector{Float64}(), [Q0], Vector{Float64}(), [η], [3.5], BitVector(ones(1)), BitVector(ones(1)), :Rps)
AnelasticAttenuationParameters(Q0::T, η::T, rmetric::Symbol) where {T<:Real} = AnelasticAttenuationParameters([0.0, Inf], Vector{Float64}(), [Q0], Vector{Float64}(), [η], [3.5], BitVector(ones(1)), BitVector(ones(1)), rmetric)
AnelasticAttenuationParameters(Q0::S, η::T) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters([0.0, Inf], [Q0], Vector{T}(), Vector{S}(), [η], [3.5], BitVector(zeros(1)), BitVector(ones(1)), :Rps)
AnelasticAttenuationParameters(Q0::S, η::T, rmetric::Symbol) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters([0.0, Inf], [Q0], Vector{T}(), Vector{S}(), [η], [3.5], BitVector(zeros(1)), BitVector(ones(1)), rmetric)
AnelasticAttenuationParameters(Q0::S, η::T) where {S<:Real,T<:Float64} = AnelasticAttenuationParameters([0.0, Inf], Vector{T}(), [Q0], [η], Vector{S}(), [3.5], BitVector(ones(1)), BitVector(zeros(1)), :Rps)
AnelasticAttenuationParameters(Q0::S, η::T, rmetric::Symbol) where {S<:Real,T<:Float64} = AnelasticAttenuationParameters([0.0, Inf], Vector{T}(), [Q0], [η], Vector{S}(), [3.5], BitVector(ones(1)), BitVector(zeros(1)), rmetric)

AnelasticAttenuationParameters(Q0, η, cQ::Float64) = AnelasticAttenuationParameters([0.0, Inf], [Q0], Vector{Float64}(), [η], Vector{Float64}(), [cQ], BitVector(zeros(1)), BitVector(zeros(1)), :Rps)
AnelasticAttenuationParameters(Q0, η, cQ::Float64, rmetric::Symbol) = AnelasticAttenuationParameters([0.0, Inf], [Q0], Vector{Float64}(), [η], Vector{Float64}(), [cQ], BitVector(zeros(1)), BitVector(zeros(1)), rmetric)

AnelasticAttenuationParameters(Rrefi::Vector{T}, Q0coni::Vector{T}) where {T<:Float64} = AnelasticAttenuationParameters{T,T,BitVector}(Rrefi, Q0coni, Vector{T}(), zeros(T, length(Q0coni)), Vector{T}(), 3.5 * ones(T, length(Q0coni)), BitVector(zeros(length(Q0coni))), BitVector(zeros(length(Q0coni))), :Rps)
AnelasticAttenuationParameters(Rrefi::Vector{T}, Q0coni::Vector{T}, rmetric::Symbol) where {T<:Float64} = AnelasticAttenuationParameters{T,T,BitVector}(Rrefi, Q0coni, Vector{T}(), zeros(T, length(Q0coni)), Vector{T}(), 3.5 * ones(T, length(Q0coni)), BitVector(zeros(length(Q0coni))), BitVector(zeros(length(Q0coni))), rmetric)
AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0vari::Vector{T}) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Vector{S}(), Q0vari, zeros(S, length(Q0vari)), Vector{S}(), 3.5 * ones(S, length(Q0vari)), BitVector(ones(length(Q0vari))), BitVector(zeros(length(Q0vari))), :Rps)
AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0vari::Vector{T}, rmetric::Symbol) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Vector{S}(), Q0vari, zeros(S, length(Q0vari)), Vector{S}(), 3.5 * ones(S, length(Q0vari)), BitVector(ones(length(Q0vari))), BitVector(zeros(length(Q0vari))), rmetric)

AnelasticAttenuationParameters(Rrefi::Vector{T}, Q0coni::Vector{T}, ηconi::Vector{T}) where {T<:Float64} = AnelasticAttenuationParameters{T,T,BitVector}(Rrefi, Q0coni, Vector{T}(), ηconi, Vector{T}(), 3.5 * ones(T, length(Q0coni)), BitVector(zeros(length(Q0coni))), BitVector(zeros(length(ηconi))), :Rps)
AnelasticAttenuationParameters(Rrefi::Vector{T}, Q0coni::Vector{T}, ηconi::Vector{T}, rmetric::Symbol) where {T<:Float64} = AnelasticAttenuationParameters{T,T,BitVector}(Rrefi, Q0coni, Vector{T}(), ηconi, Vector{T}(), 3.5 * ones(T, length(Q0coni)), BitVector(zeros(length(Q0coni))), BitVector(zeros(length(ηconi))), rmetric)
AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0vari::Vector{T}, ηvari::Vector{T}) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Vector{S}(), Q0vari, zeros(S, length(Q0vari)), ηvari, 3.5 * ones(S, length(Q0vari)), BitVector(ones(length(Q0vari))), BitVector(ones(length(ηvari))), :Rps)
AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0vari::Vector{T}, ηvari::Vector{T}, rmetric::Symbol) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Vector{S}(), Q0vari, zeros(S, length(Q0vari)), ηvari, 3.5 * ones(S, length(Q0vari)), BitVector(ones(length(Q0vari))), BitVector(ones(length(Q0vari))), rmetric)

AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0coni::Vector{S}, ηvari::Vector{T}) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Q0coni, Vector{T}(), Vector{S}(), ηvari, 3.5 * ones(S, length(Q0coni)), BitVector(zeros(length(Q0coni))), BitVector(ones(length(ηvari))), :Rps)
AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0coni::Vector{S}, ηvari::Vector{T}, rmetric::Symbol) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Q0coni, Vector{T}(), Vector{S}(), ηvari, 3.5 * ones(S, length(Q0coni)), BitVector(zeros(length(Q0coni))), BitVector(ones(length(ηvari))), rmetric)
AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0vari::Vector{T}, ηconi::Vector{S}) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Vector{S}(), Q0vari, ηconi, Vector{T}(), 3.5 * ones(S, length(Q0vari)), BitVector(ones(length(Q0vari))), BitVector(zeros(length(ηconi))), :Rps)
AnelasticAttenuationParameters(Rrefi::Vector{S}, Q0vari::Vector{T}, ηconi::Vector{S}, rmetric::Symbol) where {S<:Float64,T<:Real} = AnelasticAttenuationParameters{S,T,BitVector}(Rrefi, Vector{S}(), Q0vari, ηconi, Vector{T}(), 3.5 * ones(S, length(Q0vari)), BitVector(ones(length(Q0vari))), BitVector(zeros(length(ηconi))), rmetric)



"""
	get_parametric_type(anelastic::AnelasticAttenuationParameters{S,T,U}) where {S,T,U}

Extract type most elaborate type `S` or `T` from parametric `AnelasticAttenuationParameters` struct
"""
function get_parametric_type(anelastic::AnelasticAttenuationParameters{S,T,U}) where {S,T,U}
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

The base constructor is: `PathParameters(geo::G, sat::S, ane::A) where {G<:GeometricSpreadingParameters, S<:NearSourceSaturationParameters, A<:AnelasticAttenuationParameters}`

See also: [`FourierParameters`](@ref)
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



@doc raw"""
	SiteParameters

Custom type defining the site parameters of a Fourier spectrum

- `κ0::T where T<:Real` is the site kappa in units of s
- `ζ0::T where T<:Real` is the Haendel et al. (2020) ζ parameter (for a reference frequency of ``f_0=1`` Hz)
- `η::T where T<:Real` is the Haendel et al. (2020) η parameter
- `model::Symbol` is a symbol identifying the impedance function

The argument `model` is currently one of:
- `:Unit` for a generic unit amplification
- `:Boore2016` for the Boore (2016) amplification for ``V_{S,30}=760`` m/s
- `:AlAtik2021_ask14_620` for the Al Atik & Abrahamson (2021) inversion of ASK14 for ``V_{S,30}=620`` m/s
- `:AlAtik2021_ask14_760` for the Al Atik & Abrahamson (2021) inversion of ASK14 for ``V_{S,30}=760`` m/s
- `:AlAtik2021_ask14_1100` for the Al Atik & Abrahamson (2021) inversion of ASK14 for ``V_{S,30}=1100`` m/s
- `:AlAtik2021_bssa14_620` for the Al Atik & Abrahamson (2021) inversion of BSSA14 for ``V_{S,30}=620`` m/s
- `:AlAtik2021_bssa14_760` for the Al Atik & Abrahamson (2021) inversion of BSSA14 for ``V_{S,30}=760`` m/s
- `:AlAtik2021_bssa14_1100` for the Al Atik & Abrahamson (2021) inversion of BSSA14 for ``V_{S,30}=1100`` m/s
- `:AlAtik2021_cb14_620` for the Al Atik & Abrahamson (2021) inversion of CB14 for ``V_{S,30}=620`` m/s
- `:AlAtik2021_cb14_760` for the Al Atik & Abrahamson (2021) inversion of CB14 for ``V_{S,30}=760`` m/s
- `:AlAtik2021_cb14_1100` for the Al Atik & Abrahamson (2021) inversion of CB14 for ``V_{S,30}=1100`` m/s
- `:AlAtik2021_cy14_620` for the Al Atik & Abrahamson (2021) inversion of CY14 for ``V_{S,30}=620`` m/s
- `:AlAtik2021_cy14_760` for the Al Atik & Abrahamson (2021) inversion of CY14 for ``V_{S,30}=760`` m/s
- `:AlAtik2021_cy14_1100` for the Al Atik & Abrahamson (2021) inversion of CY14 for ``V_{S,30}=1100`` m/s

See also: [`FourierParameters`](@ref), [`site_amplification`](@ref)
"""
struct SiteParameters{T<:Real,U<:Real,V<:Real}
    # site parameters
    κ0::T           # site kappa
    ζ0::U           # zeta parameter for f0=1 Hz
    η::V            # eta parameter 
    model::Symbol   # site amplification model
end

SiteParameters(κ0::T) where {T} = SiteParameters(κ0, NaN, NaN, :AlAtik2021_cy14_760)
SiteParameters(κ0::T, model::Symbol) where {T} = SiteParameters(κ0, NaN, NaN, model)
SiteParameters(ζ0::T, η::U) where {T<:Real,U<:Real} = SiteParameters(NaN, ζ0, η, :AlAtik2021_cy14_760)
SiteParameters(ζ0::T, η::U, model::Symbol) where {T<:Real,U<:Real} = SiteParameters(NaN, ζ0, η, model)


"""
	get_parametric_type(site::SiteParameters{T}) where {T} = T

Extract type of `T` from parametric `SiteParameters` struct
"""
get_parametric_type(site::SiteParameters{T}) where {T} = T


"""
    FourierParameters

Custom type for the parameters the Fourier amplitude spectrum.
This type is comprised of source, path and site types, and so has a base constructor of: `FourierParameters(src::S, path::T, site::U) where {S<:SourceParameters, T<:PathParameters, U<:SiteParameters}`

See also: [`SourceParameters`](@ref), [`PathParameters`](@ref), [`SiteParameters`](@ref)
"""
struct FourierParameters{S<:SourceParameters,T<:PathParameters,U<:SiteParameters}
    source::S   # source parameters
    path::T     # path parameters
    site::U     # site parameters
end

# initialiser focussing upon source and path response only. Uses zero kappa and unit amplification
FourierParameters(src::S, path::T) where {S<:SourceParameters,T<:PathParameters} = FourierParameters(src, path, SiteParameters(0.0, :Unit))

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
