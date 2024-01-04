
"""
	site_amplification(f::T, amp_model::S; fmin=1e-3, fmax=999.99) where {T<:Real,S<:SiteAmplification}

Computes the site amplification (impedance) for a given frequency `f`. The argument `amp_model` is a subtype of the abstract type `SiteAmplification`.

# Examples
```julia-repl
	f = 5.0
	# returns the amplification from AlAtik (2021) in both cases
	Af = site_amplification(f)
    Af = site_amplification(f, SiteAmpAlAtikAbrahamson2021_cy14_760())
    # returns the Boore (2016) amplification
	Af = site_amplification(f, SiteAmpBoore2016_760())
	# returns 1.0
	Af = site_amplification(f, SiteAmpUnit())
```
"""
function site_amplification(f::T, amp_model::S; fmin=1e-3, fmax=999.99) where {T<:Real,S<:SiteAmplification}
    if f < fmin
        f = fmin
    elseif f > fmax
        f = fmax
    end
    return amp_model.amplification(f)::T
end

"""
	site_amplification(f::Vector{T}, amp_model::S; fmin=1e-3, fmax=999.99) where {T<:Real,S<:SiteAmplification}

Computes the site amplification (impedance) for a given frequency `f`. The argument `amp_model` is a subtype of the abstract type `SiteAmplification`.

# Examples
```julia-repl
	f = 5.0
	# returns the amplification from AlAtik (2021) in both cases
	Af = site_amplification(f)
    Af = site_amplification(f, SiteAmpAlAtikAbrahamson2021_cy14_760())
    # returns the Boore (2016) amplification
	Af = site_amplification(f, SiteAmpBoore2016_760())
	# returns 1.0
	Af = site_amplification(f, SiteAmpUnit())
```
"""
function site_amplification(f::Vector{T}, amp_model::S; fmin=1e-3, fmax=999.99) where {T<:Real,S<:SiteAmplification}
    clamp!(f, fmin, fmax)
    return amp_model.amplification.(f)::Vector{T}
end


"""
	site_amplification(f, site::SiteParameters)

Computes the site amplification (impedance) for a given frequency `f`.
"""
site_amplification(f, site::SiteParameters) = site_amplification(f, site.model)

"""
	site_amplification(f, fas::FourierParameters)

Computes the site amplification (impedance) for a given frequency `f`.
"""
site_amplification(f, fas::FourierParameters) = site_amplification(f, fas.site)


"""
    kappa_filter(f, site::SiteParameters)

Kappa filter for a given frequency `f`
"""
function kappa_filter(f, site::SiteParameters)
    if isnan(site.κ0)
        # Eq. 20 of Haendel et al. (2020) (frequency dependent kappa)
        # specified for a reference frequency of f0=1 Hz
        return exp(-π * site.ζ0 * (1.0 - site.η) * f^(1.0 - site.η))
    else
        return exp(-π * f * site.κ0)
    end
end
