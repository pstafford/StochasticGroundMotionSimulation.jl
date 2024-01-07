
"""
	geometric_spreading(r_ps::T, m::S, geo::GeometricSpreadingParameters, sat::NearSourceSaturationParameters) where {S<:Real, T<:Real}

Geometric spreading function, switches between different approaches on `path.geo_model`.
"""
function geometric_spreading(r_ps::T, m::S, geo::GeometricSpreadingParameters, sat::NearSourceSaturationParameters) where {S<:Real,T<:Real}
    if geo.model == :Piecewise
        return geometric_spreading_piecewise(r_ps, geo)
    elseif geo.model == :CY14
        return geometric_spreading_cy14(r_ps, geo)
    elseif geo.model == :CY14mod
        return geometric_spreading_cy14mod(r_ps, m, geo, sat)
    else
        return NaN * oneunit(get_parametric_type(geo))
    end
end

geometric_spreading(r_ps::T, m::S, path::PathParameters) where {S<:Real,T<:Real} = geometric_spreading(r_ps, m, path.geometric, path.saturation)
geometric_spreading(r_ps::T, m::S, fas::FourierParameters) where {S<:Real,T<:Real} = geometric_spreading(r_ps, m, fas.path)
# special variants that don't require a magnitude input
geometric_spreading(r_ps::T, path::PathParameters) where {T<:Real} = geometric_spreading(r_ps, NaN, path.geometric, path.saturation)
geometric_spreading(r_ps::T, fas::FourierParameters) where {T<:Real} = geometric_spreading(r_ps, fas.path)



"""
	geometric_spreading_piecewise(r_ps::W, geo::GeometricSpreadingParameters{S,T,U,V}) where {S<:Real, T<:Real, U<:Real, V<:AbstractVector{Bool}, W<:Real}

Piecewise linear (in log-log space) geometric spreading function.
Makes use of the reference distances `Rrefi` and spreading rates `γi` in `path`.
"""
function geometric_spreading_piecewise(r_ps::W, geo::GeometricSpreadingParameters{S,T,U,V}) where {S<:Real,T<:Real,U<:Real,V<:AbstractVector{Bool},W<:Real}
    z_r = oneunit(get_parametric_type(geo))
    j = 1
    k = 1
    for i = 1:length(geo.γfree)
        @inbounds Rr0 = geo.Rrefi[i]
        @inbounds Rr1 = geo.Rrefi[i+1]
        γ_r = oneunit(get_parametric_type(geo))

        if geo.γfree[i] == 0
            @inbounds γ_r *= geo.γconi[j]
            j += 1
        else
            @inbounds γ_r *= geo.γvari[k]
            k += 1
        end

        if r_ps < Rr1
            z_r *= (Rr0 / r_ps)^γ_r
            return z_r
        else
            z_r *= (Rr0 / Rr1)^γ_r
        end
    end
    return z_r
end

geometric_spreading_piecewise(r_ps, path::PathParameters) = geometric_spreading_piecewise(r_ps, path.geometric)
geometric_spreading_piecewise(r_ps, fas::FourierParameters) = geometric_spreading_piecewise(r_ps, fas.path)


"""
	geometric_spreading_cy14(r_ps::W, geo::GeometricSpreadingParameters{S,T,U,V}) where {S<:Real, T<:Real, U<:Real, V<:AbstractVector{Bool}, W<:Real}

Geometric spreading function from Chiou & Youngs (2014).
Defines a smooth transition from one rate `γi[1]` to another `γi[2]`, with a spreading bandwidth of `Rrefi[2]` km.
"""
function geometric_spreading_cy14(r_ps::W, geo::GeometricSpreadingParameters{S,T,U,V}) where {S<:Real,T<:Real,U<:Real,V<:AbstractVector{Bool},W<:Real}
    unit = oneunit(get_parametric_type(geo))
    j = 1
    k = 1
    if geo.γfree[1] == 0
        γ1 = geo.γconi[j] * unit
        j += 1
    else
        γ1 = geo.γvari[k]
        k += 1
    end
    if geo.γfree[2] == 0
        γ2 = geo.γconi[j] * unit
    else
        γ2 = geo.γvari[k]
    end
    R0sq = (geo.Rrefi[1])^2
    Rrsq = (geo.Rrefi[2])^2
    ln_z_r = -γ1 * log(r_ps) + (-γ2 + γ1) / 2 * log((r_ps^2 + Rrsq) / (R0sq + Rrsq))
    z_r = exp(ln_z_r)
    return z_r
end


geometric_spreading_cy14(r_ps, path::PathParameters) = geometric_spreading_cy14(r_ps, path.geometric)
geometric_spreading_cy14(r_ps, fas::FourierParameters) = geometric_spreading_cy14(r_ps, fas.path)


"""
	geometric_spreading_cy14mod(r_ps::W, m::X, geo::GeometricSpreadingParameters{S,T,U,V}, sat::NearSourceSaturationParameters) where {S<:Real, T<:Real, U<:Real, V<:AbstractVector{Bool}, W<:Real, X<:Real}

Geometric spreading function from Chiou & Youngs (2014).
Modified to make use of both `r_ps` and `r_rup` so that only the first saturation term contaminates the source amplitudes.
Defines a smooth transition from one rate `γi[1]` to another `γi[2]`, with a spreading bandwidth of `Rrefi[2]` km.
"""
function geometric_spreading_cy14mod(r_ps::W, m::X, geo::GeometricSpreadingParameters{S,T,U,V}, sat::NearSourceSaturationParameters) where {S<:Real,T<:Real,U<:Real,V<:AbstractVector{Bool},W<:Real,X<:Real}
    unit = oneunit(get_parametric_type(geo))
    j = 1
    k = 1
    if geo.γfree[1] == 0
        γ1 = geo.γconi[j] * unit
        j += 1
    else
        γ1 = geo.γvari[k]
        k += 1
    end
    if geo.γfree[2] == 0
        γ2 = geo.γconi[j] * unit
    else
        γ2 = geo.γvari[k]
    end
    R0sq = (geo.Rrefi[1])^2
    Rrsq = (geo.Rrefi[2])^2
    r_rup = rupture_distance_from_equivalent_point_source_distance(r_ps, m, sat)
    ln_z_r = -γ1 * log(r_ps) + (-γ2 + γ1) / 2 * log((r_rup^2 + Rrsq) / (R0sq + Rrsq))
    z_r = exp(ln_z_r)
    return z_r
end


geometric_spreading_cy14mod(r_ps, m, path::PathParameters) = geometric_spreading_cy14mod(r_ps, m, path.geometric, path.saturation)
geometric_spreading_cy14mod(r_ps, m, fas::FourierParameters) = geometric_spreading_cy14mod(r_ps, m, fas.path)



"""
	near_source_saturation(m, sat::NearSourceSaturationParameters)

Near-source saturation term. Used to create equivalent point-source distance.
Switches methods based upon `sat.model`.

# Arguments
Options for `sat.model` are:
- `:BT15` for Boore & Thompson (2015) finite fault factor
- `:YA15` for Yenier & Atkinson (2014) finite fault factor
- `:CY14` for a model fitted to the Chiou & Youngs (2014) saturation lengths (over all periods)
- `:SEA21` for the Stafford et al. (2021) saturation model obtained from inversion of Chiou & Youngs (2014)
- `:None` for zero saturation length
- `:ConstantConstrained` for a constant value, `sat.hconi[1]`, not subject to AD operations
- `:ConstantVariable` for a constant value, `sat.hvari[1]`, that is subject to AD operations

Any other symbol passed will return `NaN`.

See also: [`near_source_saturation`](@ref)
"""
function near_source_saturation(m, sat::NearSourceSaturationParameters)
    unit = oneunit(get_parametric_type(sat))
    if sat.model == :BT15
        # use the Boore & Thompson (2015) finite fault factor
        return finite_fault_factor_bt15(m) * unit
    elseif sat.model == :YA15
        # use the Yenier & Atkinson (2015) finite fault factor
        return finite_fault_factor_ya15(m) * unit
    elseif sat.model == :CY14
        # use the fitted model averaging over Chiou & Youngs (2014)
        return finite_fault_factor_cy14(m) * unit
    elseif sat.model == :SEA21 
        # use the near-source saturation model of Stafford et al. (2021)'s inversion of CY14 
        return finite_fault_factor_sea21(m) * unit
    elseif sat.model == :None
        return zero(get_parametric_type(sat))
    elseif sat.model == :ConstantConstrained
        return sat.hconi[1] * unit
    elseif sat.model == :ConstantVariable
        return sat.hvari[1] * unit
    else
        return NaN * unit
    end
end

"""
	near_source_saturation(m, path::PathParameters)

Near-source saturation term taking a `PathParameters` struct.

See also: [`near_source_saturation`](@ref)
"""
near_source_saturation(m, path::PathParameters) = near_source_saturation(m, path.saturation)

"""
	near_source_saturation(m, fas::FourierParameters)

Near-source saturation term taking a `FourierParameters` struct.

See also: [`near_source_saturation`](@ref)
"""
near_source_saturation(m, fas::FourierParameters) = near_source_saturation(m, fas.path)



"""
	finite_fault_factor_bt15(m::T) where T<:Real

Finite fault factor from Boore & Thompson (2015) used to create equivalent point-source distance.

See also: [`near_source_saturation`](@ref), [`finite_fault_factor_ya15`](@ref), [`finite_fault_factor_cy14`](@ref), [`finite_fault_factor_sea21`](@ref)
"""
function finite_fault_factor_bt15(m::T) where {T<:Real}
    Mt1 = 5.744
    Mt2 = 7.744
    if m <= Mt1
        a1 = 0.7497
        b1 = 0.4300
        h = a1 + b1 * (m - Mt1)
    elseif m >= Mt2
        a2 = 1.4147
        b2 = 0.2350
        h = a2 + b2 * (m - Mt2)
    else
        c0 = 0.7497
        c1 = 0.4300
        c2 = -0.04875
        c3 = 0.0
        h = c0 + c1 * (m - Mt1) + c2 * (m - Mt1)^2 + c3 * (m - Mt1)^3
    end
    return 10.0^h
end

"""
	finite_fault_factor_ya15(m::T) where T<:Real

Finite fault factor from Yenier & Atkinson (2015) used to create equivalent point-source distance.

See also: [`near_source_saturation`](@ref), [`finite_fault_factor_bt15`](@ref), [`finite_fault_factor_cy14`](@ref), [`finite_fault_factor_sea21`](@ref)
"""
function finite_fault_factor_ya15(m::T) where {T<:Real}
    return 10.0^(-0.405 + 0.235 * m)
end


"""
	finite_fault_factor_cy14(m::T) where T<:Real

Finite fault factor for Chiou & Youngs (2014) used to create equivalent point-source distance.
This is a model developed to match the average saturation length over the full period range.

See also: [`near_source_saturation`](@ref), [`finite_fault_factor_bt15`](@ref), [`finite_fault_factor_ya15`](@ref), [`finite_fault_factor_sea21`](@ref)
"""
function finite_fault_factor_cy14(m::T) where {T<:Real}
    hα = 7.308
    hβ = 0.4792
    hγ = 3.556
    return hα * cosh(hβ * max(m - hγ, 0.0))
end


"""
	finite_fault_factor_sea21(m::T) where T<:Real

Finite fault factor for Stafford et al. (2021) used to create equivalent point-source distance.

See also: [`near_source_saturation`](@ref), [`finite_fault_factor_bt15`](@ref), [`finite_fault_factor_ya15`](@ref), [`finite_fault_factor_cy14`](@ref)
"""
function finite_fault_factor_sea21(m::T) where {T<:Real}
    hα = -0.8712
    hβ = 0.4451
    hγ = 1.1513
    hδ = 5.0948
    hε = 7.2725
    lnh = hα + hβ * m + ((hβ - hγ)/hδ) * log( 1.0 + exp(-hδ * (m - hε)) )
    return exp(lnh)
end


"""
	equivalent_point_source_distance(r, m, sat::NearSourceSaturationParameters)

Compute equivalent point source distance
- `r` is ``r_{hyp}`` or ``r_{rup}`` (depending upon the size of the event -- but is nominally ``r_{rup}``)
- `m` is magnitude
- `sat` contains the `NearSourceSaturationParameters`

See also: [`near_source_saturation`](@ref)
"""
function equivalent_point_source_distance(r, m, sat::NearSourceSaturationParameters)
    h = near_source_saturation(m, sat)
    n = sat.exponent
    return (r^n + h^n)^(1 / n)
end

equivalent_point_source_distance(r, m, path::PathParameters) = equivalent_point_source_distance(r, m, path.saturation)
equivalent_point_source_distance(r, m, fas::FourierParameters) = equivalent_point_source_distance(r, m, fas.path)


"""
	rupture_distance_from_equivalent_point_source_distance(r_ps, m, sat::NearSourceSaturationParameters)

Compute rupture distance from equivalent point source distance
- `r_ps` is the equivalent point source distance
- `m` is magnitude
- `sat` contains the `NearSourceSaturationParameters`
"""
function rupture_distance_from_equivalent_point_source_distance(r_ps, m, sat::NearSourceSaturationParameters)
    h = near_source_saturation(m, sat)
    n = sat.exponent
    return (r_ps^n - h^n)^(1 / n)
end

rupture_distance_from_equivalent_point_source_distance(r_ps, m, path::PathParameters) = rupture_distance_from_equivalent_point_source_distance(r_ps, m, path.saturation)
rupture_distance_from_equivalent_point_source_distance(r_ps, m, fas::FourierParameters) = rupture_distance_from_equivalent_point_source_distance(r_ps, m, fas.path)


"""
	anelastic_attenuation(f::S, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real}

Anelastic attenuation filter, computed using equivalent point source distance metric or a standard rupture distance.
"""
function anelastic_attenuation(f::S, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real}
    PT = get_parametric_type(anelastic)
    Qfilt = oneunit(PT)
    jq = jη = 1
    kq = kη = 1
    nsegs = length(anelastic.Qfree)
    @inbounds for i = 1:nsegs
        Rr0 = anelastic.Rrefi[i]
        Rr1 = anelastic.Rrefi[i+1]
        cQ_r = anelastic.cQ[i]

        # get the relevant constrained or free quality factor for this path segment
        if anelastic.Qfree[i] == 0
            Q0_r = anelastic.Q0coni[jq]
            jq += 1
        else
            Q0_r = anelastic.Q0vari[kq]
            kq += 1
        end
        # get the relevant constrained of free quality exponent for this path segment
        if anelastic.ηfree[i] == 0
            η_r = anelastic.ηconi[jη]
            jη += 1
        else
            η_r = anelastic.ηvari[kη]
            kη += 1
        end

        if nsegs == 1
            fpow = f^(1.0 - η_r)
        else
            # this avoids f^(1-η)
            fpow = exp((1.0 - η_r) * log(f))
        end
        if r < Rr1
            Qfilt *= exp(-π * fpow * (r - Rr0) / (Q0_r * cQ_r))
            return Qfilt
        else
            Qfilt *= exp(-π * fpow * (Rr1 - Rr0) / (Q0_r * cQ_r))
        end
    end
    return PT(NaN)
end

anelastic_attenuation(f, r, path::PathParameters) = anelastic_attenuation(f, r, path.anelastic)
anelastic_attenuation(f, r, fas::FourierParameters) = anelastic_attenuation(f, r, fas.path)

"""
	anelastic_attenuation(f::Vector{S}, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real}

Anelastic attenuation filter, computed using equivalent point source distance metric or a standard rupture distance.
"""
function anelastic_attenuation(f::Vector{S}, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real}
    PT = promote_type(get_parametric_type(anelastic), T)
    nf = length(f)
    Qfilt = ones(PT, nf)
    jq = jη = 1
    kq = kη = 1
    nsegs = length(anelastic.Qfree)
    @inbounds for i = 1:nsegs
        Rr0 = anelastic.Rrefi[i]
        Rr1 = anelastic.Rrefi[i+1]
        cQ_r = anelastic.cQ[i]

        # get the relevant constrained or free quality factor for this path segment
        if anelastic.Qfree[i] == 0
            Q0_r = anelastic.Q0coni[jq]
            jq += 1
        else
            Q0_r = anelastic.Q0vari[kq]
            kq += 1
        end
        # get the relevant constrained of free quality exponent for this path segment
        if anelastic.ηfree[i] == 0
            η_r = anelastic.ηconi[jη]
            jη += 1
        else
            η_r = anelastic.ηvari[kη]
            kη += 1
        end

        if nsegs == 1
            fpow = f .^ (1.0 - η_r)
        else
            # this avoids f^(1-η)
            fpow = @. exp((1.0 - η_r) * log(f))
        end
        if r < Rr1
            for i in 1:nf
                Qfilt[i] *= exp(-π * fpow[i] * (r - Rr0) / (Q0_r * cQ_r))
            end
            # Qfilt .*= @. exp(-π * fpow * (r - Rr0) / (Q0_r * cQ_r))
            return Qfilt
        else
            for i in 1:nf
                Qfilt[i] *= exp(-π * fpow[i] * (Rr1 - Rr0) / (Q0_r * cQ_r))
            end
            # Qfilt .*= @. exp(-π * fpow * (Rr1 - Rr0) / (Q0_r * cQ_r))
        end
    end
    return PT(NaN)
end

"""
	anelastic_attenuation!(Qfilt::Vector{U}, f::Vector{S}, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real,U<:Real}

Anelastic attenuation filter, computed using equivalent point source distance metric or a standard rupture distance.
"""
function anelastic_attenuation!(Qfilt::Vector{U}, f::Vector{S}, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real,U<:Real}
    length(Qfilt) == length(f) || error("length of frequency vector must match the filter vector length")
    jq = jη = 1
    kq = kη = 1
    nsegs = length(anelastic.Qfree)
    for i = 1:nsegs
        @inbounds Rr0 = anelastic.Rrefi[i]
        @inbounds Rr1 = anelastic.Rrefi[i+1]
        @inbounds cQ_r = anelastic.cQ[i]

        # get the relevant constrained or free quality factor for this path segment
        if anelastic.Qfree[i] == 0
            @inbounds Q0_r = anelastic.Q0coni[jq]
            jq += 1
        else
            @inbounds Q0_r = anelastic.Q0vari[kq]
            kq += 1
        end
        # get the relevant constrained of free quality exponent for this path segment
        if anelastic.ηfree[i] == 0
            @inbounds η_r = anelastic.ηconi[jη]
            jη += 1
        else
            @inbounds η_r = anelastic.ηvari[kη]
            kη += 1
        end

        if nsegs == 1
            fpow = @. f^(1.0 - η_r)
        else
            # this avoids f^(1-η)
            fpow = @. exp((1.0 - η_r) * log(f))
        end
        if r < Rr1
            if i == 1
                Qfilt .= @. exp(-π * fpow * (r - Rr0) / (Q0_r * cQ_r))
            else
                Qfilt .*= @. exp(-π * fpow * (r - Rr0) / (Q0_r * cQ_r))
            end
            return Qfilt
        else
            if i == 1
                Qfilt .= @. exp(-π * fpow * (Rr1 - Rr0) / (Q0_r * cQ_r))
            else
                Qfilt .*= @. exp(-π * fpow * (Rr1 - Rr0) / (Q0_r * cQ_r))
            end
        end
    end
end

"""
	apply_anelastic_attenuation!(Af::Vector{T}, f::Vector{U}, r::V, anelastic::AnelasticAttenuationParameters) where {T<:Real,U<:Real,V<:Real}

Apply an anelastic attenuation filter to a FAS, computed using equivalent point source distance metric or a standard rupture distance.
"""
function apply_anelastic_attenuation!(Af::Vector{T}, f::Vector{U}, r::V, anelastic::AnelasticAttenuationParameters) where {T<:Real,U<:Real,V<:Real}
    numAf = length(Af)
    numf = length(f)
    numAf == numf || error("length of vector `f` must match the length of vector `Af`")
    jq = jη = 1
    kq = kη = 1
    nsegs = length(anelastic.Qfree)
    exp_arg = zeros(T, numf)
    for i = 1:nsegs
        @inbounds Rr0 = anelastic.Rrefi[i]
        @inbounds Rr1 = anelastic.Rrefi[i+1]
        @inbounds cQ_r = anelastic.cQ[i]

        # get the relevant constrained or free quality factor for this path segment
        @inbounds if anelastic.Qfree[i] == 0
            @inbounds Q0_r = anelastic.Q0coni[jq]
            jq += 1
        else
            @inbounds Q0_r = anelastic.Q0vari[kq]
            kq += 1
        end
        # get the relevant constrained of free quality exponent for this path segment
        @inbounds if anelastic.ηfree[i] == 0
            @inbounds η_r = anelastic.ηconi[jη]
            jη += 1
        else
            @inbounds η_r = anelastic.ηvari[kη]
            kη += 1
        end

        if η_r ≈ 0.0
            fpow = f
        else
            if nsegs == 1
                fpow = @. f^(1.0 - η_r)
            else
                # this avoids f^(1-η)
                fpow = @. exp((1.0 - η_r) * log(f))
            end
        end
        if r < Rr1
            rlim = r
        else
            rlim = Rr1
        end
        for (j, fp) in pairs(fpow)
            @inbounds exp_arg[j] += -π * fp * (rlim - Rr0) / (Q0_r * cQ_r)
        end
        # for (j, fp) in pairs(fpow)
        #     @inbounds Af[j] *= exp(-π * fp * (rlim - Rr0) / (Q0_r * cQ_r))
        # end
    end
    for (i, ea) in pairs(exp_arg)
        @inbounds Af[i] *= exp(ea)
    end
    return nothing
end

"""
	effective_quality_parameters(r::T, anelastic::AnelasticAttenuationParameters) where {T<:Real}

Distance weighted quality parameters for segmented model.
Returns a tuple of effective Q0, η & cQ values from a weighted sum, weighted by relative distance in each path segment.
"""
function effective_quality_parameters(r::T, anelastic::AnelasticAttenuationParameters) where {T<:Real}
    PT = get_parametric_type(anelastic)
    Q_eff = zero(PT)
    η_eff = zero(PT)
    cQ_eff = zero(PT)
    jq = jη = 1
    kq = kη = 1
    for i = 1:length(anelastic.Qfree)
        @inbounds Rr0 = anelastic.Rrefi[i]
        @inbounds Rr1 = anelastic.Rrefi[i+1]

        Q0_r = oneunit(PT)
        η_r = oneunit(PT)
        cQ_r = anelastic.cQ[i]

        # get the relevant constrained or free quality factor for this path segment
        if anelastic.Qfree[i] == 0
            @inbounds Q0_r *= anelastic.Q0coni[jq]
            jq += 1
        else
            @inbounds Q0_r *= anelastic.Q0vari[kq]
            kq += 1
        end
        # get the relevant constrained of free quality exponent for this path segment
        if anelastic.ηfree[i] == 0
            @inbounds η_r *= anelastic.ηconi[jη]
            jη += 1
        else
            @inbounds η_r *= anelastic.ηvari[kη]
            kη += 1
        end

        if r < Rr1
            rfrac = (r - Rr0) / r
            Q_eff += rfrac * Q0_r
            η_eff += rfrac * η_r
            cQ_eff += rfrac * cQ_r
            return (Q_eff, η_eff, cQ_eff)
        else
            rfrac = (Rr1 - Rr0) / r
            Q_eff += rfrac * Q0_r
            η_eff += rfrac * η_r
            cQ_eff += rfrac * cQ_r
        end
    end
end
