
"""
	geometric_spreading(r_ps::T, m::S, geo::GeometricSpreadingParameters, sat::NearSourceSaturationParameters) where {S<:Float64, T<:Real}

Geometric spreading function, switches between different approaches on `path.geo_model`.
"""
function geometric_spreading(r_ps::T, m::S, geo::GeometricSpreadingParameters, sat::NearSourceSaturationParameters) where {S<:Float64, T<:Real}
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

geometric_spreading(r_ps::T, m::S, path::PathParameters) where {S<:Float64, T<:Real} = geometric_spreading(r_ps, m, path.geometric, path.saturation)
geometric_spreading(r_ps::T, m::S, fas::FourierParameters) where {S<:Float64, T<:Real} = geometric_spreading(r_ps, m, fas.path)
# special variants that don't require a magnitude input
geometric_spreading(r_ps::T, path::PathParameters) where T<:Real = geometric_spreading(r_ps, NaN, path.geometric, path.saturation)
geometric_spreading(r_ps::T, fas::FourierParameters) where T<:Real = geometric_spreading(r_ps, fas.path)



"""
	geometric_spreading_piecewise(r_ps::S, geo::GeometricSpreadingParameters{T,T,U}) where {S<:Real, T<:Float64, U<:AbstractVector{Bool}}

Piecewise linear (in log-log space) geometric spreading function.
Makes use of the reference distances `Rrefi` and spreading rates `γi` in `path`.
"""
function geometric_spreading_piecewise(r_ps::S, geo::GeometricSpreadingParameters{T,T,U}) where {S<:Real, T<:Float64, U<:AbstractVector{Bool}}
	z_r = oneunit(T)
    for i in 1:length(geo.γfree)
		@inbounds Rr0 = geo.Rrefi[i]
		@inbounds Rr1 = geo.Rrefi[i+1]
		@inbounds γ_r = geo.γconi[i]
    	if r_ps < Rr1
        	z_r *= (Rr0 / r_ps)^γ_r
        	return z_r
      	else
        	z_r *= (Rr0 / Rr1)^γ_r
      	end
    end
    return z_r
end

"""
	geometric_spreading_piecewise(r_ps::V, geo::GeometricSpreadingParameters{S,T,U}) where {S<:Float64, T<:Real, U<:AbstractVector{Bool}, V<:Real}

Piecewise linear (in log-log space) geometric spreading function.
Makes use of the reference distances `Rrefi` and spreading rates `γi` in `path`.
"""
function geometric_spreading_piecewise(r_ps::V, geo::GeometricSpreadingParameters{S,T,U}) where {S<:Float64, T<:Real, U<:AbstractVector{Bool}, V<:Real}
	z_r = oneunit(T)
	j = 1
	k = 1
    for i in 1:length(geo.γfree)
		@inbounds Rr0 = geo.Rrefi[i]
		@inbounds Rr1 = geo.Rrefi[i+1]
		γ_r = oneunit(T)

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
	geometric_spreading_cy14(r_ps::V, geo::GeometricSpreadingParameters{S,T,U}) where {S<:Float64, T<:Real, U<:AbstractVector{Bool}, V<:Real}

Geometric spreading function from Chiou & Youngs (2014).
Defines a smooth transition from one rate `γi[1]` to another `γi[2]`, with a spreading bandwidth of `Rrefi[2]` km.
"""
function geometric_spreading_cy14(r_ps::V, geo::GeometricSpreadingParameters{S,T,U}) where {S<:Float64, T<:Real, U<:AbstractVector{Bool}, V<:Real}
	unit = oneunit(T)
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
	ln_z_r = -γ1*log(r_ps) + (-γ2+γ1)/2*log( (r_ps^2 + Rrsq) / (R0sq + Rrsq) )
	z_r = exp(ln_z_r)
	return z_r
end


geometric_spreading_cy14(r_ps, path::PathParameters) = geometric_spreading_cy14(r_ps, path.geometric)
geometric_spreading_cy14(r_ps, fas::FourierParameters) = geometric_spreading_cy14(r_ps, fas.path)


"""
	geometric_spreading_cy14mod(r_ps::V, m::S, geo::GeometricSpreadingParameters{S,T,U}, sat::NearSourceSaturationParameters) where {S<:Float64, T<:Real, U<:AbstractVector{Bool}, V<:Real}

Geometric spreading function from Chiou & Youngs (2014).
Modified to make use of both `r_ps` and `r_rup` so that only the first saturation term contaminates the source amplitudes.
Defines a smooth transition from one rate `γi[1]` to another `γi[2]`, with a spreading bandwidth of `Rrefi[2]` km.
"""
function geometric_spreading_cy14mod(r_ps::V, m::S, geo::GeometricSpreadingParameters{S,T,U}, sat::NearSourceSaturationParameters) where {S<:Float64, T<:Real, U<:AbstractVector{Bool}, V<:Real}
	unit = oneunit(T)
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
	ln_z_r = -γ1*log(r_ps) + (-γ2+γ1)/2*log( (r_rup^2 + Rrsq) / (R0sq + Rrsq) )
	z_r = exp(ln_z_r)
	return z_r
end


geometric_spreading_cy14mod(r_ps, m, path::PathParameters) = geometric_spreading_cy14mod(r_ps, m, path.geometric, path.saturation)
geometric_spreading_cy14mod(r_ps, m, fas::FourierParameters) = geometric_spreading_cy14mod(r_ps, m, fas.path)



"""
	near_source_saturation(m::Float64, sat::NearSourceSaturationParameters)

Near-source saturation term. Used to create equivalent point-source distance. Switches methods based upon `sat.model`.
"""
function near_source_saturation(m::Float64, sat::NearSourceSaturationParameters)
	unit = oneunit(get_parametric_type(sat))
	if sat.model == :BT15
		# use the Boore & Thompson (2015) finite fault factor
		return finite_fault_factor_bt15(m) * unit
	elseif sat.model == :YA15
		# use the Yenier & Atkinson (2015) finite fault factor
		return finite_fault_factor_ya15(m) * unit
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

near_source_saturation(m::Float64, path::PathParameters) = near_source_saturation(m, path.saturation)
near_source_saturation(m::Float64, fas::FourierParameters) = near_source_saturation(m, fas.path)



"""
	finite_fault_factor_bt15(m::Float64)

Finite fault factor from Boore & Thompson (2015) used to create equivalent point-source distance.
"""
function finite_fault_factor_bt15(m::Float64)
	Mt1 = 5.744
	Mt2 = 7.744
	if m <= Mt1
		a1 = 0.7497
		b1 = 0.4300
		h = a1 + b1*(m - Mt1)
	elseif m >= Mt2
		a2 = 1.4147
		b2 = 0.2350
		h = a2 + b2*(m - Mt2)
	else
		c0 = 0.7497
		c1 = 0.4300
		c2 = -0.04875
		c3 = 0.0
		h = c0 + c1*(m - Mt1) + c2*(m - Mt1)^2 + c3*(m - Mt1)^3
	end
	return 10.0^h
end

"""
	finite_fault_factor_ya15(m::Float64)

Finite fault factor from Yenier & Atkinson (2015) used to create equivalent point-source distance.
"""
function finite_fault_factor_ya15(m::Float64)
	return 10.0^(-0.405 + 0.235*m)
end


"""
	equivalent_point_source_distance(r, m, sat::NearSourceSaturationParameters)

Compute equivalent point source distance
- `r` is r_hyp or r_rup (depending upon the size of the event -- but is nominally r_rup)
- `m` is magnitude
- `sat` contains the `NearSourceSaturationParameters`
"""
function equivalent_point_source_distance(r, m, sat::NearSourceSaturationParameters)
	h = near_source_saturation(m, sat)
	n = sat.exponent
	return ( r^n + h^n )^(1/n)
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
	return ( r_ps^n - h^n )^(1/n)
end

rupture_distance_from_equivalent_point_source_distance(r_ps, m, path::PathParameters) = rupture_distance_from_equivalent_point_source_distance(r_ps, m, path.saturation)
rupture_distance_from_equivalent_point_source_distance(r_ps, m, fas::FourierParameters) = rupture_distance_from_equivalent_point_source_distance(r_ps, m, fas.path)


"""
	anelastic_attenuation(f::S, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real}

Anelastic attenuation filter, computed using equivalent point source distance metric or a standard rupture distance.
"""
function anelastic_attenuation(f::S, r::T, anelastic::AnelasticAttenuationParameters) where {S<:Float64,T<:Real}
	return exp( -π*f*r / ( anelastic.Q0 * f^anelastic.η * anelastic.cQ ) )
end

anelastic_attenuation(f, r, path::PathParameters) = anelastic_attenuation(f, r, path.anelastic)
anelastic_attenuation(f, r, fas::FourierParameters) = anelastic_attenuation(f, r, fas.path)
