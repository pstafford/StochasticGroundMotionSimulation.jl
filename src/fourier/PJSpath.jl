
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

# path = PathParameters([1.0, 50.0, Inf], Union{Float64, Dual{Float64}}[Dual{Float64}(1.0), 0.5], 200.0)
# typeof(path.γi)
# typeof(path.γi[1])
# typeof(path.γi[2])
# path = PathParameters([1.0, 50.0, Inf], [1.0, Dual(0.5)], 200.0)
# typeof(path.γi[1])
# typeof(path.γi[2])
# path.γi[1].value



"""
	geometric_spreading(r::Float64, path::PathParameters)

Geometric spreading function, switches between different approaches on `path.geo_model`.
"""
function geometric_spreading(r::Float64, path::PathParameters)
	if path.geo_model == :Piecewise
		return geometric_spreading_piecewise(r, path)
	elseif path.geo_model == :CY14
		return geometric_spreading_cy14(r, path)
	else
		return NaN
	end
end

geometric_spreading(r::Float64, fas::FourierParameters) = geometric_spreading(r, fas.path)


"""
	geometric_spreading_piecewise(r::Float64, path::PathParameters)

Piecewise linear (in log-log space) geometric spreading function.
Makes use of the reference distances `Rrefi` and spreading rates `γi` in `path`.
"""
function geometric_spreading_piecewise(r::Float64, path::PathParameters)
	z_r = 1.0
    for i = 2:length(path.Rrefi)
      @inbounds if r < path.Rrefi[i]
        @inbounds z_r *= (path.Rrefi[i-1] / r)^path.γi[i-1]
        return z_r
      else
        @inbounds z_r *= (path.Rrefi[i-1] / path.Rrefi[i])^path.γi[i-1]
      end
    end
    return z_r
end

geometric_spreading_piecewise(r::Float64, fas::FourierParameters) = geometric_spreading_piecewise(r, fas.path)

"""
	geometric_spreading_cy14(r::Float64, path::PathParameters)

Geometric spreading function from Chiou & Youngs (2014).
Defines a smooth transition from one rate `γi[1]` to another `γi[2]`, with a spreading bandwidth of 50 km.
"""
function geometric_spreading_cy14(r::Float64, path::PathParameters)
    ln_z_r = -path.γi[1]*log(r) + (-path.γi[2]+path.γi[1])*log(sqrt(r^2 + 50.0^2)) - (-path.γi[2]+path.γi[1])*log(sqrt(1.0^2 + 50.0^2))
	z_r = exp(ln_z_r)
	return z_r
end

geometric_spreading_cy14(r::Float64, fas::FourierParameters) = geometric_spreading_cy14(r, fas.path)


"""
	near_source_saturation(m::Float64, path::PathParameters)

Near-source saturation term. Used to create equivalent point-source distance. Switches methods based upon `path.sat_model`.
"""
function near_source_saturation(m::Float64, path::PathParameters)
	if path.sat_model == :BT15
		# use the Boore & Thompson (2015) finite fault factor
		return finite_fault_factor_bt15(m)
	elseif path.sat_model == :YA15
		# use the Yenier & Atkinson (2015) finite fault factor
		return finite_fault_factor_ya15(m)
	elseif path.sat_model == :None
		return 0.0
	elseif path.sat_model == :Constant
		return path.heff
	else
		return NaN
	end
end

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
