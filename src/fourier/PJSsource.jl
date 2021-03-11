

"""
	magnitude_to_moment(m::Real)

Converts moment magnitude to seismic moment (in dyne-cm).

# Examples
```julia-repl
	m = 6.0
	M0 = magnitude_to_moment(m)
```
"""
function magnitude_to_moment(m::Real)
	# equivalent to: return 10.0^( 1.5*( m + 10.7 ) )
    return 10.0^( 1.5m + 16.05 )
end

"""
	corner_frequency_brune(m::Real, Δσ::Real, β::Real=3.5)

Computes the corner frequency using the Brune model.
	- `m` is the moment magnitude
	- `Δσ` is the stress drop in units of bars
	- `β` is the shear-wave velocity at the source in units of km/s

# Examples
```julia-repl
	m = 6.0
	Δσ = 100.0
	β = 3.5
	fc = corner_frequency_brune(m, Δσ, β)
```
"""
function corner_frequency_brune(m::Real, Δσ::Real, β::Real=3.5)
    Mo = magnitude_to_moment(m)
    return 4.9058e6 * β * ( Δσ / Mo )^(1/3)
end


"""
    corner_frequency_atkinson_silva_2000(m::Real)

Computes the corner frequencies, `fa`, `fb`, and the mixing parameter `ε` from the Atkinson & Silva (2000) double corner frequency model.
This is the default source corner frequency model used by Boore & Thompson (2014) to define their source duration. But note that they just use fa. This function returns fb and ε also

# Examples
```julia-repl
    m = 6.0
    fa, fb, ε = corner_frequency_atkinson_silva_2000(m)
```
"""
function corner_frequency_atkinson_silva_2000(m::Real)
	fa = 10.0^( 2.181 - 0.496*m )
	fb = 10.0^( 2.410 - 0.408*m )
	ε = 10.0^( 0.605 - 0.255*m )
    return fa, fb, ε
end


"""
    corner_frequency(m::Real, fas::Union{PJSfasPararms,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::String="Brune")

Computes a 3-tuple of corner frequency components depending upon source spectrum type. By default the single-corner Brune spectrum is considered, but it `fc_fun` equals `"Atkinson_Silva_2000"` then the components of the double-corner spectrum are returned. If some other string is passed then a 3-tuple of NaN::Real values is returned.

# Examples
```julia-repl
    m = 6.0
    Δσ = 100.0
    κ0 = 0.035
    fas = FASParams(Δσ, κ0)
    # compute single corner frequency
    fc, tmp1, tmp2 = corner_frequency(m, fas)
    fc, tmp1, tmp2 = corner_frequency(m, fas; fc_fun="Brune")
    # compute double corner frequencies
    fa, fb, ε = corner_frequency(m, fas; fc_fun="Atkinson_Silva_2000")
```
"""
function corner_frequency(m::Real, fas::Union{FASParams,FASParamsGeo,FASParamsQr,FASParamsGeoQr}; fc_fun::Symbol=:Brune)
    if fc_fun == :Brune
        return corner_frequency_brune(m, fas.Δσ, fas.β)::Real, NaN::Real, NaN::Real
    elseif fc_fun == :Atkinson_Silva_2000
        return corner_frequency_atkinson_silva_2000(m)
    else
        return NaN::Real, NaN::Real, NaN::Real
    end
end
