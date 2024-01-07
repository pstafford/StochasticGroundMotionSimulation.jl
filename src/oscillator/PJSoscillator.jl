
"""
    Oscillator{T<:Float64}

Custom type to represent a SDOF oscillator. The type has two fields:
  - `f_n` is the natural frequency of the oscillator
  - `ζ_n` is the damping ratio

# Examples
```julia-repl
    sdof = Oscillator( 1.0, 0.05 )
```
"""
struct Oscillator{T<:Float64}
    f_n::T
    ζ_n::T
end


"""
    Oscillator(f_n)

Default initializer setting the damping ratio to 5% of critical
  - `f_n` is the natural frequency of the oscillator (in Hz)

# Examples
```julia-repl
    f_n = 2.0
    sdof = Oscillator(f_n)
```
"""
Oscillator(f_n::T) where {T<:Float64} = Oscillator(f_n, 0.05)


"""
	period(sdof::Oscillator)

Natural period (s) of the sdof Oscillator.
"""
function period(sdof::Oscillator)
	return 1.0/sdof.f_n
end

@doc raw"""
    transfer(f::T, sdof::Oscillator) where {T<:Real}

Compute the modulus of the transfer function for a SDOF system.

The transfer function is defined as:
```math
|H(f,f_n,\zeta_n)| = \frac{1}{\sqrt{ \left(1 - \beta^2 \right)^2 + \left(2\zeta_n\beta\right)^2 }}
```
where ``\beta`` is the tuning ratio defined by ``f/f_n``.

# Examples
```julia-repl
    f = 2.0
    sdof = Oscillator(1.0, 0.05)
    Hf = transfer(f, sdof)
```

See also: [`squared_transfer`](@ref)
"""
function transfer(f, sdof::Oscillator)
    # tuning ratio
    β = f / sdof.f_n
    return 1.0 / sqrt( (1.0 - β^2)^2 + (2sdof.ζ_n*β)^2 )
end


"""
	squared_transfer(f, sdof::Oscillator)

Compute the square of the transfer function for a SDOF system, `sdof`, at frequency `f`.

# Examples
```julia-repl
	f = 2.0
	# create sdof with natural frequency f_n=1.0 and damping ζ=0.05
	sdof = Oscillator( 1.0, 0.05 )
	Hf2 = squared_transfer( f, sdof )
```

See also: [`transfer`](@ref)
"""
function squared_transfer(f, sdof::Oscillator)
    # tuning ratio
    β = f / sdof.f_n
    return 1.0 / ( (1.0 - β^2)^2 + (2sdof.ζ_n*β)^2 )
end


"""
    transfer(f::Vector{T}, sdof::Oscillator) where T<:Real

Computes the modulus of the transfer function of a SDOF for a vector of frequencies
  - `f::Vector` is the vector of frequencies
  - `sdof::Oscillator` is the oscillator instance

# Examples
```julia-repl
  f = collect(range(0.1, stop=10.0, step=0.01))
  sdof = Oscillator(1.0)
  Hf = transfer(f, sdof)
```
"""
function transfer(f::Vector{T}, sdof::Oscillator) where T<:Real
    # tuning ratio
    Hf = similar(f)
    for (i, fi) in pairs(f)
    # for i in 1:lastindex(f)
        @inbounds Hf[i] = transfer(fi, sdof)
    end
    return Hf
end


"""
    transfer!(Hf::Vector{T}, f::Vector{T}, sdof::Oscillator) where T<:Real

Computes the modulus of the transfer function of a SDOF for a vector of frequencies in place
	- `Hf::Vector` is the pre-allocated vector into which the results are stored
    - `f::Vector` is the vector of frequencies
    - `sdof::Oscillator` is the oscillator instance

# Examples
```julia-repl
    f = collect(range(0.1, stop=10.0, step=0.01))
    sdof = Oscillator(1.0)
    Hf = similar(f)
    transfer!(Hf, f, sdof)
```
"""
function transfer!(Hf::Vector{T}, f::Vector{T}, sdof::Oscillator) where T<:Real
    for (i, fi) in pairs(f)
    # for i in 1:lastindex(f)
        @inbounds Hf[i] = transfer(fi, sdof)
	end
	return nothing
end

"""
    squared_transfer!(Hf2::Vector{T}, f::Vector{T}, sdof::Oscillator) where T<:Real

Computes the square of the modulus of the transfer function of a SDOF for a vector of frequencies in place:
	- `Hf2::Vector` is the pre-allocated vector into which the results are stored
    - `f::Vector` is the vector of frequencies
    - `sdof::Oscillator` is the oscillator instance
Inputs derive from the `Real` type and so are differentiable.

# Examples
```julia-repl
    f = collect(range(0.1, stop=10.0, step=0.01))
    sdof = Oscillator(1.0)
    Hf2 = similar(f)
    squared_transfer!(Hf2, f, sdof)
```
"""
function squared_transfer!(Hf2::Vector{T}, f::Vector{U}, sdof::Oscillator) where {T<:Real,U<:Real}
    for (i, fi) in pairs(f)
    # for i in 1:lastindex(f)
        @inbounds Hf2[i] = squared_transfer(fi, sdof)
	end
	return nothing
end

