
"""
    simpsons_rule(x::Vector, y::Vector)

Numerical integration via Simpson's rule.

Lengths of vectors must be equal to each other, and be odd (to have an even number of intervals).
Values in `x` represent abcissa values, while `y` are the integrand counterparts.
"""
function simpsons_rule(x::Vector, y::Vector)
    n = length(y)-1
    n % 2 == 0 || error("`y` length (number of intervals) must be odd")
    length(x)-1 == n || error("`x` and `y` length must be equal")
    h = (x[end]-x[1])/n
    @inbounds @views s = sum(y[1:2:n] .+ 4y[2:2:n] .+ y[3:2:n+1])
    return h/3 * s
end


"""
    trapezoidal_rule(x::Vector, y::Vector)

Numerical integration via the Trapezoidal rule.

Lengths of vectors must be equal to each other.
Values in `x` represent equally-spaced abcissa values, while `y` are the integrand counterparts.
"""
function trapezoidal_rule(x::Vector, y::Vector)
    return (x[2] - x[1]) * ( sum(y) - (y[1] + y[end])/2 )
end

function trapezoidal(fun::Function, n, flim...)
    ii = 0.0
    for i in 2:length(flim)
        xi = collect(range(flim[i-1], stop=flim[i], length=n))
        yi = fun.(xi)
        ii += trapezoidal_rule(xi, yi)
    end
    return ii
end


function gauss_interval(integrand::Function, n, fmin, fmax)
    xi, wi = gausslegendre(n)
    ifi = @. integrand( (fmax-fmin)/2 * xi + (fmin+fmax)/2 )
    return (fmax-fmin)/2 * dot( wi, ifi )
end

function gauss_intervals(fun::Function, n, flim...)
    xi, wi = gausslegendre(n)
    ii = 0.0
    for i in 2:length(flim)
        ii += (flim[i]-flim[i-1])/2 * dot( wi, fun.( (flim[i]-flim[i-1])/2 * xi .+ (flim[i]+flim[i-1])/2 ) )
    end
    return ii
end
