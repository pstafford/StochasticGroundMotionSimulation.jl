
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
    # Δx = x[2] - x[1]
    # return Δx * ( sum(y) - (y[1] + y[end])/2 )
    return (x[2] - x[1]) * ( sum(y) - (y[1] + y[end])/2 )
end
