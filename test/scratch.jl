using StochasticGroundMotionSimulation
using Plots
using QuadGK
using BenchmarkTools
using Test
# using ForwardDiff
using Polynomials
using PolynomialRoots

fi = exp10.(range(-2.0, stop=2.0, length=2001))

m = 6.0
r = 100.0
fas = FASParams(100.0, [1.0, 50.0, Inf], [1.0, 0.5], 200.0, 0.4, 0.039)
sdof = Oscillator(1.0)

Afi = map(f -> fourier_spectral_ordinate(f, m, r, fas), fi)
Hfi = transfer(fi, sdof)

fac = 1.5

plot(fi, Afi.^2, lab="Af^2", box=:true, grid=:on, xscale=:log10)
plot!(fi, (Hfi .* Afi).^2, lab="Hf^2Af^2")
vline!([sdof.f_n, sdof.f_n/fac, sdof.f_n*fac])
xlims!(sdof.f_n/ceil(fac),sdof.f_n*ceil(fac))

fun(f) = squared_transfer(f, sdof)

fac = 1.5

fj = collect(range(sdof.f_n/fac, stop=sdof.f_n*fac, length=101))
H2j = map(f -> fun(f), fj)

@time igk = quadgk(fun, sdof.f_n/fac, sdof.f_n*fac)[1]
@time isr = simpsons_rule(fj, H2j)
@time itr = trapezoidal_rule(fj, H2j)

@test igk ≈ isr rtol=1e-5
@test igk ≈ itr rtol=1e-5


fun(f) = squared_transfer(f, sdof) * ( fourier_spectral_ordinate(f, m, r, fas)^2 )

fj = collect(range(sdof.f_n/fac, stop=sdof.f_n*fac, length=101))
H2j = map(f -> fun(f), fj)

@time igk = quadgk(fun, sdof.f_n/fac, sdof.f_n*fac)[1]
@time isr = simpsons_rule(fj, H2j)
@time itr = trapezoidal_rule(fj, H2j)

@test igk ≈ isr rtol=1e-5
@test igk ≈ itr rtol=1e-4

plot(fi, Afi.^2, lab="Af^2", box=:true, grid=:on)
plot!(fi, (Hfi .* Afi).^2, lab="Hf^2Af^2")
# vline!([sdof.f_n, sdof.f_n/fac, sdof.f_n*fac])
xlims!(1e-3, sdof.f_n/fac)

fj = collect(range(1e-3, stop=sdof.f_n/fac, length=101))
H2j = map(f -> fun(f), fj)

@time igk = quadgk(fun, 1e-3, sdof.f_n/fac)[1]
@time isr = simpsons_rule(fj, H2j)
@time itr = trapezoidal_rule(fj, H2j)

@test igk ≈ isr rtol=1e-5
@test igk ≈ itr rtol=1e-3

fj = collect(range(sdof.f_n*fac, stop=1e1, length=101))
H2j = map(f -> fun(f), fj)

@time igk = quadgk(fun, sdof.f_n*fac, 1e1)[1]
@time isr = simpsons_rule(fj, H2j)
@time itr = trapezoidal_rule(fj, H2j)

@test igk ≈ isr rtol=1e-3
@test igk ≈ itr rtol=1e-3



fun(f) = (2*π*f)^2 * squared_transfer(f, sdof) * ( fourier_spectral_ordinate(f, m, r, fas)^2 )

fj = collect(range(sdof.f_n/fac, stop=sdof.f_n*fac, length=101))
H2j = map(f -> fun(f), fj)
fk = collect(range(sdof.f_n/fac, stop=sdof.f_n*fac, length=71))
H2k = map(f -> fun(f), fk)

@time igk = quadgk(fun, sdof.f_n/fac, sdof.f_n*fac)[1]
@time isr = simpsons_rule(fj, H2j)
@time itr = trapezoidal_rule(fk, H2k)

@test igk ≈ isr rtol=1e-5
@test igk ≈ itr rtol=1e-4


function fourier_derivative(f, m, r, fas)
    βQ = π * r / ( fas.Q0 * fas.cQ )
    βκ = π * fas.κ0
    fc, ~, ~ = corner_frequency(m, fas)

    return -βκ*f^4 - βQ*(1.0-fas.η)*f^(4-fas.η) + 2*(1.0-fc^2)*f^3 - βκ*fc^2*f^2 - βQ*(1.0-fas.η)*fc^2*f^(2-fas.η) + 2*fc^2*f
end

fas = FASParams(100.0, [1.0, 50.0, Inf], [1.0, 0.5], 200.0, 0.0, 0.039)

fi = exp10.(range(-2, stop=2, length=101))
Afi = fourier_spectrum(fi, m, r, fas)
dAfi = map(f -> fourier_derivative(f, m, r, fas), fi)

# FAS(f) = fourier_spectral_ordinate(f, m, r, fas)
# dFAS(f) = ForwardDiff.derivative(FAS, f)
#
# gAfi = map(f -> dFAS(f), fi)

βQ = π * r / ( fas.Q0 * fas.cQ )
βκ = π * fas.κ0
fc, ~, ~ = corner_frequency(m, fas)

dAfp = Polynomial([ 0, 2*fc^2, -fc^2*(βQ+βκ), 2*(1.0-fc^2), -(βQ+βκ)])
ri = PolynomialRoots.roots([ 0, 2*fc^2, -fc^2*(βQ+βκ), 2*(1.0-fc^2), -(βQ+βκ)])

real.(ri)
imag.(ri)


plot(fi, Afi, lab="FAS")
vline!([abs(ri[1])])


plot(fi, dAfi, lab="dFAS")
ylims!(-1e2,1e2)
