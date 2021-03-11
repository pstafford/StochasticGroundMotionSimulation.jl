module StochasticGroundMotionSimulation

export Oscillator,
        period,
        transfer,
        transfer!,
        squared_transfer,
        squared_transfer!,
        site_amplification

# Write your package code here.
include("oscillator/PJSoscillator.jl")
include("fourier/PJSsite.jl")

end
