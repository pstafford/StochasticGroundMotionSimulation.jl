```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# StochasticGroundMotionSimulation

Documentation for the Julia package `StochasticGroundMotionSimulation.jl`.
The main module `StochasticGroundMotionSimulation` provides an interface to the stochastic method for the simulation of response spectral ordinates via Random Vibration Theory.

## Example
A number of default parameters are already set within the package.

```@example
using StochasticGroundMotionSimulation

Δσ = 100.0
Rrefi = [ 1.0, 50.0, Inf ]
γi = [ 1.0, 0.5 ]
Q0 = 200.0
η = 0.5
κ0 = 0.039

src = SourceParameters(Δσ)
geo = GeometricSpreadingParameters(Rrefi, γi, :CY14)
sat = NearSourceSaturationParameters(:BT15)
ane = AnelasticAttenuationParameters(Q0, η)
path = PathParameters(geo, sat, ane)
site = SiteParameters(κ0, :Boore2016)

fas = FourierParameters(src, path, site)

rvt = RandomVibrationParameters()

T = 1.0
m = 6.0
r_rup = 10.0
r_ps = equivalent_point_source_distance(r_rup, m, fas)

Sa = rvt_response_spectral_ordinate(T, m, r_ps, fas, rvt)
println("Sa = $(Sa)g, for m = $m, r_rup = $r_rup, and a period of T = $(T) s")

```


```@contents
Pages = ["fourier_parameters.md","random_vibration_parameters.md","sdof_parameters.md]
Depth = 4
```

```@index
```

```@autodocs
Modules = [StochasticGroundMotionSimulation]
```
