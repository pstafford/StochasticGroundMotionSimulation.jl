```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# StochasticGroundMotionSimulation

Documentation for the Julia package `StochasticGroundMotionSimulation.jl`.
The main module `StochasticGroundMotionSimulation` provides an interface to the stochastic method for the simulation of response spectral ordinates via Random Vibration Theory.

The package makes use of three main components:
- `FourierParameters`: defining the properties of the Fourier amplitude spectrum  
- `RandomVibrationParameters`: defining the properties of the random vibration theory calculations. Specifically, defining the duration model(s) to be used along with the peak factor method
- `Oscillator`: defines the properties of the single degree-of-freedom oscillator for which response spectral ordinates are computed.

To compute a response spectrum, or a response spectral ordinate, the above components along with a definition of a magnitude-distance scenario are required.

The package is written to enable automatic differentiation operations to be applied to the principle parameters defining the Fourier amplitude spectrum.
This is done in order to facilitate the use of this package for gradient-based inversions of observed ground-motions, or inversions of published empirical ground-motion models.

## Contents
```@contents
Pages = ["fourier_parameters.md","random_vibration_parameters.md","sdof_parameters.md]
Depth = 4
```

## Example
A number of default parameters are already set within the package.

```@example
using StochasticGroundMotionSimulation

# specify some parameters defining the Fourier amplitude spectrum
Δσ = 100.0                  # the stress parameter (in bar)
Rrefi = [ 1.0, 50.0, Inf ]  # reference distances for the geometric spreading
γi = [ 1.0, 0.5 ]           # geometric spreading rates for distances between the references distances
Q0 = 200.0                  # quality factor ``Q_0 \in Q(f) = Q_0 f^\eta``
η = 0.5                     # quality exponent  ``\eta \in Q(f) = Q_0 f^\eta``
κ0 = 0.039                  # site kappa value

# construct a `SourceParameters` instance
src = SourceParameters(Δσ)
# construct a `GeometricSpreadingParameters` instance using the reference distances, spreading rates, and spreading model
geo = GeometricSpreadingParameters(Rrefi, γi, :CY14)
# define the near-source saturation model
sat = NearSourceSaturationParameters(:BT15)
# define the anelastic attenuation properties
ane = AnelasticAttenuationParameters(Q0, η)
# use the `geo`, `sat` and `ane` instances to construct a `PathParameters` instance
path = PathParameters(geo, sat, ane)
# define the `SiteParameters`
site = SiteParameters(κ0, SiteAmpBoore2016_760())

# combine `src`, `path` and `site` instances to define the overall `FourierParameters`
fas = FourierParameters(src, path, site)

# use default properties for the `RandomVibrationParameters`
rvt = RandomVibrationParameters()

# define the response period, and magnitude-distance scenario of interest
T = 1.0
m = 6.0
r_rup = 10.0
# compute the equivalent point-source distance for this scenario
r_ps = equivalent_point_source_distance(r_rup, m, fas)

# compute the response spectral ordinate
Sa = rvt_response_spectral_ordinate(T, m, r_ps, fas, rvt)

# write out the results
println("Sa = $(round(Sa, sigdigits=4)) g, for m = $m, r_rup = $r_rup km, and a period of T = $(T) s")

```
