```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Fourier Amplitude Spectrum
The Fourier amplitude spectrum (FAS) can be represented as the product of source, path and site contributions.

Specifically, the Fourier amplitude spectrum ``|A(f)|`` of acceleration (in units of m/s) is defined as:

```math
|A(f; \bm{\theta})| = E(f; \bm{\theta}_E)\times P(f; \bm{\theta}_P) \times S(f; \bm{\theta}_S)
```

where ``f`` is a frequency in Hz, and ``\bm{\theta}`` holds all of the relevant model parameters and predictor variables.
The [Fourier Source Spectrum](@ref), ``E(f; \bm{\theta}_E)`` is a function of the earthquake magnitude ``m``, as well as other properties of the source.
The [Path Scaling](@ref), ``P(f; \bm{\theta}_P)`` accounts for the effects of both geometric spreading and anelastic attenuation.
The [Site Scaling](@ref), ``S(f; \bm{\theta}_S)`` includes the effects of near-surface impedance as well as damping (``\kappa_0``) effects.
