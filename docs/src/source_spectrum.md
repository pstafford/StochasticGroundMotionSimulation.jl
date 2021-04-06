```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Fourier Source Spectrum

The source spectrum can be computed through interaction with the `SourceParameters` type, or the higher level `FourierParameters` type that holds a `SourceParameters` instance as a property.

The source spectrum ``E(f; \bm{\theta}E)`` is most commonly written in terms of:
```math
  E(f; \bm{\theta}_E) = \mathcal{C} M_0 E_s(f; \bm{\theta}_E)
```
where ``\mathcal{C}`` is a constant term, to be defined shortly, ``M_0`` is the seismic moment, and ``E_s(f; \bm{\theta}_E)`` is the source spectral shape.

The most commonly adopted source spectral shape is the ``\omega^2`` model that has the form:
```math
 E_s(f) = \frac{1}{1 + \left(\frac{f}{f_c}\right)^2}
```

## Functionality

```@docs
fourier_constant
fourier_source_shape
fourier_source
corner_frequency
magnitude_to_moment
```
