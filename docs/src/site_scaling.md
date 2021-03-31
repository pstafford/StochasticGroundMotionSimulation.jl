```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Site Scaling
Site response is defined in terms of site amplification, or impedance effects, as well as damping -- via a _kappa_ filter.

The overall site model is therefore written as:
`` S(f; \bm{\theta}_S) = S_I(f) \times S_K(f) ``
with ``S_I(f)`` representing the impedance effects, and ``S_K(f)`` being the kappa filter.
