```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Site Scaling

Site response is defined in terms of site amplification, or impedance effects, as well as damping -- via a _kappa_ filter.

The overall site model is therefore written as:
`` S(f; \bm{\theta}_S) = S_I(f) \times S_K(f) ``
with ``S_I(f)`` representing the impedance effects, and ``S_K(f)`` being the kappa filter.

## Impedance functions

Currently, three impedance functions are implemented:
- `:Boore2016`: is the Boore (2016) impedance function for a Western US generic rock profile with ``V_{S,30}=760`` m/s
- `:AlAtik2021_cy14`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Chiou & Youngs (2014) GMM. The reference profile has a ``V_{S,30}=760`` m/s
- `:Unit`: simply provides a unit impedance for all frequencies, _i.e._, ``S_I(f)=1.0``

```@docs
site_amplification
```

## Kappa filter

In addition to the impedance effects, the near surface damping is represented by a generic kappa filter:
```math
S_K(f) = \exp\left( -\pi \kappa_0 f \right)
```

```@docs
kappa_filter
```
