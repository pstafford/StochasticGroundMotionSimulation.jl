```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Site Scaling

Site response is defined in terms of site amplification, or impedance effects, as well as damping -- via a _kappa_ filter.

The overall site model is therefore written as:
`` S(f; \bm{\theta}_S) = S_I(f) \times S_K(f) ``
with ``S_I(f)`` representing the impedance effects, and ``S_K(f)`` being the kappa filter.

## Impedance functions

Impedance effects are represented using custom types that define an `amplification` function. 
These functions take a frequency as an argument and return the corresponding amplification level.
Each of these custom types is a subtype of the abstract type `SiteAmplification`

Currently, the following impedance functions (custom types) are implemented:
- `SiteAmpBoore2016_760`: is the Boore (2016) impedance function for a Western US generic rock profile with ``V_{S,30}=760`` m/s
- `SiteAmpAlAtikAbrahamson2021_ask14_620`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Abrahamson _et al._ (2014) GMM. The reference profile has a ``V_{S,30}=620`` m/s
- `SiteAmpAlAtikAbrahamson2021_ask14_760`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Abrahamson _et al._ (2014) GMM. The reference profile has a ``V_{S,30}=760`` m/s
- `SiteAmpAlAtikAbrahamson2021_ask14_1100`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Abrahamson _et al._ (2014) GMM. The reference profile has a ``V_{S,30}=1100`` m/s
- `SiteAmpAlAtikAbrahamson2021_bssa14_620`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Boore _et al._ (2014) GMM. The reference profile has a ``V_{S,30}=620`` m/s
- `SiteAmpAlAtikAbrahamson2021_bssa14_760`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Boore _et al._ (2014) GMM. The reference profile has a ``V_{S,30}=760`` m/s
- `SiteAmpAlAtikAbrahamson2021_bssa14_1100`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Boore _et al._ (2014) GMM. The reference profile has a ``V_{S,30}=1100`` m/s
- `SiteAmpAlAtikAbrahamson2021_cb14_620`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Campbell & Bozorgnia (2014) GMM. The reference profile has a ``V_{S,30}=620`` m/s
- `SiteAmpAlAtikAbrahamson2021_cb14_760`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Campbell & Bozorgnia (2014) GMM. The reference profile has a ``V_{S,30}=760`` m/s
- `SiteAmpAlAtikAbrahamson2021_cb14_1100`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Campbell & Bozorgnia (2014) GMM. The reference profile has a ``V_{S,30}=1100`` m/s
- `SiteAmpAlAtikAbrahamson2021_cy14_620`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Chiou & Youngs (2014) GMM. The reference profile has a ``V_{S,30}=620`` m/s
- `SiteAmpAlAtikAbrahamson2021_cy14_760`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Chiou & Youngs (2014) GMM. The reference profile has a ``V_{S,30}=760`` m/s
- `SiteAmpAlAtikAbrahamson2021_cy14_1100`: is the Al Atik & Abrahamson (2021) impedance function obtained by inverting the Chiou & Youngs (2014) GMM. The reference profile has a ``V_{S,30}=1100`` m/s
- `SiteAmpUnit`: simply provides a unit impedance for all frequencies, _i.e._, ``S_I(f)=1.0``

```@docs
site_amplification
SiteAmpUnit
SiteAmpBoore2016_760
SiteAmpAlAtikAbrahamson2021_ask14_620
SiteAmpAlAtikAbrahamson2021_ask14_760
SiteAmpAlAtikAbrahamson2021_ask14_1100
SiteAmpAlAtikAbrahamson2021_bssa14_620
SiteAmpAlAtikAbrahamson2021_bssa14_760
SiteAmpAlAtikAbrahamson2021_bssa14_1100
SiteAmpAlAtikAbrahamson2021_cb14_620
SiteAmpAlAtikAbrahamson2021_cb14_760
SiteAmpAlAtikAbrahamson2021_cb14_1100
SiteAmpAlAtikAbrahamson2021_cy14_620
SiteAmpAlAtikAbrahamson2021_cy14_760
SiteAmpAlAtikAbrahamson2021_cy14_1100
```

## Kappa filter

In addition to the impedance effects, the near surface damping is represented by a generic kappa filter:
```math
S_K(f) = \exp\left( -\pi \kappa_0 f \right)
```

```@docs
kappa_filter
```
