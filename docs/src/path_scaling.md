```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Path Scaling
The path scaling can be broken into geometric spreading -- including the effects of near-source saturation -- and anelastic attenuation.
Within `StochasticGroundMotionSimulation` the `PathParameters` type holds custom structs that relate to each of these three components:
- `GeometricSpreadingParameters` defines the geometric spreading model (spreading rates, transition distances, and functional scaling)
- `NearSourceSaturationParameters` defines the near-source saturation model, or the finite fault factors
- `AnelasticAttenuationParameters` defines the properties of the anelastic attenuation model.





## Geometric Spreading
Parameters for representing geometric spreading are contained within a `GeometricSpreadingParameters` instance.
To compute the actual geometric spreading for a given distance we make use of the `geometric_spreading` function:

```@docs
geometric_spreading
```

This function takes different options that define different spreading functions.
For example, the `:CY14` option uses the functional form of Chiou & Youngs (2014), but uses an equivalent point-source distance throughout.

```math
  \ln g(r_{ps}) = -\gamma_1 \ln(r_{ps}) + \frac{\left(\gamma_1 -\gamma_f\right)}{2} \ln\left( \frac{ r_{ps}^2 + r_t^2 }{r_{0}^2 + r_t^2} \right)
```

The alternative `:CY14mod` option combines a point-source distance in the near field, with `r_rup` scaling in the far field.
```math
\ln g(r_{ps},r_{rup}) = -\gamma_1 \ln(r_{ps}) + \frac{\left(\gamma_1 -\gamma_f\right)}{2} \ln\left( \frac{ r_{rup}^2 + r_t^2 }{r_{0}^2 + r_t^2} \right)
```

In both of the above cases, the ``r_{0}`` term is the reference distance that is used to define the source spectral amplitude.


## Near Source Saturation

The `NearSourceSaturationParameters` are crucial for computing the equivalent point-source distance metric.
Generally, the equivalent point-source distance can be computed via:
```math
  r_{ps} = \left( r_{rup}^n + h(\bm{M})^n \right)^{1/n}
```
and it is most common to follow Boore & Thompson (2015) and to use ``n=2`` so that:
```math
  r_{ps} = \sqrt{ r_{rup}^2 + h(\bm{M})^2 }
```


## Anelastic Attenuation

The anelastic attenuation filter has the general form:
```math
  \exp\left[ -\frac{\pi f r}{Q(f) c_Q} \right]
```
where, normally, ``Q(f)=Q_0 f^\eta`` such that:
```math
  \exp\left[ -\frac{\pi f^{1-\eta} r}{Q_0 c_Q} \right]
```

The `AnelasticAttenuationParameters` type therefore holds the values of ``Q0``, ``\eta``, and ``c_Q``.
In addition, it holds a field `rmetric` that can take values of `:Rrup` and `:Rps` depending upon whether one wishes to interpret the distance within the exponential function as the rupture distance, `:Rrup`, or the equivalent point-source distance, `:Rps`.
