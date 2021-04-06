```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Fourier Parameters
Definition of the various custom types within the `StochasticGroundMotionSimulation` module.
Types to store properties related to source, path, and site components of the Fourier spectral model are provided.

```@docs
FourierParameters
```

## Source Parameters
The type `SourceParameters` holds the properties required to define the source spectrum of the Fourier Amplitude Spectrum.

```@docs
SourceParameters
```

## Path Parameters
The type `PathParameters` holds the properties required to define the path scaling of the Fourier Amplitude Spectrum.

```@docs
PathParameters
```

This type also hold instances of three other custom types that define aspects of the path scaling:
- `GeometricSpreadingParameters` defined in [Geometric Spreading](@ref)
- `NearSourceSaturationParameters` defined in [Near Source Saturation](@ref)
- `AnelasticAttenuationParameters` defined in [Anelastic Attenuation](@ref)

### Geometric Spreading

```@docs
GeometricSpreadingParameters
```

### Near-Source Saturation

Near source saturation models are represented within the `NearSourceSaturationParameters` type.
This type can simply identify existing models that are implemented, such as:
- Yenier & Atkinson (2015)
- Boore & Thompson (2015)
- Chiou & Youngs (2014) (the average of their ``h(\bm{M})`` term over all periods)

But, specific fixed values can also be provided as well as parameters that are subsequently operated upon:

```@docs
NearSourceSaturationParameters
```

Consider defining a new saturation model that was a simply bilinear model in ``\ln h(\bm{M})-\bm{M}`` space.

We simply pass in the various parameters that would be required for our saturation model into the available fields of `NearSourceSaturationParameters`, and then define a custom function that operates upon these fields.

```@example
m_min = 3.0
h_min = 0.5
m_hinge = 6.0
h_hinge = 5.0
m_max = 8.0
h_max = 30.0

sat = NearSourceSaturationParameters([m_min, m_hinge, m_max], [h_min, h_hinge, h_max])

function bilinear_saturation(m, sat)
  if m <= sat.mRefi[1]
    return sat.hconi[1]
  elseif m <= sat.mRefi[2]
    return sat.hconi[1] + (m - sat.mRefi[1])/(sat.mRefi[2]-sat.mRefi[1])*(sat.hconi[2] - sat.hconi[1])
  elseif m <= sat.mRefi[3]
    return sat.hconi[2] + (m - sat.mRefi[2])/(sat.mRefi[3]-sat.mRefi[2])*(sat.hconi[3] - sat.hconi[2])
  else
    return sat.hconi[3]
  end
end

```

Any subsequent calculation for a particular magnitude could then make use of this function along with a new `NearSourceSaturationParameters` instance that just contains a fixed saturation length.
```@example
m = 5.0
h_m = bilinear_saturation(m, sat)
new_sat = NearSourceSaturationParameters(h_m)
```


### Anelastic Attenuation

```@docs
AnelasticAttenuationParameters
```

## Site Parameters

The type `SiteParameters` holds information related to the site response -- both impedance effects and damping.

```@docs
SiteParameters
```
