```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Fourier Parameters
Definition of the various custom types within the `StochasticGroundMotionSimulation` module.
Types to store properties related to source, path, and site components of the Fourier spectral model are provided.

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

### Near Source Saturation

```docs
NearSourceSaturationParameters
```

### Anelastic Attenuation

```docs
AnelasticAttenuationParameters
```

## Site Parameters
The type `SiteParameters` holds information related to the site response -- both impedance effects and damping.

```docs
SiteParameters
```
