# StochasticGroundMotionSimulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pstafford.github.io/StochasticGroundMotionSimulation.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pstafford.github.io/StochasticGroundMotionSimulation.jl/dev)
[![Build Status](https://github.com/pstafford/StochasticGroundMotionSimulation.jl/workflows/CI/badge.svg)](https://github.com/pstafford/StochasticGroundMotionSimulation.jl/actions)
[![codecov](https://codecov.io/gh/pstafford/StochasticGroundMotionSimulation.jl/branch/master/graph/badge.svg?token=EDEF06FN61)](https://codecov.io/gh/pstafford/StochasticGroundMotionSimulation.jl)
[![DOI](https://zenodo.org/badge/338342369.svg)](https://zenodo.org/badge/latestdoi/338342369)

Julia package to simulate response spectral ordinates via random vibration theory.

Package defines new custom types:
- `FourierParameters`: representing the parameters of the Fourier amplitude spectrum
- `Oscillator`: representing a single degree-of-freedom oscillator, and
- `RandomVibrationParameters`: defining methods/models used for random vibration calculations

The `FourierParameters` type is constructed from three components:
- `SourceParameters`: representing source properties, such as stress parameter, source velocity and density, _etc_
- `PathParameters`: representing the path scaling. This component is itself comprised of three components:
  - `GeometricSpreadingParameters`: defines the geometric spreading model
  - `NearSourceSaturationParameters`: defines the near-source saturation model, and
  - `AnelasticAttenuationParameters`: defines the anelastic attenuation
- `SiteParameters`: defines both the site amplification/impedance function and the site damping (via a kappa filter)

The package is developed in a manner to enable automatic differentiation operations to be performed via `ForwardDiff.jl`.
This makes the package suitable for gradient-based inversions of ground-motion data, as well as inversions of published ground-motion models.


## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
