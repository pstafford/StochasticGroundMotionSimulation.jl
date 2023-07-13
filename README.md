# StochasticGroundMotionSimulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pstafford.github.io/StochasticGroundMotionSimulation.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pstafford.github.io/StochasticGroundMotionSimulation.jl/dev)
[![Build Status](https://github.com/pstafford/StochasticGroundMotionSimulation.jl/workflows/CI/badge.svg)](https://github.com/pstafford/StochasticGroundMotionSimulation.jl/actions)
[![codecov](https://codecov.io/gh/pstafford/StochasticGroundMotionSimulation.jl/branch/master/graph/badge.svg?token=EDEF06FN61)](https://codecov.io/gh/pstafford/StochasticGroundMotionSimulation.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4667333.svg)](https://doi.org/10.5281/zenodo.4667333)

[Julia](http://www.julialang.org) package to simulate response spectral ordinates via random vibration theory.
The package also provides general functionality for working with Fourier amplitude spectra and duration models.

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

## Installation

First, a working version of [Julia](http://www.julialang.org) needs to be installed.
The relevant binary (or source code) can be downloaded from the [Julia Downloads Page](https://julialang.org/downloads/).

`StochasticGroundMotionSimulation.jl` is a registered package and can be installed directly from the package manager.

Within a Julia REPL session, access the package manager via `]`, and then at the `pkg>` prompt type (below the `pkg>` component is part of the prompt, so only the `add ...` portion is necessary).
```julia
pkg> add StochasticGroundMotionSimulation
```

## Usage

Within a Julia session, bring the functionality of `StochasticGroundMotionSimulation.jl` into scope by typing (here the `julia>` component represents the prompt within a REPL session, within a text editor, simply type `using StochasticGroundMotionSimulation`):
```julia
julia> using StochasticGroundMotionSimulation
```

## Accessing Help

Aside from the [documentation](https://pstafford.github.io/StochasticGroundMotionSimulation.jl/stable) accessible from the links at the top of this page, descriptions of methods and types within the package can be accessed within REPL sessions (or within Juno).

Within, a REPL session, enter `?` to access the help prompt.
Then, type the relevant item.
For example:

```julia
help?> FourierParameters
```

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
