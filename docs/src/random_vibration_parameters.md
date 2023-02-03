```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Random Vibration Theory Parameters
Definition of the custom type `RandomVibrationParameters` to represent the model components/approaches used for the random vibration theory calculations.
In particular, the type stores symbols to define the:
- `excitation_duration` model to use,
- `rms_duration` model to use, and
- `peak_factor` model to use.

```@docs
  RandomVibrationParameters
```

Note that the default specification is:
```@example
RandomVibrationParameters() = RandomVibrationParameters(:DK80, :BT14, :BT15, :ACR)
```
However, an alternative constructor exists that takes a `pf_method` as an argument.
For this constructor, the `rms_duration` model is linked to the peak factor method:
- `DK80` is paired with `:BT15` as a default
- `CL56` is paired with `:BT12` as a default

As these are currently the only two `rms_duration` models implemented, the constructor is specified as:
```@example
RandomVibrationParameters(pf) = RandomVibrationParameters(pf, :BT14, ((pf == :DK80) ? :BT15 : :BT12), :ACR)
```

In all cases, the Boore & Thompson (2014) excitation duration model is employed as that is the only model currently implemented.


## Functionality

The overall goal of these random vibration methods is to compute:
```math
S_a = \psi \sqrt{ \frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment computed from `spectral_moment`, and ``D_{rms}`` is the root-mean-square duration computed from `rms_duration`.

The main methods used to interact with `RandomVibrationParameters` are:

```@docs
spectral_moment
spectral_moments
excitation_duration
rms_duration
peak_factor
rvt_response_spectral_ordinate
rvt_response_spectrum
```
