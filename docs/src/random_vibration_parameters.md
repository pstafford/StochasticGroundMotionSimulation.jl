```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Random Vibration Theory Parameters
Definition of the custom type `RandomVibrationParameters` to represent the model components/approaches used for the random vibration theory calculations.
In particular, the type stores symbols to define the:
- `pf_method` the peak factor model/method to use.
- `dur_ex` specifies the excitation duration model to use,
- `dur_rms` specifies the model to use for converting excitation to RMS duration, and
- `dur_region`: the region or tectonic setting for excitation duration predictions, _i.e._ `:ACR` or `:SCR` for active and stable crustal regions

```@docs
  RandomVibrationParameters
```

Note that the default specification is:
```@example
RandomVibrationParameters() = RandomVibrationParameters(:DK80, :BT14, :BT15, :ACR)
```
However, an alternative constructor exists that takes a `pf_method` as a single argument, or that take a `pf_method` and `dur_region` specification.
For these constructors, the `dur_rms` model is linked to the `pf_method` peak factor method:
- `DK80` is paired with `:BT15`, and is the default
- `CL56` is paired with `:BT12` 

As these are currently the only two `dur_rms` models implemented, the constructor is specified as:
```@example
RandomVibrationParameters(pf) = RandomVibrationParameters(pf, :BT14, ((pf == :DK80) ? :BT15 : :BT12), :ACR)
```

In all cases, the Boore & Thompson (2014, 2015) excitation duration model is employed as that is the only model currently implemented.
This model uses a standard source duration definition, related to the source corner frequencies, and then applies a path duration model that is purely a function of the equivalent point-source distance.
The path duration model depends upon the region, as in active crustal regions or stable crustal regions.


## Functionality

The overall goal of these random vibration methods is to compute:
```math
S_a = \psi \sqrt{ \frac{m_0}{D_{rms}}}
```
where ``\psi`` is the peak factor computed from `peak_factor`, ``m_0`` is the zeroth order spectral moment computed from `spectral_moment`, and ``D_{rms}`` is the root-mean-square duration computed from `dur_rms`.

The main methods used to interact with `RandomVibrationParameters` are:

```@docs
spectral_moment
spectral_moments
spectral_moments_gk
excitation_duration
rms_duration
peak_factor
rvt_response_spectral_ordinate
rvt_response_spectrum
rvt_response_spectrum!
```
