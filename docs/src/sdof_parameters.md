```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Single Degree of Freedom Oscillator Parameters
Definition of the custom type, `Oscillator` to represent a single degree of freedom (SDOF) oscillator.
Type simply stores the oscillator frequency and damping ratio.

```@docs
  Oscillator
```

## Functionality

The main methods that are used to interact with `Oscillator` instances are:

```@docs
  period
  transfer
  squared_transfer
```
