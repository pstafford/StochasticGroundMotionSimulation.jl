# Custom Models

The package now supports custom FAS and duration models for maximum flexibility.

## Quick Start
```julia
using StochasticGroundMotionSimulation

# Define a custom FAS function
my_fas = (f, m, r, params) -> begin
    stress_drop, Q0, kappa = params
    # ... your model implementation
end

# Create a functional model
fas_model = FunctionalFASModel(my_fas, [100.0, 200.0, 0.04])

# Use it in RVT calculations 
Sa = rvt_response_spectral_ordinate(1.0, 6.0, 10.0, fas_model, duration_model)
```

## Custom Type Models

For more complex models, define your own types:
```julia
struct MyFASModel <: AbstractFASModel
    # your parameters
end

function compute_fas(model::MyFASModel, freq, mag, dist)
    # your implementation
end
```

## ForwardDiff Compatibility

All custom models must be ForwardDiff-compatible. Use the validation functions:
```julia
validate_fas_model(my_model)
validate_duration_model(my_model)
```

---

## API Reference

### Abstract Types
```@docs
StochasticGroundMotionSimulation.AbstractFASModel
StochasticGroundMotionSimulation.AbstractDurationModel
```

### FAS Model Types
```@docs
StochasticGroundMotionSimulation.CustomModels.FunctionalFASModel
StochasticGroundMotionSimulation.CustomModels.CustomFASModel
StochasticGroundMotionSimulation.CustomModels.HybridFASModel
```

### Duration Model Types
```@docs
StochasticGroundMotionSimulation.CustomModels.FunctionalDurationModel
StochasticGroundMotionSimulation.CustomModels.CustomDurationModel
```

### Core Interface Functions
```@docs
StochasticGroundMotionSimulation.CustomModels.compute_fas
StochasticGroundMotionSimulation.CustomModels.compute_duration
```

### Validation Functions
```@docs
StochasticGroundMotionSimulation.CustomModels.validate_fas_model
StochasticGroundMotionSimulation.CustomModels.validate_duration_model
```

### RVT Computation with Custom Models
```@docs
StochasticGroundMotionSimulation.CustomModels.rvt_response_spectral_ordinate_custom
```

### Utilities
```@docs
StochasticGroundMotionSimulation.CustomModels.FourierParametersWrapper
StochasticGroundMotionSimulation.CustomModels.ExistingDurationWrapper
```
