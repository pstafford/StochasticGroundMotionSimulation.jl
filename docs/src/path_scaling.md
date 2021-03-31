```@meta
CurrentModule = StochasticGroundMotionSimulation
```

# Path Scaling
The path scaling can be broken into geometric spreading -- including the effects of near-source saturation -- and anelastic attenuation.
Within `StochasticGroundMotionSimulation` the `PathParameters` type holds custom structs that relate to each of these three components:
- `GeometricSpreadingParameters` defines the geometric spreading model (spreading rates, transition distances, and functional scaling)
- `NearSourceSaturationParameters` defines the near-source saturation model, or the finite fault factors
- `AnelasticAttenuationParameters` defines the properties of the anelastic attenuation model.
