
module StochasticGroundMotionSimulation

using Interpolations
using Roots
using ForwardDiff
using ForwardDiff: Dual
using QuadGK
using FastGaussQuadrature
using LinearAlgebra


export Oscillator,
	FourierParameters,
	SourceParameters,
	GeometricSpreadingParameters,
	NearSourceSaturationParameters,
	AnelasticAttenuationParameters,
	PathParameters,
	SiteParameters,
	RandomVibrationParameters,
    period,
    transfer,
    transfer!,
    squared_transfer,
    squared_transfer!,
    site_amplification,
	kappa_filter,
    corner_frequency,
	geometric_spreading,
	near_source_saturation,
	equivalent_point_source_distance,
	anelastic_attenuation,
	fourier_constant,
	fourier_source,
	fourier_source_shape,
	fourier_path,
	fourier_attenuation,
	fourier_site,
	combined_kappa_frequency,
	fourier_spectral_ordinate,
	fourier_spectrum,
	fourier_spectrum!,
	spectral_moment,
	spectral_moments,
	excitation_duration,
	rms_duration,
	peak_factor,
	rvt_response_spectral_ordinate,
	rvt_response_spectrum,
	rvt_response_spectrum!


# Write your package code here.
include("oscillator/PJSoscillator.jl")
include("fourier/PJSfourierParameters.jl")
include("rvt/PJSrandomVibrationParameters.jl")
include("fourier/PJSsite.jl")
include("fourier/PJSsource.jl")
include("fourier/PJSpath.jl")
include("fourier/PJSfourierSpectrum.jl")
include("duration/PJSduration.jl")
include("rvt/PJSintegration.jl")
include("rvt/PJSpeakFactor.jl")
include("rvt/PJSrandomVibration.jl")


end
