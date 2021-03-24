
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
	get_parametric_type,
    period,
    transfer,
    transfer!,
    squared_transfer,
    squared_transfer!,
    site_amplification,
	kappa_filter,
	magnitude_to_moment,
    corner_frequency,
	corner_frequency_brune,
	corner_frequency_atkinson_silva_2000,
	geometric_spreading,
	geometric_spreading_piecewise,
	geometric_spreading_cy14,
	geometric_spreading_cy14mod,
	near_source_saturation,
	equivalent_point_source_distance,
	rupture_distance_from_equivalent_point_source_distance,
	anelastic_attenuation,
	fourier_constant,
	fourier_source,
	fourier_source_shape,
	fourier_path,
	fourier_attenuation,
	fourier_site,
	fourier_spectral_ordinate,
	squared_fourier_spectral_ordinate,
	fourier_spectrum,
	squared_fourier_spectrum,
	fourier_spectrum!,
	squared_fourier_spectrum!,
	combined_kappa_frequency,
	boore_thompson_2014,
	simpsons_rule,
	trapezoidal_rule,
	spectral_moment,
	spectral_moments,
	spectral_moments_ln,
	spectral_moments_gk,
	zeros_extrema_frequencies,
	zeros_extrema_numbers,
	boore_thompson_2012_coefs,
	boore_thompson_2012_base,
	boore_thompson_2012,
	boore_thompson_2015_coefs,
	boore_thompson_2015_base,
	boore_thompson_2015,
	excitation_duration,
	rms_duration,
	peak_factor,
	peak_factor_cl56,
	peak_factor_cl56_gk,
	peak_factor_dk80,
	peak_factor_dk80_gk,
	find_integrand_value,
	peak_factor_integrand,
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
