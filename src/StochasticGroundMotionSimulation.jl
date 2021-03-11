module StochasticGroundMotionSimulation

using Interpolations
using Roots
using ForwardDiff
using ForwardDiff: Dual


export Oscillator,
        FASParams,
        FASParamsGeo,
        FASParamsQr,
        FASParamsGeoQr,
        period,
        transfer,
        transfer!,
        squared_transfer,
        squared_transfer!,
        site_amplification,
        finite_fault_factor,
    	geometric_spreading,
    	geometric_spreading_cy,
    	fourier_constant,
    	fourier_source,
    	fourier_source_shape,
    	fourier_path,
	fourier_path_cy,
    	fourier_attenuation,
    	fourier_site,
    	fourier_spectral_ordinate,
    	fourier_spectrum,
    	fourier_spectrum!,
    	squared_fourier_spectrum!,
	squared_fourier_spectrum_cy!,
    	boore_thompson_2014,
    	combined_kappa_frequency,
    	spectral_moment,
    	spectral_moments,
	spectral_moments_cy,
    	zeros_extrema_frequencies,
    	zeros_extrema_numbers,
    	boore_thompson_2012_coefs,
    	boore_thompson_2012_base,
    	boore_thompson_2012,
    	peak_factor,
    	find_integrand_value,
    	peak_factor_integrand,
    	rvt_response_spectral_ordinate,
	rvt_response_spectral_ordinate_cy,
    	rvt_response_spectrum,
	rvt_response_spectrum_cy,
    	rvt_response_spectrum!,
	rvt_response_spectrum_cy!


# Write your package code here.
include("oscillator/PJSoscillator.jl")
include("fourier/PJSfourierParameters.jl")
include("fourier/PJSsite.jl")
include("fourier/PJSsource.jl")
include("fourier/PJSpath.jl")
include("duration/PJSduration.jl")
include("rvt/PJSrandomVibration.jl")


end
