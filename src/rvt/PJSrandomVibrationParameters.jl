
"""
	RandomVibrationParameters

Struct holding parameters/methods for Random Vibration Theory.

- `pf_method` is the method used for peak factor computation
	- `:DK80` (default) is Der Kiureghian (1980), building on Vanmarcke (1975)
	- `:CL56` is Cartwright Longuet-Higgins (1956)
- `dur_ex` is the model for excitation duration
	- `:BT14` (default) is the Boore & Thompson (2014) model for ACRs - note that this is adpated to work with `r_ps`
	- `:BT15` is the Boore & Thompson (2015) model for SCRs
	- `:BE23` is an excitation duration model from Ben Edwards (2023), suggested for a South African NPP project
- `dur_rms` is the model for rms duration
	- `:BT12` is the Boore & Thompson (2012) model
	- `:BT15` (default) is the Boore & Thompson (2015) model
	- `:LP99` is the Liu & Pezeshk (1999) model linking rms, excitation and oscillator durations
- `dur_region` is the region specified for the duration model
	- `:ACR` (default) is active crustal regions (like western North America)
	- `:SCR` is stable crustal regions (like eastern North America)

Note that only certain combinations are meaningful:
- `:CL56` peak factor method should be paired with `:BT12` or `:LP99` rms duration 
- `:DK80` peak factor method should be paired with `:BT15` rms duration
Constructors that take only the peak factor as input, or the peak factor and duration region automatically assign the appropriate rms duration method.
"""
struct RandomVibrationParameters
    pf_method::Symbol
    dur_ex::Symbol
    dur_rms::Symbol
    dur_region::Symbol
end

RandomVibrationParameters() = RandomVibrationParameters(:DK80, :BT15, :BT15, :ACR)
RandomVibrationParameters(pf) = RandomVibrationParameters(pf, :BT15, ((pf == :DK80) ? :BT15 : :BT12), :ACR)
RandomVibrationParameters(pf, region) = RandomVibrationParameters(pf, :BT15, ((pf == :DK80) ? :BT15 : :BT12), region)


