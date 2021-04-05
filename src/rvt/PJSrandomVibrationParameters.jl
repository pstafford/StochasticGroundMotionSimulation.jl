
"""
	RandomVibrationParameters

Struct holding parameters/methods for Random Vibration Theory.

- `pf_method` is the method used for peak factor computation
	- `:DK80` (default) is Der Kiureghian (1980), building on Vanmarcke (1975)
	- `:CL56` is Cartwright Longuet-Higgins (1956)
- `dur_ex` is the model for excitation duration
	- `:BT14` (default) is the Boore & Thompson (2014) model - note that this is adpated to work with `r_ps`
- `dur_rms` is the model for rms duration
	- `:BT12` is the Boore & Thompson (2012) model
	- `:BT15` (default) is the Boore & Thompson (2015) model
- `dur_region` is the region specified for the duration model
	- `:WNA` (default) is western North America
	- `:ENA` is eastern North America
"""
struct RandomVibrationParameters
	pf_method::Symbol
	dur_ex::Symbol
	dur_rms::Symbol
	dur_region::Symbol
end

RandomVibrationParameters() = RandomVibrationParameters(:DK80, :BT14, :BT15, :WNA)
RandomVibrationParameters(pf) = RandomVibrationParameters(pf, :BT14, ((pf == :DK80) ? :BT15 : :BT12), :WNA)
