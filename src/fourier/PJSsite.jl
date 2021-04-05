

# Boore (2016) WUS generic rock amplification function (relative to β=3.5km/s and ρ=2.72g/cm^3 source properties) for Vs30 of 760 m/s
const fii_b16 = [ 0.010, 0.015, 0.021, 0.031, 0.045, 0.065, 0.095, 0.138, 0.200, 0.291, 0.423, 0.615, 0.894, 1.301, 1.892, 2.751, 4.000, 5.817, 8.459, 12.301, 17.889, 26.014, 37.830, 55.012, 80.000 ]
const Aii_b16 = [ 1.00, 1.01, 1.02, 1.02, 1.04, 1.06, 1.09, 1.13, 1.18, 1.25, 1.32, 1.41, 1.51, 1.64, 1.80, 1.99, 2.18, 2.38, 2.56, 2.75, 2.95, 3.17, 3.42, 3.68, 3.96 ]

"""
	boore_2016_generic_amplification(f)

Computes the generic crustal amplification for a WUS velocity profile with Vs30 of 760 m/s using the Boore (2016) velocity model.

# Examples
```julia-repl
	f = 5.0
	Af = boore_2016_generic_amplification(f)
```
"""
function boore_2016_generic_amplification(f::T) where T<:Real
    if isnan(f)
        return T(NaN)
    else
        if f <= 0.01
            return 1.0
        elseif f >= 80.0
            return 3.96
        else
            for i = 1:length(fii_b16)
                @inbounds if fii_b16[i] > f
                    j = i-1
                    @inbounds if fii_b16[j] == f
                        @inbounds amp = Aii_b16[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_b16[j]) + (f-fii_b16[j])*log(Aii_b16[i]/Aii_b16[j])/(fii_b16[i]-fii_b16[j])
                        return exp(lnAmp)
                    end
                end
            end
        end
    end
end


# Al Atik & Abrahamson (2021) generic rock amplification function for Vs30 of 760 m/s obtained from inversion of the Chiou & Youngs (2014) response spectral model
# note that the ordinate at (0.01, 1.0) has been added to the amp function that was provided
const fii_aa21_cy14 = [ 0.01, 0.1, 0.102329, 0.104713, 0.107152, 0.109648, 0.112202, 0.114815, 0.11749, 0.120226, 0.123027, 0.125893, 0.128825, 0.131826, 0.134896, 0.138038, 0.141254, 0.144544, 0.147911, 0.151356, 0.154882, 0.158489, 0.162181, 0.165959, 0.169824, 0.17378, 0.177828, 0.18197, 0.186209, 0.190546, 0.194984, 0.199526, 0.204174, 0.20893, 0.213796, 0.218776, 0.223872, 0.229087, 0.234423, 0.239883, 0.245471, 0.251189, 0.25704, 0.263027, 0.269153, 0.275423, 0.281838, 0.288403, 0.295121, 0.301995, 0.30903, 0.316228, 0.323594, 0.331131, 0.338844, 0.346737, 0.354813, 0.363078, 0.371535, 0.380189, 0.389045, 0.398107, 0.40738, 0.416869, 0.42658, 0.436516, 0.446684, 0.457088, 0.467735, 0.47863, 0.489779, 0.501187, 0.512861, 0.524807, 0.537032, 0.549541, 0.562341, 0.57544, 0.588844, 0.60256, 0.616595, 0.630957, 0.645654, 0.660693, 0.676083, 0.691831, 0.707946, 0.724436, 0.74131, 0.758578, 0.776247, 0.794328, 0.81283, 0.831764, 0.851138, 0.870964, 0.891251, 0.912011, 0.933254, 0.954993, 0.977237, 1.0, 1.023293, 1.047129, 1.071519, 1.096478, 1.122018, 1.148153, 1.174897, 1.202264, 1.230269, 1.258926, 1.28825, 1.318257, 1.348963, 1.380384, 1.412537, 1.44544, 1.479108, 1.513561, 1.548816, 1.584893, 1.62181, 1.659587, 1.698244, 1.737801, 1.778279, 1.819701, 1.862087, 1.905461, 1.949844, 1.995262, 2.041738, 2.089296, 2.137962, 2.187761, 2.238721, 2.290868, 2.344229, 2.398833, 2.454709, 2.511886, 2.570396, 2.630268, 2.691535, 2.754228, 2.818383, 2.884031, 2.951209, 3.019952, 3.090296, 3.162278, 3.235937, 3.311311, 3.388441, 3.467368, 3.548134, 3.63078, 3.715352, 3.801893, 3.890451, 3.981071, 4.073803, 4.168694, 4.265795, 4.365158, 4.466835, 4.570881, 4.677351, 4.7863, 4.897787, 5.011872, 5.128613, 5.248074, 5.370318, 5.495409, 5.623413, 5.754399, 5.888436, 6.025596, 6.165949, 6.309573, 6.456542, 6.606934, 6.760828, 6.918308, 7.079456, 7.24436, 7.413103, 7.585776, 7.762471, 7.943282, 8.128304, 8.317636, 8.511379, 8.709635, 8.912507, 9.120107, 9.332541, 9.549923, 9.772372, 10.0, 10.23293, 10.471284, 10.715192, 10.96478, 11.220183, 11.481534, 11.748973, 12.022642, 12.302684, 12.589251, 12.882492, 13.182563, 13.489624, 13.80384, 14.12537, 14.454392, 14.79108, 15.135614, 15.48817, 15.848933, 16.218101, 16.59587, 16.98244, 17.37801, 17.782793, 18.19701, 18.62087, 19.05461, 19.498443, 19.952621, 20.41738, 20.89296, 21.37962, 21.877611, 22.38721, 22.908672, 23.442283, 23.988321, 24.4608, 25.032, 25.6166, 26.2148, 26.827, 27.4534, 28.0945, 28.7506, 29.422, 30.109, 30.8121, 31.5317, 32.268, 33.0215, 33.7926, 34.5818, 35.3893, 36.2157, 37.0614, 37.9269, 38.8126, 39.7189, 40.6465, 41.5956, 42.567, 43.561, 44.5782, 45.6192, 46.6845, 47.7747, 48.8903, 50.032, 51.2004, 52.396, 53.6196, 54.8717, 56.1531, 57.4643, 58.8062, 60.1795, 61.5848, 63.023, 64.4947, 66.0007, 67.542, 69.1192, 70.7333, 72.3851, 74.0754, 75.8052, 77.5754, 79.387, 81.2409, 83.138, 85.0794, 87.0662, 89.0994, 91.18, 93.3093, 95.4883, 97.7181, 100.0 ]
const Aii_aa21_cy14 = [ 1.0, 1.264743, 1.271871, 1.279014, 1.286166, 1.293327, 1.300495, 1.307668, 1.314844, 1.322022, 1.3292, 1.33638, 1.343563, 1.350753, 1.357952, 1.365166, 1.3724, 1.379659, 1.386946, 1.394264, 1.401614, 1.408999, 1.416419, 1.423873, 1.431361, 1.438883, 1.446437, 1.454019, 1.461627, 1.469255, 1.476897, 1.484547, 1.492199, 1.499845, 1.507478, 1.515089, 1.52267, 1.530214, 1.537712, 1.545156, 1.552541, 1.55986, 1.567111, 1.574289, 1.581393, 1.588421, 1.595371, 1.602243, 1.609038, 1.615756, 1.622398, 1.628966, 1.635461, 1.641884, 1.64824, 1.654529, 1.660755, 1.666921, 1.673029, 1.679083, 1.685087, 1.691043, 1.696955, 1.702826, 1.708661, 1.714462, 1.720234, 1.725979, 1.731701, 1.737405, 1.743093, 1.748769, 1.754436, 1.7601, 1.765763, 1.771428, 1.777101, 1.782784, 1.788483, 1.7942, 1.79994, 1.805707, 1.811506, 1.81734, 1.823215, 1.829134, 1.835103, 1.841123, 1.847199, 1.853331, 1.85952, 1.865766, 1.872069, 1.878427, 1.884838, 1.891302, 1.897816, 1.904376, 1.910981, 1.917627, 1.924312, 1.931032, 1.937785, 1.944567, 1.951374, 1.958204, 1.965053, 1.971919, 1.978799, 1.98569, 1.992588, 1.999493, 2.0064, 2.013308, 2.020215, 2.02712, 2.034019, 2.040913, 2.047799, 2.054677, 2.061546, 2.068404, 2.075252, 2.082088, 2.088913, 2.095727, 2.102529, 2.109321, 2.116101, 2.122872, 2.129632, 2.136385, 2.14313, 2.149868, 2.156601, 2.163329, 2.170056, 2.176782, 2.183509, 2.190238, 2.196973, 2.203714, 2.210465, 2.217228, 2.224004, 2.230798, 2.237612, 2.244448, 2.25131, 2.258201, 2.265124, 2.272083, 2.279082, 2.286124, 2.293212, 2.300353, 2.307549, 2.314805, 2.322127, 2.329519, 2.336986, 2.344534, 2.352169, 2.359894, 2.367718, 2.375645, 2.383683, 2.391837, 2.400114, 2.408522, 2.417069, 2.425761, 2.434606, 2.443614, 2.452794, 2.462154, 2.471704, 2.481454, 2.491416, 2.501599, 2.512017, 2.522681, 2.533593, 2.544623, 2.555702, 2.566822, 2.577981, 2.589182, 2.600422, 2.611704, 2.623026, 2.634389, 2.645793, 2.657238, 2.668725, 2.680254, 2.691824, 2.703436, 2.71509, 2.726786, 2.738525, 2.750306, 2.76213, 2.773996, 2.785906, 2.797859, 2.809855, 2.821894, 2.833977, 2.846104, 2.858274, 2.870489, 2.882748, 2.895052, 2.9074, 2.919793, 2.93223, 2.944713, 2.957241, 2.969814, 2.982433, 2.995097, 3.007808, 3.020564, 3.033367, 3.046216, 3.059111, 3.072054, 3.085043, 3.098079, 3.111162, 3.124293, 3.137471, 3.150697, 3.163971, 3.177292, 3.190662, 3.204081, 3.217547, 3.231063, 3.242526, 3.256161, 3.269851, 3.28359, 3.29738, 3.311217, 3.325106, 3.339046, 3.353036, 3.367075, 3.381166, 3.395309, 3.409502, 3.423746, 3.438042, 3.452391, 3.466789, 3.48124, 3.495744, 3.5103, 3.524909, 3.539569, 3.554285, 3.56905, 3.583871, 3.598744, 3.613671, 3.628651, 3.643686, 3.658775, 3.673917, 3.689114, 3.704366, 3.719671, 3.735033, 3.750447, 3.765917, 3.781441, 3.797022, 3.812659, 3.82835, 3.844099, 3.859904, 3.875764, 3.891682, 3.907655, 3.923684, 3.939771, 3.955913, 3.972112, 3.988367, 4.004681, 4.021051, 4.037478, 4.053962, 4.070504, 4.087104, 4.103761, 4.120476, 4.13725, 4.154081, 4.170967 ]



const fii_aa21_cy14_a = [ 0.01, 0.1, 0.102329, 0.104713, 0.107152, 0.109648, 0.112202, 0.114815, 0.11749, 0.120226, 0.123027, 0.125893, 0.128825, 0.131826, 0.134896, 0.138038, 0.141254, 0.144544, 0.147911, 0.151356, 0.154882, 0.158489, 0.162181, 0.165959, 0.169824, 0.17378, 0.177828, 0.18197, 0.186209, 0.190546, 0.194984, 0.199526 ]
const Aii_aa21_cy14_a = [ 1.0, 1.264743, 1.271871, 1.279014, 1.286166, 1.293327, 1.300495, 1.307668, 1.314844, 1.322022, 1.3292, 1.33638, 1.343563, 1.350753, 1.357952, 1.365166, 1.3724, 1.379659, 1.386946, 1.394264, 1.401614, 1.408999, 1.416419, 1.423873, 1.431361, 1.438883, 1.446437, 1.454019, 1.461627, 1.469255, 1.476897, 1.484547 ]

const fii_aa21_cy14_b = [ 0.199526, 0.204174, 0.20893, 0.213796, 0.218776, 0.223872, 0.229087, 0.234423, 0.239883, 0.245471, 0.251189, 0.25704, 0.263027, 0.269153, 0.275423, 0.281838, 0.288403, 0.295121, 0.301995, 0.30903, 0.316228, 0.323594, 0.331131, 0.338844, 0.346737, 0.354813, 0.363078, 0.371535, 0.380189, 0.389045, 0.398107 ]
const Aii_aa21_cy14_b = [ 1.484547, 1.492199, 1.499845, 1.507478, 1.515089, 1.52267, 1.530214, 1.537712, 1.545156, 1.552541, 1.55986, 1.567111, 1.574289, 1.581393, 1.588421, 1.595371, 1.602243, 1.609038, 1.615756, 1.622398, 1.628966, 1.635461, 1.641884, 1.64824, 1.654529, 1.660755, 1.666921, 1.673029, 1.679083, 1.685087, 1.691043 ]

const fii_aa21_cy14_c = [ 0.398107, 0.40738, 0.416869, 0.42658, 0.436516, 0.446684, 0.457088, 0.467735, 0.47863, 0.489779, 0.501187, 0.512861, 0.524807, 0.537032, 0.549541, 0.562341, 0.57544, 0.588844, 0.60256, 0.616595, 0.630957, 0.645654, 0.660693, 0.676083, 0.691831, 0.707946, 0.724436, 0.74131, 0.758578, 0.776247, 0.794328 ]
const Aii_aa21_cy14_c = [ 1.691043, 1.696955, 1.702826, 1.708661, 1.714462, 1.720234, 1.725979, 1.731701, 1.737405, 1.743093, 1.748769, 1.754436, 1.7601, 1.765763, 1.771428, 1.777101, 1.782784, 1.788483, 1.7942, 1.79994, 1.805707, 1.811506, 1.81734, 1.823215, 1.829134, 1.835103, 1.841123, 1.847199, 1.853331, 1.85952, 1.865766 ]

const fii_aa21_cy14_d = [ 0.794328, 0.81283, 0.831764, 0.851138, 0.870964, 0.891251, 0.912011, 0.933254, 0.954993, 0.977237, 1, 1.023293, 1.047129, 1.071519, 1.096478, 1.122018, 1.148153, 1.174897, 1.202264, 1.230269, 1.258926, 1.28825, 1.318257, 1.348963, 1.380384, 1.412537, 1.44544, 1.479108, 1.513561, 1.548816, 1.584893 ]
const Aii_aa21_cy14_d = [ 1.865766, 1.872069, 1.878427, 1.884838, 1.891302, 1.897816, 1.904376, 1.910981, 1.917627, 1.924312, 1.931032, 1.937785, 1.944567, 1.951374, 1.958204, 1.965053, 1.971919, 1.978799, 1.98569, 1.992588, 1.999493, 2.0064, 2.013308, 2.020215, 2.02712, 2.034019, 2.040913, 2.047799, 2.054677, 2.061546, 2.068404 ]

const fii_aa21_cy14_e = [ 1.584893, 1.62181, 1.659587, 1.698244, 1.737801, 1.778279, 1.819701, 1.862087, 1.905461, 1.949844, 1.995262, 2.041738, 2.089296, 2.137962, 2.187761, 2.238721, 2.290868, 2.344229, 2.398833, 2.454709, 2.511886, 2.570396, 2.630268, 2.691535, 2.754228, 2.818383, 2.884031, 2.951209, 3.019952, 3.090296, 3.162278 ]
const Aii_aa21_cy14_e = [ 2.068404, 2.075252, 2.082088, 2.088913, 2.095727, 2.102529, 2.109321, 2.116101, 2.122872, 2.129632, 2.136385, 2.14313, 2.149868, 2.156601, 2.163329, 2.170056, 2.176782, 2.183509, 2.190238, 2.196973, 2.203714, 2.210465, 2.217228, 2.224004, 2.230798, 2.237612, 2.244448, 2.25131, 2.258201, 2.265124, 2.272083 ]

const fii_aa21_cy14_f = [ 3.162278, 3.235937, 3.311311, 3.388441, 3.467368, 3.548134, 3.63078, 3.715352, 3.801893, 3.890451, 3.981071, 4.073803, 4.168694, 4.265795, 4.365158, 4.466835, 4.570881, 4.677351, 4.7863, 4.897787, 5.011872, 5.128613, 5.248074, 5.370318, 5.495409, 5.623413, 5.754399, 5.888436, 6.025596, 6.165949, 6.309573 ]
const Aii_aa21_cy14_f = [ 2.272083, 2.279082, 2.286124, 2.293212, 2.300353, 2.307549, 2.314805, 2.322127, 2.329519, 2.336986, 2.344534, 2.352169, 2.359894, 2.367718, 2.375645, 2.383683, 2.391837, 2.400114, 2.408522, 2.417069, 2.425761, 2.434606, 2.443614, 2.452794, 2.462154, 2.471704, 2.481454, 2.491416, 2.501599, 2.512017, 2.522681 ]

const fii_aa21_cy14_g = [ 6.309573, 6.456542, 6.606934, 6.760828, 6.918308, 7.079456, 7.24436, 7.413103, 7.585776, 7.762471, 7.943282, 8.128304, 8.317636, 8.511379, 8.709635, 8.912507, 9.120107, 9.332541, 9.549923, 9.772372, 10.0, 10.23293, 10.471284, 10.715192, 10.96478, 11.220183, 11.481534, 11.748973, 12.022642, 12.302684, 12.589251 ]
const Aii_aa21_cy14_g = [ 2.522681, 2.533593, 2.544623, 2.555702, 2.566822, 2.577981, 2.589182, 2.600422, 2.611704, 2.623026, 2.634389, 2.645793, 2.657238, 2.668725, 2.680254, 2.691824, 2.703436, 2.71509, 2.726786, 2.738525, 2.750306, 2.76213, 2.773996, 2.785906, 2.797859, 2.809855, 2.821894, 2.833977, 2.846104, 2.858274, 2.870489 ]

const fii_aa21_cy14_h = [ 12.589251, 12.882492, 13.182563, 13.489624, 13.80384, 14.12537, 14.454392, 14.79108, 15.135614, 15.48817, 15.848933, 16.218101, 16.59587, 16.98244, 17.37801, 17.782793, 18.19701, 18.62087, 19.05461, 19.498443, 19.952621, 20.41738, 20.89296, 21.37962, 21.877611, 22.38721, 22.908672, 23.442283, 23.988321, 24.4608, 25.032 ]
const Aii_aa21_cy14_h = [ 2.870489, 2.882748, 2.895052, 2.9074, 2.919793, 2.93223, 2.944713, 2.957241, 2.969814, 2.982433, 2.995097, 3.007808, 3.020564, 3.033367, 3.046216, 3.059111, 3.072054, 3.085043, 3.098079, 3.111162, 3.124293, 3.137471, 3.150697, 3.163971, 3.177292, 3.190662, 3.204081, 3.217547, 3.231063, 3.242526, 3.256161 ]

const fii_aa21_cy14_i = [ 25.032, 25.6166, 26.2148, 26.827, 27.4534, 28.0945, 28.7506, 29.422, 30.109, 30.8121, 31.5317, 32.268, 33.0215, 33.7926, 34.5818, 35.3893, 36.2157, 37.0614, 37.9269, 38.8126, 39.7189, 40.6465, 41.5956, 42.567, 43.561, 44.5782, 45.6192, 46.6845, 47.7747, 48.8903, 50.032 ]
const Aii_aa21_cy14_i = [ 3.256161, 3.269851, 3.28359, 3.29738, 3.311217, 3.325106, 3.339046, 3.353036, 3.367075, 3.381166, 3.395309, 3.409502, 3.423746, 3.438042, 3.452391, 3.466789, 3.48124, 3.495744, 3.5103, 3.524909, 3.539569, 3.554285, 3.56905, 3.583871, 3.598744, 3.613671, 3.628651, 3.643686, 3.658775, 3.673917, 3.689114 ]

const fii_aa21_cy14_j = [ 50.032, 51.2004, 52.396, 53.6196, 54.8717, 56.1531, 57.4643, 58.8062, 60.1795, 61.5848, 63.023, 64.4947, 66.0007, 67.542, 69.1192, 70.7333, 72.3851, 74.0754, 75.8052, 77.5754, 79.387, 81.2409, 83.138, 85.0794, 87.0662, 89.0994, 91.18, 93.3093, 95.4883, 97.7181, 100.0 ]
const Aii_aa21_cy14_j = [ 3.689114, 3.704366, 3.719671, 3.735033, 3.750447, 3.765917, 3.781441, 3.797022, 3.812659, 3.82835, 3.844099, 3.859904, 3.875764, 3.891682, 3.907655, 3.923684, 3.939771, 3.955913, 3.972112, 3.988367, 4.004681, 4.021051, 4.037478, 4.053962, 4.070504, 4.087104, 4.103761, 4.120476, 4.13725, 4.154081, 4.170967 ]


"""
	alatik_2021_cy14_inverted_amplification_seg(f::T) where T<:Real

Computes the generic crustal amplification for a velocity profile with Vs30 of 760 m/s obtained from inverting the Chiou & Youngs (2014) ground-motion model.

Amplifications provided ahead of Al Atik & Abrahamson (2021) publication.

Same as alatik_2021_cy14_inverted_amplification, just segmented to speed up the indexing.

# Examples
```julia-repl
	f = 5.0
	Af = alatik_2021_cy14_inverted_amplification_seg(f)
```
"""
function alatik_2021_cy14_inverted_amplification_seg(f::T) where T<:Real
    if isnan(f)
        return T(NaN)
    else
        if f <= 0.01
            return 1.0
        elseif f <= 0.199526
            for i in 1:length(fii_aa21_cy14_a)
                @inbounds if fii_aa21_cy14_a[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_a[j] == f
                        @inbounds amp = Aii_aa21_cy14_a[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_a[j]) + (f-fii_aa21_cy14_a[j])*log(Aii_aa21_cy14_a[i]/Aii_aa21_cy14_a[j])/(fii_aa21_cy14_a[i]-fii_aa21_cy14_a[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_a[end]
        elseif f <= 0.398107
            for i in 1:length(fii_aa21_cy14_b)
                @inbounds if fii_aa21_cy14_b[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_b[j] == f
                        @inbounds amp = Aii_aa21_cy14_b[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_b[j]) + (f-fii_aa21_cy14_b[j])*log(Aii_aa21_cy14_b[i]/Aii_aa21_cy14_b[j])/(fii_aa21_cy14_b[i]-fii_aa21_cy14_b[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_b[end]
        elseif f <= 0.794328
            for i in 1:length(fii_aa21_cy14_c)
                @inbounds if fii_aa21_cy14_c[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_c[j] == f
                        @inbounds amp = Aii_aa21_cy14_c[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_c[j]) + (f-fii_aa21_cy14_c[j])*log(Aii_aa21_cy14_c[i]/Aii_aa21_cy14_c[j])/(fii_aa21_cy14_c[i]-fii_aa21_cy14_c[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_c[end]
        elseif f <= 1.584893
            for i in 1:length(fii_aa21_cy14_d)
                @inbounds if fii_aa21_cy14_d[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_d[j] == f
                        @inbounds amp = Aii_aa21_cy14_d[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_d[j]) + (f-fii_aa21_cy14_d[j])*log(Aii_aa21_cy14_d[i]/Aii_aa21_cy14_d[j])/(fii_aa21_cy14_d[i]-fii_aa21_cy14_d[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_d[end]
        elseif f <= 3.162278
            for i in 1:length(fii_aa21_cy14_e)
                @inbounds if fii_aa21_cy14_e[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_e[j] == f
                        @inbounds amp = Aii_aa21_cy14_e[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_e[j]) + (f-fii_aa21_cy14_e[j])*log(Aii_aa21_cy14_e[i]/Aii_aa21_cy14_e[j])/(fii_aa21_cy14_e[i]-fii_aa21_cy14_e[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_e[end]
        elseif f <= 6.309573
            for i in 1:length(fii_aa21_cy14_f)
                @inbounds if fii_aa21_cy14_f[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_f[j] == f
                        @inbounds amp = Aii_aa21_cy14_f[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_f[j]) + (f-fii_aa21_cy14_f[j])*log(Aii_aa21_cy14_f[i]/Aii_aa21_cy14_f[j])/(fii_aa21_cy14_f[i]-fii_aa21_cy14_f[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_f[end]
        elseif f <= 12.589251
            for i in 1:length(fii_aa21_cy14_g)
                @inbounds if fii_aa21_cy14_g[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_g[j] == f
                        @inbounds amp = Aii_aa21_cy14_g[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_g[j]) + (f-fii_aa21_cy14_g[j])*log(Aii_aa21_cy14_g[i]/Aii_aa21_cy14_g[j])/(fii_aa21_cy14_g[i]-fii_aa21_cy14_g[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_g[end]
        elseif f <= 25.032
            for i in 1:length(fii_aa21_cy14_h)
                @inbounds if fii_aa21_cy14_h[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_h[j] == f
                        @inbounds amp = Aii_aa21_cy14_h[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_h[j]) + (f-fii_aa21_cy14_h[j])*log(Aii_aa21_cy14_h[i]/Aii_aa21_cy14_h[j])/(fii_aa21_cy14_h[i]-fii_aa21_cy14_h[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_h[end]
        elseif f <= 50.032
            for i in 1:length(fii_aa21_cy14_i)
                @inbounds if fii_aa21_cy14_i[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_i[j] == f
                        @inbounds amp = Aii_aa21_cy14_i[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_i[j]) + (f-fii_aa21_cy14_i[j])*log(Aii_aa21_cy14_i[i]/Aii_aa21_cy14_i[j])/(fii_aa21_cy14_i[i]-fii_aa21_cy14_i[j])
                        return exp(lnAmp)
                    end
                end
            end
            return Aii_aa21_cy14_i[end]
        elseif f < 100.0
            for i in 1:length(fii_aa21_cy14_j)
                @inbounds if fii_aa21_cy14_j[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14_j[j] == f
                        @inbounds amp = Aii_aa21_cy14_j[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14_j[j]) + (f-fii_aa21_cy14_j[j])*log(Aii_aa21_cy14_j[i]/Aii_aa21_cy14_j[j])/(fii_aa21_cy14_j[i]-fii_aa21_cy14_j[j])
                        return exp(lnAmp)
                    end
                end
            end
        else
            return 4.170967
        end
    end
end

"""
	alatik_2021_cy14_inverted_amplification(f::T) where T<:Real

Computes the generic crustal amplification for a velocity profile with Vs30 of 760 m/s obtained from inverting the Chiou & Youngs (2014) ground-motion model.

Amplifications provided ahead of Al Atik & Abrahamson (2021) publication.

# Examples
```julia-repl
	f = 5.0
	Af = alatik_2021_cy14_inverted_amplification(f)
```
"""
function alatik_2021_cy14_inverted_amplification(f::T) where T<:Real
    if isnan(f)
        return T(NaN)
    else
        if f <= 0.01
            return 1.0
        elseif f >= 100.0
            return 4.170967
        else
            for i in 1:length(fii_aa21_cy14)
                @inbounds if fii_aa21_cy14[i] > f
                    j = i-1
                    @inbounds if fii_aa21_cy14[j] == f
                        @inbounds amp = Aii_aa21_cy14[j]
                        return amp
                    else
                        @inbounds lnAmp = log(Aii_aa21_cy14[j]) + (f-fii_aa21_cy14[j])*log(Aii_aa21_cy14[i]/Aii_aa21_cy14[j])/(fii_aa21_cy14[i]-fii_aa21_cy14[j])
                        return exp(lnAmp)
                    end
                end
            end
        end
    end
end



"""
	unit_generic_amplification()

Function used to represent null impedance effects. Use this function to represent site effects when no amplification is desired. Returns unity.

# Examples
```julia-repl
	Af = unit_generic_amplification()
	# Af will equal 1.0
```
"""
function unit_generic_amplification()
    return 1.0
end


"""
	site_amplification(f::Real, amp_model::Symbol)

Computes the site amplification (impedance) for a given frequency `f`. Requires the keyword argument `amp_model` as a `String` and defaults to the Boore (2016) model. Currently, any other string passed to the function will return the unit amplification

# Examples
```julia-repl
	f = 5.0
	# returns the amplification from AlAtik (2021) in both cases
	Af = site_amplification(f)
    Af = site_amplification(f; amp_model=:AlAtik2021_cy14)
    # returns the Boore (2016) amplification
	Af = site_amplification(f; amp_model=:Boore2016)
	# returns 1.0
	Af = site_amplification(f; amp_model=:Unit)
```
"""
function site_amplification(f::T, model::Symbol) where T<:Real
    if model == :Unit
        return unit_generic_amplification() * oneunit(T)
    elseif model == :Boore2016
        return boore_2016_generic_amplification(f)
    elseif model == :AlAtik2021_cy14
        return alatik_2021_cy14_inverted_amplification_seg(f)
	else
		return T(NaN)
    end
end

"""
	site_amplification(f, site::SiteParameters)

Computes the site amplification (impedance) for a given frequency `f`.
"""
site_amplification(f, site::SiteParameters) = site_amplification(f, site.model)

"""
	site_amplification(f, fas::FourierParameters)

Computes the site amplification (impedance) for a given frequency `f`.
"""
site_amplification(f, fas::FourierParameters) = site_amplification(f, fas.site)


"""
    kappa_filter(f, site::SiteParameters)

Kappa filter for a given frequency `f`
"""
function kappa_filter(f, site::SiteParameters)
    return exp(-π*f*site.κ0)
end