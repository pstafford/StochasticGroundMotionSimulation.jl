
"""
  boore_thompson_2014_path_duration(r_ps::T) where {T<:Real}

Boore & Thompson (2014) path duration model (for ACRs).
Note that this model is the same as the Boore & Thompson (2015) model for ACRs.
"""
function boore_thompson_2014_path_duration(r_ps::T) where {T<:Real}
  # path duration
  if r_ps == 0.0
    return zero(T)
  elseif r_ps > 0.0 && r_ps <= 7.0
    return r_ps / 7.0 * 2.4
  elseif r_ps > 7.0 && r_ps <= 45.0
    return 2.4 + (r_ps - 7.0) / (45.0 - 7.0) * (8.4 - 2.4)
  elseif r_ps > 45.0 && r_ps <= 125.0
    return 8.4 + (r_ps - 45.0) / (125.0 - 45.0) * (10.9 - 8.4)
  elseif r_ps > 125.0 && r_ps <= 175.0
    return 10.9 + (r_ps - 125.0) / (175.0 - 125.0) * (17.4 - 10.9)
  elseif r_ps > 175.0 && r_ps <= 270.0
    return 17.4 + (r_ps - 175.0) / (270.0 - 175.0) * (34.2 - 17.4)
  elseif r_ps > 270.0
    return 34.2 + 0.156 * (r_ps - 270.0)
  else
    return T(NaN)
  end
end


"""
  boore_thompson_2014(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64, T<:Real, U<:Real}

Boore & Thompson (2014) excitation duration model (for ACRs).
Note that this model is the same as the Boore & Thompson (2015) model for ACRs.
"""
function boore_thompson_2014(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64,T<:Real,U<:Real}
  # source duration
  fa, fb, ε = corner_frequency(m, src)
  if src.model == :Atkinson_Silva_2000
    Ds = 0.5 / fa + 0.5 / fb
  else
    Ds = 1.0 / fa
  end
  # path duration
  Dp = boore_thompson_2014_path_duration(r_ps)
  # checks on return type for type stability
  if isnan(Dp)
    if T <: Dual
      return T(NaN)
    elseif U <: Dual
      return U(NaN)
    else
      return S(NaN)
    end
  else
    return Ds + Dp
  end
end

boore_thompson_2014(m, r_ps, fas::FourierParameters) = boore_thompson_2014(m, r_ps, fas.source)



"""
  boore_thompson_2015_path_duration_acr(r_ps::T) where {T<:Real}

Boore & Thompson (2015) path duration model (for ACRs).
Note that this model is the same as the Boore & Thompson (2014) model for ACRs.
"""
boore_thompson_2015_path_duration_acr(r_ps::T) where {T<:Real} = boore_thompson_2014_path_duration(r_ps)


"""
  boore_thompson_2015_path_duration_scr(r_ps::T) where {T<:Real}

Boore & Thompson (2015) path duration model (for SCRs).
"""
function boore_thompson_2015_path_duration_scr(r_ps::T) where {T<:Real}
  # path duration
  if r_ps == 0.0
    return zero(T)
  elseif r_ps > 0.0 && r_ps <= 15.0
    return r_ps / 15.0 * 2.6
  elseif r_ps > 15.0 && r_ps <= 35.0
    return 2.6 + (r_ps - 15.0) / (35.0 - 15.0) * (17.5 - 2.6)
  elseif r_ps > 35.0 && r_ps <= 50.0
    return 17.5 + (r_ps - 35.0) / (50.0 - 35.0) * (25.1 - 17.5)
  elseif r_ps > 50.0 && r_ps <= 125.0
    return 25.1
  elseif r_ps > 125.0 && r_ps <= 200.0
    return 25.1 + (r_ps - 125.0) / (200.0 - 125.0) * (28.5 - 25.1)
  elseif r_ps > 200.0 && r_ps <= 392.0
    return 28.5 + (r_ps - 200.0) / (392.0 - 200.0) * (46.0 - 28.5)
  elseif r_ps > 392.0 && r_ps <= 600.0
    return 46.0 + (r_ps - 392.0) / (600.0 - 392.0) * (69.1 - 46.0)
  elseif r_ps > 600.0
    return 69.1 + 0.111 * (r_ps - 600.0)
  else
    return T(NaN)
  end
end




"""
  boore_thompson_2015(m, r_ps::U, src::SourceParameters{S,T}, region::Symbol) where {S<:Float64, T<:Real, U<:Real}

Boore & Thompson (2015) excitation duration model (for ACRs or SCRs).
"""
function boore_thompson_2015(m, r_ps::U, src::SourceParameters{S,T}, region::Symbol) where {S<:Float64,T<:Real,U<:Real}
  # source duration
  fa, fb, ε = corner_frequency(m, src)
  if src.model == :Atkinson_Silva_2000
    Ds = 0.5 / fa + 0.5 / fb
  else
    Ds = 1.0 / fa
  end
  # path duration
  if region == :ACR
    Dp = boore_thompson_2015_path_duration_acr(r_ps)
  elseif region == :SCR
    Dp = boore_thompson_2015_path_duration_scr(r_ps)
  else 
    Dp = NaN
  end
  # checks on return type for type stability
  if isnan(Dp)
    if T <: Dual
      return T(NaN)
    elseif U <: Dual
      return U(NaN)
    else
      return S(NaN)
    end
  else
    return Ds + Dp
  end
end

boore_thompson_2015(m, r_ps, fas::FourierParameters, region) = boore_thompson_2015(m, r_ps, fas.source, region)
boore_thompson_2015(m, r_ps::U, src::SourceParameters{S,T}, rvt::RandomVibrationParameters) where {S<:Float64, T<:Real, U<:Real} = boore_thompson_2015(m, r_ps, src, rvt.dur_region)
boore_thompson_2015(m, r_ps::U, fas::FourierParameters, rvt::RandomVibrationParameters) where {U<:Real} = boore_thompson_2015(m, r_ps, fas.source, rvt.dur_region)


"""
  edwards_2023_path_duration(r_ps::T) where {T<:Real}

Edwards (2023) path duration model using within South African NPP study.
Note that the model was specified in terms of R_hyp -- so interpret this to be the effective point-source distance.
Further note that this model is apparently based upon work of Afshari & Stewart (2016) who found
a dependence of 5-75% significant duration on distance of 0.07s/km. Edwards then found for the 
Groningen project in the Netherlands (for induced seismicity) that the excitation duration was 
equivalent to the 5-75% duration divided by 0.55. 
However, in Edwards' definition, this 0.55 scale factor is supposed to apply to the 5-75% duration
in its entirety, but he only scales the path duration (not the source duration).
"""
function edwards_2023_path_duration(r_ps::T) where {T<:Real}
  # path duration
  if r_ps <= 0.0
    return zero(T)
  else
    return 0.13 * r_ps
  end
end


"""
  edwards_2023(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64, T<:Real, U<:Real}

Edwards (2023) excitation duration model proposed for use within South African NPP study.
"""
function edwards_2023(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64,T<:Real,U<:Real}
  # source duration
  fa, fb, ε = corner_frequency(m, src)
  # this is not explicitly considered in the Edwards model, but keep the same functionality as other models for consistency
  if src.model == :Atkinson_Silva_2000
    Ds = 0.5 / fa + 0.5 / fb
  else
    Ds = 1.0 / fa
  end
  # path duration
  Dp = edwards_2023_path_duration(r_ps)
  return Ds + Dp
end


"""
  uk_path_duration_free(r_ps::T) where {T<:Real}

Path duration model for the UK allowing the plateau segment to have a free slope
"""
function uk_path_duration_free(r_ps::T) where {T<:Real}
  r_ref = [ 0.0, 21.8, 140.4, 254.2 ]
  dg_ref = [ 0.667, 0.002, 0.128, 0.192 ]
  d_path = zero(T)

  if r_ps <= r_ref[1]
    return d_path
  elseif r_ps <= r_ref[2]
    return d_path + dg_ref[1] * r_ps
  elseif r_ps <= r_ref[3]
    return d_path + dg_ref[1] * r_ref[2] + dg_ref[2] * (r_ps - r_ref[2])
  elseif r_ps <= r_ref[4]
    return d_path + dg_ref[1] * r_ref[2] + dg_ref[2] * (r_ref[3] - r_ref[2]) + dg_ref[3] * (r_ps - r_ref[3])
  else
    return d_path + dg_ref[1] * r_ref[2] + dg_ref[2] * (r_ref[3] - r_ref[2]) + dg_ref[3] * (r_ref[4] - r_ref[3]) + dg_ref[4] * (r_ps - r_ref[4])
  end
end


"""
  uk_path_duration_fixed(r_ps::T) where {T<:Real}

Path duration model for the UK forcing the plateau segment to have a flat slope
"""
function uk_path_duration_fixed(r_ps::T) where {T<:Real}
  r_ref = [0.0, 28.3, 140.5, 254.1]
  dg_ref = [0.519, 0.0, 0.13, 0.193]
  d_path = zero(T)

  if r_ps <= r_ref[1]
    return d_path
  elseif r_ps <= r_ref[2]
    return d_path + dg_ref[1] * r_ps
  elseif r_ps <= r_ref[3]
    return d_path + dg_ref[1] * r_ref[2] + dg_ref[2] * (r_ps - r_ref[2])
  elseif r_ps <= r_ref[4]
    return d_path + dg_ref[1] * r_ref[2] + dg_ref[2] * (r_ref[3] - r_ref[2]) + dg_ref[3] * (r_ps - r_ref[3])
  else
    return d_path + dg_ref[1] * r_ref[2] + dg_ref[2] * (r_ref[3] - r_ref[2]) + dg_ref[3] * (r_ref[4] - r_ref[3]) + dg_ref[4] * (r_ps - r_ref[4])
  end
end

"""
  uk_duration_free(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64,T<:Real,U<:Real}

Excitation duration for the UK. Combines standard reciprocal corner frequency source duration with a UK path duration. 
"""
function uk_duration_free(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64,T<:Real,U<:Real}
  # source duration
  fa, fb, ε = corner_frequency(m, src)
  # this is not explicitly considered in the Edwards model, but keep the same functionality as other models for consistency
  if src.model == :Atkinson_Silva_2000
    Ds = 0.5 / fa + 0.5 / fb
  else
    Ds = 1.0 / fa
  end
  # path duration
  Dp = uk_path_duration_free(r_ps)
  # combine source and path durations
  return Ds + Dp
end

"""
  uk_duration_fixed(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64,T<:Real,U<:Real}

Excitation duration for the UK. Combines standard reciprocal corner frequency source duration with a UK path duration. 
"""
function uk_duration_fixed(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64,T<:Real,U<:Real}
  # source duration
  fa, fb, ε = corner_frequency(m, src)
  # this is not explicitly considered in the Edwards model, but keep the same functionality as other models for consistency
  if src.model == :Atkinson_Silva_2000
    Ds = 0.5 / fa + 0.5 / fb
  else
    Ds = 1.0 / fa
  end
  # path duration
  Dp = uk_path_duration_fixed(r_ps)
  # combine source and path durations
  return Ds + Dp
end


"""
  excitation_duration(m, r_ps::U, src::SourceParameters{S,T}, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real,U<:Real}

Generic function implementing excitation duration models.

Currently, only the general Boore & Thompson (2014, 2015) models are implemented along with some project-specific models from Edwards (2023) (for South Africa), 
and two UK duration models. 
The first two are both represented within the Boore & Thompson (2015) paper, so just switch path duration based upon `rvt.dur_region`
To activate the Edwards model, use `rvt.dur_ex == :BE23`
For the UK models use `rvt.dur_ex == :UKfree` or `rvt.dur_ex == :UKfixed`
"""
function excitation_duration(m, r_ps::U, src::SourceParameters{S,T}, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real,U<:Real}
  if rvt.dur_ex == :BE23
    # use the BT15 saturation model to obtain an r_rup from the specified r_ps
    # Edwards actually suggested using R_{EFF}, but this metric is heavily dependent upon a particular source/site geometry
    # sat = NearSourceSaturationParameters(:BT15)
    # compute r_rup from r_ps
    # r_rup = rupture_distance_from_equivalent_point_source_distance(r_ps, m, sat)
    return edwards_2023(m, r_ps, src)
  elseif rvt.dur_ex == :UKfree 
    return uk_duration_free(m, r_ps, src)
  elseif rvt.dur_ex == :UKfixed 
    return uk_duration_fixed(m, r_ps, src)
  elseif (rvt.dur_ex == :BT14) & (rvt.dur_region == :ACR)
    return boore_thompson_2014(m, r_ps, src)
  elseif (rvt.dur_ex == :BT15)
    return boore_thompson_2015(m, r_ps, src, rvt)
  else
    if T <: Dual
      return T(NaN)
    elseif U <: Dual
      return U(NaN)
    else
      return S(NaN)
    end
  end
end

excitation_duration(m, r_ps, fas::FourierParameters, rvt::RandomVibrationParameters) = excitation_duration(m, r_ps, fas.source, rvt)


# Definition of Boore & Thompson (2012) constant values to be used within subsequent functions
include("coefficients/PJSbooreThompson2012.jl")


"""
  boore_thompson_2012_coefs(idx_m::T, idx_r::T; region::Symbol=:ACR) where T<:Int

Base function to extract the coefficients of the Boore & Thompson (2012) rms duration model.
"""
function boore_thompson_2012_coefs(idx_m::T, idx_r::T; region::Symbol=:ACR) where {T<:Int}
  idx = (idx_r - 1) * num_m_ii_bt12 + idx_m
  if region == :SCR
    @inbounds c = coefs_ena_bt12[idx, 3:9]
    return c
  else
    @inbounds c = coefs_wna_bt12[idx, 3:9]
    return c
  end
end


"""
  boore_thompson_2012_base(η::S, c::Vector{T}, ζ::T=0.05) where {S<:Real,T<:Float64}

Base function to compute the Boore & Thompson (2012) rms duration model for known magnitude and distance.
"""
function boore_thompson_2012_base(η::S, c::Vector{T}, ζ::T=0.05) where {S<:Real,T<:Float64}
  @inbounds ηc3 = η^c[3]
  @inbounds ratio = (c[1] + c[2] * ((1.0 - ηc3) / (1.0 + ηc3))) * (1.0 + c[4] / (2π * ζ) * (η / (1.0 + c[5] * η^c[6]))^c[7])
  return ratio
end


"""
  boore_thompson_2012(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}

Boore & Thompson (2012) rms duration model. Also outputs the excitation duration (given that its required within the rms duration calculations).

"""
function boore_thompson_2012(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}
  # for magnitude and distance that don't match coefficient tables we need to use bilinear interpolation
  # get the excitation duration (as recommended by Boore & Thompson, 2012)
  Dex = excitation_duration(m, r_ps, src, rvt)
  # get the oscillator period
  T_n = period(sdof)
  ζ = sdof.ζ_n
  # define the η parameter as T_n/Dex
  η = T_n / Dex

  # impose limits on the magnitudes and distances
  m = (m < 4.0) ? 4.0 : m
  m = (m > 8.0) ? 8.0 : m
  r_ps = (r_ps < 2.0) ? 2.0 : r_ps
  r_ps = (r_ps > 1262.0) ? 1262.0 : r_ps

  # get the bounding indicies
  i_lo = findlast(m_ii_bt12 .<= m)
  i_hi = findfirst(m_ii_bt12 .>= m)
  j_lo = findlast(r_jj_bt12 .<= r_ps)
  j_hi = findfirst(r_jj_bt12 .>= r_ps)

  # get the region for the rms duration model
  region = rvt.dur_region

  # check for situations in which the coefficients are known for the given m,r
  if i_lo == i_hi # we have coefficients for this magnitude
    if j_lo == j_hi # we have coefficients for this distance
      c = boore_thompson_2012_coefs(i_lo, j_lo; region=region)
      Dratio = boore_thompson_2012_base(η, c, ζ)
    else # we need to interpolate the distance values only
      r_lo = r_jj_bt12[j_lo]
      r_hi = r_jj_bt12[j_hi]
      c_lo = boore_thompson_2012_coefs(i_lo, j_lo; region=region)
      c_hi = boore_thompson_2012_coefs(i_lo, j_hi; region=region)
      D_lo = boore_thompson_2012_base(η, c_lo, ζ)
      D_hi = boore_thompson_2012_base(η, c_hi, ζ)
      lnDratio = log(D_lo) + (r_ps - r_lo) / (r_hi - r_lo) * log(D_hi / D_lo)
      Dratio = exp(lnDratio)
    end
  else # we need to interpolate the magnitudes
    if j_lo == j_hi # we have coefficients for this distance
      m_lo = m_ii_bt12[i_lo]
      m_hi = m_ii_bt12[i_hi]
      c_lo = boore_thompson_2012_coefs(i_lo, j_lo; region=region)
      c_hi = boore_thompson_2012_coefs(i_hi, j_lo; region=region)
      D_lo = boore_thompson_2012_base(η, c_lo, ζ)
      D_hi = boore_thompson_2012_base(η, c_hi, ζ)
      lnDratio = log(D_lo) + (m - m_lo) / (m_hi - m_lo) * log(D_hi / D_lo)
      Dratio = exp(lnDratio)
    else # we need to interpolate for this magnitude and distance
      # create the combinations of (m,r) from these indices
      m_lo = m_ii_bt12[i_lo]
      m_hi = m_ii_bt12[i_hi]
      r_lo = r_jj_bt12[j_lo]
      r_hi = r_jj_bt12[j_hi]

      # get the corresponding coefficients
      c_ll = boore_thompson_2012_coefs(i_lo, j_lo; region=region)
      c_lh = boore_thompson_2012_coefs(i_lo, j_hi; region=region)
      c_hl = boore_thompson_2012_coefs(i_hi, j_lo; region=region)
      c_hh = boore_thompson_2012_coefs(i_hi, j_hi; region=region)

      # get the corresponding ratios
      Dr_ll = boore_thompson_2012_base(η, c_ll, ζ)
      Dr_lh = boore_thompson_2012_base(η, c_lh, ζ)
      Dr_hl = boore_thompson_2012_base(η, c_hl, ζ)
      Dr_hh = boore_thompson_2012_base(η, c_hh, ζ)

      # bilinear interpolation
      # used repeated linear interpolation as the fastest method
      Δm = m_hi - m_lo
      Δr = r_hi - r_lo
      fac = 1.0 / (Δm * Δr)
      lnDr_ll = log(Dr_ll)
      lnDr_lh = log(Dr_lh)
      lnDr_hl = log(Dr_hl)
      lnDr_hh = log(Dr_hh)
      lnDratio = fac * ((lnDr_ll * (r_hi - r_ps) + lnDr_lh * (r_ps - r_lo)) * (m_hi - m) + (lnDr_hl * (r_hi - r_ps) + lnDr_hh * (r_ps - r_lo)) * (m - m_lo))
      Dratio = exp(lnDratio)
    end
  end
  Drms = Dratio * Dex
  return (Drms, Dex, Dratio)
end

boore_thompson_2012(m, r_ps, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) = boore_thompson_2012(m, r_ps, fas.source, sdof, rvt)





# Definition of Boore & Thompson (2015) constant values to be used within subsequent functions
include("coefficients/PJSbooreThompson2015.jl")


"""
  boore_thompson_2015_coefs(idx_m::T, idx_r::T; region::Symbol=:ACR) where T<:Int

Base function to extract the coefficients of the Boore & Thompson (2015) rms duration model.
"""
function boore_thompson_2015_coefs(idx_m::T, idx_r::T; region::Symbol=:ACR) where {T<:Int}
  idx = (idx_r - 1) * num_m_ii_bt15 + idx_m
  if region == :SCR
    @inbounds c = coefs_ena_bt15[idx, 3:9]
    return c
  else
    @inbounds c = coefs_wna_bt15[idx, 3:9]
    return c
  end
end


"""
  boore_thompson_2015_base(η::S, c::Vector{T}, ζ::T=0.05) where {S<:Real,T<:Float64}

Base function to compute the Boore & Thompson (2015) rms duration model for known magnitude and distance.
"""
function boore_thompson_2015_base(η::S, c::Vector{T}, ζ::T=0.05) where {S<:Real,T<:Float64}
  @inbounds ηc3 = η^c[3]
  @inbounds ratio = (c[1] + c[2] * ((1.0 - ηc3) / (1.0 + ηc3))) * (1.0 + c[4] / (2π * ζ) * (η / (1.0 + c[5] * η^c[6]))^c[7])
  return ratio
end


"""
  boore_thompson_2015(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}

Boore & Thompson (2015) rms duration model. Also outputs the excitation duration (given that its required within the rms duration calculations).

"""
function boore_thompson_2015(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}
  # for magnitude and distance that don't match coefficient tables we need to use bilinear interpolation of the log Drms values
  # get the excitation duration (as recommended by Boore & Thompson, 2015)
  Dex = excitation_duration(m, r_ps, src, rvt)
  # get the oscillator period
  T_n = period(sdof)
  ζ = sdof.ζ_n
  # define the η parameter as T_n/Dex
  η = T_n / Dex

  # impose limits on the magnitudes and distances
  m = (m < 2.0) ? 2.0 : m
  m = (m > 8.0) ? 8.0 : m
  r_ps = (r_ps < 2.0) ? 2.0 : r_ps
  r_ps = (r_ps > 1262.0) ? 1262.0 : r_ps

  # get the bounding indicies
  i_lo = findlast(m_ii_bt15 .<= m)
  i_hi = findfirst(m_ii_bt15 .>= m)
  j_lo = findlast(r_jj_bt15 .<= r_ps)
  j_hi = findfirst(r_jj_bt15 .>= r_ps)

  # get the region for the rms duration model
  region = rvt.dur_region

  # check for situations in which the coefficients are known for the given m,r
  if i_lo == i_hi # we have coefficients for this magnitude
    if j_lo == j_hi # we have coefficients for this distance
      c = boore_thompson_2015_coefs(i_lo, j_lo; region=region)
      Dratio = boore_thompson_2015_base(η, c, ζ)
    else # we need to interpolate the distance values only
      r_lo = r_jj_bt15[j_lo]
      r_hi = r_jj_bt15[j_hi]
      c_lo = boore_thompson_2015_coefs(i_lo, j_lo; region=region)
      c_hi = boore_thompson_2015_coefs(i_lo, j_hi; region=region)
      D_lo = boore_thompson_2015_base(η, c_lo, ζ)
      D_hi = boore_thompson_2015_base(η, c_hi, ζ)
      lnDratio = log(D_lo) + (r_ps - r_lo) / (r_hi - r_lo) * log(D_hi / D_lo)
      Dratio = exp(lnDratio)
    end
  else # we need to interpolate the magnitudes
    if j_lo == j_hi # we have coefficients for this distance
      m_lo = m_ii_bt15[i_lo]
      m_hi = m_ii_bt15[i_hi]
      c_lo = boore_thompson_2015_coefs(i_lo, j_lo; region=region)
      c_hi = boore_thompson_2015_coefs(i_hi, j_lo; region=region)
      D_lo = boore_thompson_2015_base(η, c_lo, ζ)
      D_hi = boore_thompson_2015_base(η, c_hi, ζ)
      lnDratio = log(D_lo) + (m - m_lo) / (m_hi - m_lo) * log(D_hi / D_lo)
      Dratio = exp(lnDratio)
    else # we need to interpolate for this magnitude and distance
      # create the combinations of (m,r) from these indices
      m_lo = m_ii_bt15[i_lo]
      m_hi = m_ii_bt15[i_hi]
      r_lo = r_jj_bt15[j_lo]
      r_hi = r_jj_bt15[j_hi]

      # get the corresponding coefficients
      c_ll = boore_thompson_2015_coefs(i_lo, j_lo; region=region)
      c_lh = boore_thompson_2015_coefs(i_lo, j_hi; region=region)
      c_hl = boore_thompson_2015_coefs(i_hi, j_lo; region=region)
      c_hh = boore_thompson_2015_coefs(i_hi, j_hi; region=region)

      # get the corresponding ratios
      Dr_ll = boore_thompson_2015_base(η, c_ll, ζ)
      Dr_lh = boore_thompson_2015_base(η, c_lh, ζ)
      Dr_hl = boore_thompson_2015_base(η, c_hl, ζ)
      Dr_hh = boore_thompson_2015_base(η, c_hh, ζ)

      # bilinear interpolation
      # used repeated linear interpolation as the fastest method
      Δm = m_hi - m_lo
      Δr = r_hi - r_lo
      fac = 1.0 / (Δm * Δr)
      lnDr_ll = log(Dr_ll)
      lnDr_lh = log(Dr_lh)
      lnDr_hl = log(Dr_hl)
      lnDr_hh = log(Dr_hh)
      lnDratio = fac * ((lnDr_ll * (r_hi - r_ps) + lnDr_lh * (r_ps - r_lo)) * (m_hi - m) + (lnDr_hl * (r_hi - r_ps) + lnDr_hh * (r_ps - r_lo)) * (m - m_lo))
      Dratio = exp(lnDratio)
    end
  end
  Drms = Dratio * Dex
  return (Drms, Dex, Dratio)
end

boore_thompson_2015(m, r_ps, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) = boore_thompson_2015(m, r_ps, fas.source, sdof, rvt)



"""
  liu_pezeshk_1999(m, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}

Liu & Pezeshk (1999) rms duration model. Also outputs the excitation duration (given that its required within the rms duration calculations).

"""
function liu_pezeshk_1999(m, r_ps::T, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}
  # get the excitation duration
  Dex = excitation_duration(m, r_ps, fas, rvt)

  # oscillator duration
  Dosc = 1.0 / ( 2π * sdof.f_n * sdof.ζ_n )
  # gamma parameter
  γ = Dex / Dosc
  # spectral moments, required for use within definition of α parameter
  mi = spectral_moments([0,1,2], m, r_ps, fas, sdof)
  m0 = mi.m0
  m1 = mi.m1
  m2 = mi.m2
  # define α
  α = sqrt( 2π * (1.0 - ( m1^2 )/( m0*m2 )) )
  n = 2.0
  # RMS duration
  Drms = Dex + Dosc * (γ^n / (γ^n + α))
  Dratio = Drms / Dex
  return (Drms, Dex, Dratio)
end



"""
    rms_duration(m::S, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}

Returns a 3-tuple of (Drms, Dex, Dratio), using a switch on `rvt.dur_rms`.
Default `rvt` makes use of the `:BT14` model for excitation duration, `Dex`.
- `m` is magnitude
- `r_ps` is an equivalent point source distance
"""
function rms_duration(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {T<:Real}
  if (rvt.dur_rms == :BT15) && (rvt.pf_method == :DK80)
    return boore_thompson_2015(m, r_ps, src, sdof, rvt)
  elseif (rvt.dur_rms == :BT12) && (rvt.pf_method == :CL56)
    return boore_thompson_2012(m, r_ps, src, sdof, rvt)
  elseif (rvt.dur_rms == :LP99) && (rvt.pf_method == :CL56)
    return liu_pezeshk_1999(m, r_ps, src, sdof, rvt)
    # println("Cannot pass `src::SourceParameters` for `rvt.dur_rms == :LP99`; method requires full `fas::FourierParameters` input")
    # return (T(NaN), T(NaN), T(NaN))
  else
    println("Inconsistent combination of `rvt.dur_rms` and `rvt.pf_method`")
    return (T(NaN), T(NaN), T(NaN))
  end
end

function rms_duration(m, r_ps, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) 
  if (rvt.dur_rms == :LP99) && (rvt.pf_method == :CL56)
    return liu_pezeshk_1999(m, r_ps, fas, sdof, rvt)
  else
    return rms_duration(m, r_ps, fas.source, sdof, rvt)
  end
end
