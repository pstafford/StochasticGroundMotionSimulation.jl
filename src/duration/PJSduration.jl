


# function to implement the Boore & Thompson (2014) excitation duration function
function boore_thompson_2014(m, r_ps::U, src::SourceParameters{S,T}) where {S<:Float64, T<:Real, U<:Real}
  # source duration
  fa, fb, ε = corner_frequency(m, src)
  if src.model == :Atkinson_Silva_2000
    Ds = 0.5 / fa + 0.5 / fb
  else
    Ds = 1.0 / fa
  end
  # path duration
  if r_ps == 0.0
    return Ds * oneunit(U)
  elseif r_ps > 0.0 && r_ps <= 7.0
    return Ds + r_ps / 7.0 * 2.4
  elseif r_ps > 7.0 && r_ps <= 45.0
    return Ds + 2.4 + (r_ps - 7.0)/(45.0 - 7.0)*(8.4 - 2.4)
  elseif r_ps > 45.0 && r_ps <= 125.0
    return Ds + 8.4 + (r_ps - 45.0)/(125.0 - 45.0)*(10.9 - 8.4)
  elseif r_ps > 125.0 && r_ps <= 175.0
    return Ds + 10.9 + (r_ps - 125.0)/(175.0 - 125.0)*(17.4 - 10.9)
  elseif r_ps > 175.0 && r_ps <= 270.0
    return Ds + 17.4 + (r_ps - 175.0)/(270.0 - 175.0)*(34.2 - 17.4)
  elseif r_ps > 270.0
    return Ds + 34.2 + 0.156 * (r_ps - 270.0)
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

boore_thompson_2014(m, r_ps, fas::FourierParameters) = boore_thompson_2014(m, r_ps, fas.source)


function excitation_duration(m, r_ps::U, src::SourceParameters{S,T}, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real,U<:Real}
  if rvt.dur_ex == :BT14
    return boore_thompson_2014(m, r_ps, src)
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

function boore_thompson_2012_coefs(idx_m::T, idx_r::T; region::Symbol=:WNA) where T<:Int
  idx = (idx_r - 1) * num_m_ii_bt12 + idx_m
  if region == :ENA
    @inbounds c = coefs_ena_bt12[idx,3:9]
    return c
  else
    @inbounds c = coefs_wna_bt12[idx,3:9]
    return c
  end
end


# base function for Boore & Thompson (2012)
function boore_thompson_2012_base(η::S, c::Vector{T}, ζ::T=0.05) where {S<:Real,T<:Float64}
  @inbounds ηc3 = η^c[3]
  @inbounds ratio = (c[1] + c[2]*((1.0 - ηc3)/(1.0 + ηc3))) * (1.0 + c[4]/(2π*ζ)*( η /(1.0 + c[5]*η^c[6]) )^c[7])
  return ratio
end


# function to implement the Boore & Thompson (2012) duration ratio model
function boore_thompson_2012(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real}
  # for magnitude and distance that don't match coefficient tables we need to use bilinear interpolation
  # get the excitation duration (as recommended by Boore & Thompson, 2012)
  Dex = boore_thompson_2014(m, r_ps, src)
  # get the oscillator period
  T_n = period(sdof)
  ζ = sdof.ζ_n
  # define the η parameter as T_n/Dex
  η = T_n / Dex

  # impose limits on the magnitudes and distances
  m = ( m < 4.0 ) ? 4.0 : m
  m = ( m > 8.0 ) ? 8.0 : m
  r_ps = ( r_ps < 2.0 ) ? 2.0 : r_ps
  r_ps = ( r_ps > 1262.0 ) ? 1262.0 : r_ps

  # get the bounding indicies
  i_lo = findlast( m_ii_bt12 .<= m )
  i_hi = findfirst( m_ii_bt12 .>= m )
  j_lo = findlast( r_jj_bt12 .<= r_ps )
  j_hi = findfirst( r_jj_bt12 .>= r_ps )

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
      lnDratio = log(D_lo) + (r_ps - r_lo)/(r_hi - r_lo)*log(D_hi/D_lo)
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
      lnDratio = log(D_lo) + (m - m_lo)/(m_hi - m_lo)*log(D_hi/D_lo)
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
      D = [ log(Dr_ll) log(Dr_lh); log(Dr_hl) log(Dr_hh) ]
      knots = ([m_lo, m_hi], [r_lo, r_hi])
      # create interpolant
      itp_D = interpolate(knots, D, Gridded(Linear()))
      lnDratio = itp_D(m,r_ps)
      Dratio = exp(lnDratio)
    end
  end
  Drms = Dratio * Dex
  return (Drms, Dex, Dratio)
end

boore_thompson_2012(m, r_ps, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) = boore_thompson_2012(m, r_ps, fas.source, sdof, rvt)





# Definition of Boore & Thompson (2015) constant values to be used within subsequent functions
include("coefficients/PJSbooreThompson2015.jl")

function boore_thompson_2015_coefs(idx_m::T, idx_r::T; region::Symbol=:WNA) where T<:Int
  idx = (idx_r - 1) * num_m_ii_bt15 + idx_m
  if region == :ENA
    @inbounds c = coefs_ena_bt15[idx,3:9]
    return c
  else
    @inbounds c = coefs_wna_bt15[idx,3:9]
    return c
  end
end

# function get_mr(idx_m::T, idx_r::T; region::Symbol=:WNA) where T<:Int
#   idx = (idx_r - 1) * num_m_ii_bt15 + idx_m
#   if region == :ENA
#     return ( coefs_ena_bt15[idx,1], coefs_ena_bt15[idx,2] )
#   else
#     return ( coefs_wna_bt15[idx,1], coefs_wna_bt15[idx,2] )
#   end
# end

# base function for Boore & Thompson (2015) - Note this is the same function as in BT12
function boore_thompson_2015_base(η::S, c::Vector{T}, ζ::T=0.05) where {S<:Real,T<:Float64}
  @inbounds ηc3 = η^c[3]
  @inbounds ratio = (c[1] + c[2]*((1.0 - ηc3)/(1.0 + ηc3))) * (1.0 + c[4]/(2π*ζ)*( η /(1.0 + c[5]*η^c[6]) )^c[7])
  return ratio
end


# function to implement the Boore & Thompson (2015) duration ratio model
function boore_thompson_2015(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real}
  # for magnitude and distance that don't match coefficient tables we need to use bilinear interpolation of the log Drms values
  # get the excitation duration (as recommended by Boore & Thompson, 2015)
  Dex = boore_thompson_2014(m, r_ps, src)
  # get the oscillator period
  T_n = period(sdof)
  ζ = sdof.ζ_n
  # define the η parameter as T_n/Dex
  η = T_n / Dex

  # impose limits on the magnitudes and distances
  m = ( m < 2.0 ) ? 2.0 : m
  m = ( m > 8.0 ) ? 8.0 : m
  r_ps = ( r_ps < 2.0 ) ? 2.0 : r_ps
  r_ps = ( r_ps > 1262.0 ) ? 1262.0 : r_ps

  # get the bounding indicies
  i_lo = findlast( m_ii_bt15 .<= m )
  i_hi = findfirst( m_ii_bt15 .>= m )
  j_lo = findlast( r_jj_bt15 .<= r_ps )
  j_hi = findfirst( r_jj_bt15 .>= r_ps )

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
      lnDratio = log(D_lo) + (r_ps - r_lo)/(r_hi - r_lo)*log(D_hi/D_lo)
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
      lnDratio = log(D_lo) + (m - m_lo)/(m_hi - m_lo)*log(D_hi/D_lo)
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
      D = [ log(Dr_ll) log(Dr_lh); log(Dr_hl) log(Dr_hh) ]
      knots = ([m_lo, m_hi], [r_lo, r_hi])
      # create interpolant
      itp_D = interpolate(knots, D, Gridded(Linear()))
      lnDratio = itp_D(m,r_ps)
      Dratio = exp(lnDratio)
    end
  end
  Drms = Dratio * Dex
  return (Drms, Dex, Dratio)
end

boore_thompson_2015(m, r_ps, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) = boore_thompson_2015(m, r_ps, fas.source, sdof, rvt)




"""
  rms_duration(m::S, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real}

Returns a 3-tuple of (Drms, Dex, Dratio), using a switch on `rvt.dur_rms`.
Default `:BT12` makes use of the `:BT14` model for excitation duration, `Dex`.
- `m` is magnitude
- `r_ps` is an equivalent point source distance
"""
function rms_duration(m, r_ps::T, src::SourceParameters, sdof::Oscillator, rvt::RandomVibrationParameters) where {S<:Float64,T<:Real}
  if rvt.dur_rms == :BT15
    return boore_thompson_2015(m, r_ps, src, sdof, rvt)
  elseif rvt.dur_rms == :BT12
    return boore_thompson_2012(m, r_ps, src, sdof, rvt)
  else
    return (T(NaN), T(NaN), T(NaN))
  end
end

rms_duration(m, r_ps, fas::FourierParameters, sdof::Oscillator, rvt::RandomVibrationParameters) = rms_duration(m, r_ps, fas.source, sdof, rvt)
