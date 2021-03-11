
# function to implement generic geometric spreading (piecewise log-linear)
function geometric_spreading( r::Real, fas::Union{PJSfasParams,PJSfasParamsQr} )
    z_r = 1.0
    for i = 2:length(fas.Rrefi)
      @inbounds if r < fas.Rrefi[i]
        @inbounds z_r *= (fas.Rrefi[i-1] / r)^fas.γi[i-1]
        return z_r
      else
        @inbounds z_r *= (fas.Rrefi[i-1] / fas.Rrefi[i])^fas.γi[i-1]
      end
    end
    return z_r
end

function geometric_spreading( r::Real, fas::Union{PJSfasParamsGeo,PJSfasParamsGeoQr} )
    z_r = 1.0
    for i = 2:length(fas.Rrefi)
      	@inbounds if r < fas.Rrefi[i]
        	@inbounds z_r *= (fas.Rrefi[i-1] / r)^fas.γi[i-1]
        	return z_r
      	else
        	@inbounds z_r *= (fas.Rrefi[i-1] / fas.Rrefi[i])^fas.γi[i-1]
      	end
	  	@inbounds if r >= fas.Rrefi[end]
		  	z_r *= (fas.Rrefi[end] / r)^fas.γf
	  	end
    end
    return z_r
end

# try implementing the smooth CY14 transition here
function geometric_spreading_cy( r::Real, fas::Union{PJSfasParams,PJSfasParamsQr} )
    ln_z_r = -fas.γi[1]*log(r) + (-fas.γi[2]+fas.γi[1])*log(sqrt(r^2 + 50.0^2)) - (-fas.γi[2]+fas.γi[1])*log(sqrt(1.0^2 + 50.0^2))
	z_r = exp(ln_z_r)
	return z_r
end

function geometric_spreading_cy( r::Real, fas::Union{PJSfasParamsGeo,PJSfasParamsGeoQr} )
    ln_z_r = -fas.γi[1]*log(r) + (-fas.γf+fas.γi[1])*log(sqrt(r^2 + 50.0^2)) - (-fas.γf+fas.γi[1])*log(sqrt(1.0^2 + 50.0^2))
	z_r = exp(ln_z_r)
	return z_r
end



function finite_fault_factor(m::Real)
	Mt1 = 5.744
	Mt2 = 7.744
	if m <= Mt1
		a1 = 0.7497
		b1 = 0.4300
		h = a1 + b1*(m - Mt1)
	elseif m >= Mt2
		a2 = 1.4147
		b2 = 0.2350
		h = a2 + b2*(m - Mt2)
	else
		c0 = 0.7497
		c1 = 0.4300
		c2 = -0.04875
		c3 = 0.0
		h = c0 + c1*(m - Mt1) + c2*(m - Mt1)^2 + c3*(m - Mt1)^3
	end
	return 10.0^h
end
