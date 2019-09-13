# For defining the ODEs for describing all the molecular dynamics.
# Designed for use with the Sundials ODE solver.
# Args: xdot = vector, [dCa2+/dt, dBCa/dt, dECa/dt, dImCa/dt, dBCa*/dt]
#       x = current state of continuous variables, [Ca2+, BCa, ECa, ImCa, BCa*]
#       t = time
#       params =  Dict, parameters [calcium_baseline, ca_influx, ca_outflux, b_i, f_i, indicator,
#                 b_e, f_e, endogeneous, b_im, f_im, immobile,
#                 excitation_rate, release_rate, previous_BCa*]
# Returns: nothing

function molecularODE(t::Float64, x::Array{Float64}, xdot::Array{Float64}, baseline::Int64, 
                      calcium_rate::Float64, indicator::Int64, endogeneous::Int64, 
                      immobile::Int64, b_i::Float64, f_i::Float64, b_e::Float64, f_e::Float64, 
                      b_im::Float64, f_im::Float64, excitation::Float64, release::Float64)
  B, E, Im = indicator - (x[2] + x[5]), endogeneous - x[3], immobile - x[4]
  # d[Ca2+]/dt = b_i*[BCa] + b_e*[ECa] + b_im*x[ImCa] - f_i*[Ca2+]*[B] - f_e*[Ca2+]*[E] + f_im*x[Ca2+]*[Im] + \beta_{Ca}([Ca2+]_0 - [Ca2+])
  xdot[1] = b_i*x[2] + b_e*x[3] + b_im*x[4] - f_i*x[1]*B - f_e*x[1]*E - f_im*x[1]*Im + calcium_rate * (baseline - x[1])
  # d[BCa]/dt = f_i*[Ca2+]*[B] - b_i*[BCa] + release*[BCa*] - excitation*[BCa] 
  xdot[2] = f_i*x[1]*B - b_i*x[2] + release*x[5] - excitation*x[2] 
  # d[ECa]/dt = f_e*[Ca2+]*[E] - b_e*[ECa] 
  xdot[3] = f_e*x[1]*E - b_e*x[3] 
  # d[ImCa]/dt = f_im*[Ca2+]*[Im] - b_im*[ImCa] 
  xdot[4] = f_im*x[1]*Im - b_im*x[4] 
  # d[BCa*]/dt = excitation*[BCa] - release*[BCa*]
  xdot[5] = excitation*x[2] - release*x[5]
  return nothing
end
