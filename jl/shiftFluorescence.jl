"""
For shifting the modelled trace to be non-negative.
Shift the modelled trace to ensure no negative values.

Arguments:  modelled = the modelled trace.
Returns:    Array{Float64, 1}, zscored and non-negative fluorescence trace
"""

function shiftFluorescence(modelled::Array{Int})
  z_modelled = zscore(modelled)
  adjustment = abs(minimum([0.0, minimum(z_modelled)]))
  return z_modelled .+ adjustment
end
