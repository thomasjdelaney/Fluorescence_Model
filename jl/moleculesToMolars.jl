# For converting a number of molecules to molars, given the cell volume
# Arguments:  molecules = number of molecules
#             volume = volume in metres cubed
# Returns:    Float64, the molar density

function moleculesToMolars(molecules::Float64, volume::Float64)
  N_av = 6.02214085774e+23 # avagadro's number
  volume_litres = 1000 * volume
  molars = molecules/(N_av * volume_litres)
  return molars
end

function moleculesToMolars(molecules::Array{Float64, 1}, volume::Float64)
  N_av = 6.02214085774e+23 # avagadro's number
  volume_litres = 1000 * volume
  molars = molecules/(N_av * volume_litres)
  return molars
end

function moleculesToMolars(molecules::Array{Float64, 2}, volume::Float64)
  N_av = 6.02214085774e+23 # avagadro's number
  volume_litres = 1000 * volume
  molars = molecules/(N_av * volume_litres)
  return molars
end

