# Convert the buffer forward rates from  1/(Molar sec) to 1/(molecule sec)
# Arguments:  f = Float64, the forward rate in 1/(Molar sec)
#             v = the volume in metres cubed
# Returns: Float64, the forward rate in per molecule per second

function perMolarToPerMolecule(f, v)
  N_av = 6.02214085774e+23
  v_litres = 1000 * v
  return f/(N_av * v_litres) # factor of 1000 to convert from litres to m^3
end
