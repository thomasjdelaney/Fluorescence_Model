# Convert a number of molars and a volume to a number of molecules
# Arguments:  m = the number of molars
#             v = the volume in metres cubed
# Returns:    Int, the number of molecules

function molarsToMolecules(m, v)
  N_av = 6.02214085774e+23
  v_litres = 1000 * v
  num_molecules = try
    floor(Int, N_av * m * v_litres)
  catch
    floor(N_av * m * v_litres)
  end
  return num_molecules
end
