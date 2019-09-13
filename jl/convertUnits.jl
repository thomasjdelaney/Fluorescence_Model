"""
For converting the concentrations into molecules

Arguments:  cell_radius = the radius of the modelled cell (m)
            baseline,
            indicator,
            endogeneous,
            immobile,
            peak,
            f_i,
            f_e,
            f_im
Returns:    m_baseline, m_indicator, m_endogeneous, m_immobile, m_peak, m_f_i,
            m_f_e, m_f_im
"""
function convertUnits(cell_radius::Float64, baseline::Float64,  indicator::Float64, endogeneous::Float64,
                      immobile::Float64,    peak::Float64,      f_i::Float64,       f_e::Float64,
                      f_im::Float64)
  cell_volume = sphereVolume(cell_radius)
  m_baseline = molarsToMolecules(baseline, cell_volume)
  m_indicator = molarsToMolecules(indicator, cell_volume)
  m_endogeneous = molarsToMolecules(endogeneous, cell_volume)
  m_immobile = molarsToMolecules(immobile, cell_volume)
  m_peak = molarsToMolecules(peak, cell_volume)
  m_f_i = perMolarToPerMolecule(f_i, cell_volume)
  m_f_e = perMolarToPerMolecule(f_e, cell_volume)
  m_f_im = perMolarToPerMolecule(f_im, cell_volume)
  return cell_volume, m_baseline, m_indicator, m_endogeneous, m_immobile, m_peak, m_f_i, m_f_e, m_f_im
end

