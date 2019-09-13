"""
For running the fluorescence model. Uses the 'molecularODE' fuction with the Sundials ODE solver
to run the model.
Arguments:  spike_times = the frames when a spike occurs
            cell_radius = the modelled radius of the cell (in metres)
            baseline = the baseline amount of free calcium (M)
            calcium_rate = The rate at which calcium is kicked out of the cell (s^-1)
            indicator = the concentration of indicator in the cell (M)
            endogeneous = the concentration of endogeneous mobile buffer in the cell (M)
            immobile = the concentration of endogeneous immobile buffer in the cell (M)
            b_i = the backward rate of the indicator binding (s^-1)
            f_i = the forward rate of the indicator binding (M^-1 s^-1)
            b_e = the backward rate of the endogeneous binding (s^-1)
            f_e = the forward rate of the endogeneous binding (M^-1 s^-1)
            b_im = the backward rate of the immobile binding (s^-1)
            f_im = the forward rate of the immobile binding (M^-1 s^-1)
            excitation = the proportion of BCa that gets excited at each time step (s^-1)
            release = the proportion of BCa* that releases a photon at each time step (s^-1)
            peak = the amount of free calcium in the cell after a spike (M)
            frequency = sampling frequency (s^-1)
            capture_rate = the probability of capturing a given photon
Returns:    fluorescence = Array{Float64, 1}, the modelled fluorescence
            spiking_sim = T x 5 Array{Float, 2}, the full model
            spiking_sim_time = Array{Float, 1} length T, the time points
"""
function calciumFluorescenceModel(spike_train::Array{Int64,1};  cell_radius::Float64=1e-5,
                                  baseline::Float64=0.045e-6,   calcium_rate::Float64=0.001,  
                                  indicator::Float64=1e-4,      endogeneous::Float64=1e-4,  
                                  immobile::Float64=78.7e-6,    b_i::Float64=160.0,         f_i::Float64=7.76699e8,
                                  b_e::Float64=10000.0,         f_e::Float64=100e6,         b_im::Float64=524.0,
                                  f_im::Float64=2.47e8,         excitation::Float64=0.15,   release::Float64=0.11,
                                  peak::Float64=2.9e-7,         frequency::Float64=100.0,  
                                  capture_rate::Float64=0.62)
  cell_volume, m_baseline, m_indicator, m_endogeneous, m_immobile, m_peak, m_f_i, m_f_e, m_f_im = convertUnits(cell_radius::Float64, baseline::Float64, indicator::Float64, endogeneous::Float64, immobile::Float64, peak::Float64, f_i::Float64, f_e::Float64, f_im::Float64)
  model_time, spike_indices, model_indices = getModelTime(spike_train, frequency)
  # baseline dynamics for 25000 seconds to reach baseline conditions
  # The ODEs are solved for every 100 seconds
  baseline_time = collect(0:100.0:25000)
  baseline_initial_conditions = [m_baseline, 0.0, 0.0, 0.0, 0.0]
  baseline_sim = Sundials.cvode((t, x, xdot) -> molecularODE(t, x, xdot, m_baseline, calcium_rate, m_indicator, m_endogeneous, m_immobile, b_i, m_f_i, b_e, m_f_e, b_im, m_f_im, excitation, release), baseline_initial_conditions, baseline_time)
  # pre-spiking dynamics, ODEs are solved for every millisecond
  if spike_indices[1] > 1 # if there is no spike in the first time bin
    pre_spike_time = model_time[1:spike_indices[1]-1] 
    pre_spike_initial_conditions = baseline_sim[end,:]
    pre_spike_sim = Sundials.cvode((t, x, xdot) -> molecularODE(t, x, xdot, m_baseline, calcium_rate, m_indicator, m_endogeneous, m_immobile, b_i, m_f_i, b_e, m_f_e, b_im, m_f_im, excitation, release), pre_spike_initial_conditions, pre_spike_time)
    # spiking dynamics, discontinuous jump at each spike
    spiking_sim = pre_spike_sim
    spiking_sim_time = pre_spike_time
    for sp in 1:length(spike_indices)-1
      spike_time = model_time[spike_indices[sp]:spike_indices[sp+1]-1]
      spiking_initial_conditions = [m_peak; spiking_sim[end, 2:end]]
      spike_sim = Sundials.cvode((t, x, xdot) -> molecularODE(t, x, xdot, m_baseline, calcium_rate, m_indicator, m_endogeneous, m_immobile, b_i, m_f_i, b_e, m_f_e, b_im, m_f_im, excitation, release), spiking_initial_conditions, spike_time)
      spiking_sim = [spiking_sim; spike_sim]
      spiking_sim_time = [spiking_sim_time; spike_time]
    end
  else # if there is a spike in the first time bin
    spiking_sim = zeros(Float64,(0,5))
    spiking_sim_time = zeros(Float64, 0)
    spiking_initial_conditions = [m_peak; baseline_sim[end, 2:end]]
    for sp in 1:length(spike_indices)-1
      spike_time = model_time[spike_indices[sp]:spike_indices[sp+1]-1]
      spike_sim = Sundials.cvode((t, x, xdot) -> molecularODE(t, x, xdot, m_baseline, calcium_rate, m_indicator, m_endogeneous, m_immobile, b_i, m_f_i, b_e, m_f_e, b_im, m_f_im, excitation, release), spiking_initial_conditions, spike_time)
      spiking_sim = [spiking_sim; spike_sim]
      spiking_sim_time = [spiking_sim_time; spike_time]
      spiking_initial_conditions = [m_peak; spiking_sim[end, 2:end]]
    end
  end
  # dynamics for the last spike until the end of the trail
  last_spike_time = model_time[spike_indices[end]:end]
  last_spike_initial_conditions = [m_peak; spiking_sim[end, 2:end]]
  last_spike_sim = Sundials.cvode((t, x, xdot) -> molecularODE(t, x, xdot, m_baseline, calcium_rate, m_indicator, m_endogeneous, m_immobile, b_i, m_f_i, b_e, m_f_e, b_im, m_f_im, excitation, release), last_spike_initial_conditions, last_spike_time)
  spiking_sim = [spiking_sim; last_spike_sim]
  spiking_sim_time = [spiking_sim_time; last_spike_time]
  # photon capturing
  fluorescence = [getPhotonsCollected(ebc, release, capture_rate) for ebc in spiking_sim[model_indices, end]]
  fluorescence = shiftFluorescence(fluorescence)
  spiking_sim_molars = moleculesToMolars(spiking_sim, cell_volume)
  return fluorescence, spiking_sim_molars[model_indices,:], spiking_sim_time[model_indices]
end
