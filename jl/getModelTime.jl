"""
For getting the time period required for the model, and the indices corresponding to spike times,
and the indices for resampling.

Arguments:  spike_train, 1s and 0s
            frequency, the sampling frequency
Returns:    model_time, The times at which the model will be solved
            spike_indices, the index of each spike, and the end of the trial
            model_indices, the index of each model time point
"""
function getModelTime(spike_train::Array{Int64,1}, frequency::Float64)
  num_time_points = length(spike_train)
  model_time = collect((1:num_time_points)/frequency)
  spike_indices = findall(spike_train .> 0)
  model_indices = collect(1:num_time_points)
  return model_time, spike_indices, model_indices
end

