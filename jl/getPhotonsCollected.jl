"""
get the number of photons collected, uses a normal distribution
Arguments:  ebc, Float64, excited buffered calcium molecules
            release_rate, Float64, expected proportion ebc molecules that release a photon
            capture_rate, Float64, bernoulli probability of capturing each individual photon
Returns:    number_of_photons, Int64
"""
function getPhotonsCollected(ebc::Number, release_rate::Number, capture_rate::Number)
  0 <= capture_rate <= 1 || return 0
  released = floor(Int, ebc*release_rate)
  released <= 0 && return 0
  photon_dist = Binomial(released, capture_rate)
  return rand(photon_dist)
end

