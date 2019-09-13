# module for holding function common to both pdmp dynamics and deterministic dynamics

module FluorescenceModel

using Distributions
using StatsBase
using Sundials

export calciumFluorescenceModel,
  convertUnits,
  getModelTime,
  getPhotonsCollected,
  molarsToMolecules,
  molecularODE,
  moleculesToMolars,
  perMolarToPerMolecule,
  shiftFluorescence,
  sphereVolume

include("calciumFluorescenceModel.jl")
include("convertUnits.jl")
include("getModelTime.jl")
include("getPhotonsCollected.jl")
include("molarsToMolecules.jl")
include("molecularODE.jl")
include("moleculesToMolars.jl")
include("perMolarToPerMolecule.jl")
include("shiftFluorescence.jl")
include("sphereVolume.jl")

end
