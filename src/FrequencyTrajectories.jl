module FrequencyTrajectories

using Chain
using DataFrames
using Interpolations
using Parameters
using RecipesBase
using StatsBase


import Base: iterate, eltype, length
import Base: in, isless, filter, filter!
import StatsBase: mean

export FT
export FrequencyBin, Trajectory
export filter!, fixation_probability, fixed, fixes, get_trajectories

include("objects.jl")
include("interfaces.jl")
include("methods.jl")
include("filtering.jl")
include("io.jl")

const FT = FrequencyTrajectories

end # module
