module FrequencyTrajectories

using Parameters
using RecipesBase

import Base: iterate, eltype, length
import Base: in, isless, filter!

export FT
export FrequencyBin, Trajectory
export filter!, fixation_probability

include("objects.jl")
include("interfaces.jl")
include("methods.jl")
include("filtering.jl")

const FT = FrequencyTrajectories

end # module
