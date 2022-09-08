struct FrequencyBin
	f::Float64
	df::Float64
end

@with_kw mutable struct Trajectory{T}
	t::Vector{T} = zeros(T, 0)
	f::Vector{Float64} = zeros(Float64, 0)
	# some of this could be in a container
	pos::Int = 0. # position in the genome
	state::Any = 0.
	final_state::Symbol = :missing # :lost, :fixed, :missing
	# Other information
	data::Dict = Dict()
	#
	# index where freqbin is reached
	freq_at_bin::Dict{FrequencyBin, Float64} = Dict{FT.FrequencyBin, Float64}()
	time_at_bin::Dict{FrequencyBin, T}  = Dict{FT.FrequencyBin, eltype(t)}()
	ϕpos::Vector{Union{Missing, Float64}} = ones(Missing, 0) # contribution of the position to fitness
	ϕgenotype::Vector{Union{Missing, Float64}} = ones(Missing, 0) # average fitness of genotypes carying this mutation
end

fixed(T::Trajectory) = (T.final_state == :fixed)
fixes(T::Trajectory) = fixed(T)

duration(T::Trajectory) = T.t[end] - T.t[1]


