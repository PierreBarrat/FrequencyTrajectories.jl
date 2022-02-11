struct FrequencyBin
	f :: Float64
	df :: Float64
end

@with_kw mutable struct Trajectory{T}
	t :: Vector{T} = zeros(T, 0)
	f :: Vector{Float64} = zeros(Float64, 0)
	#
	pos :: Int = 0. # position in the genome
	state :: Int = 0.
	final_state :: Symbol = :missing # :lost, :fixed, :missing
	#
	# index where freqbin is reached
	timetofreq :: Dict{FrequencyBin, T}  = Dict{ST.FrequencyBin, eltype(t)}()
	ϕpos :: Vector{Union{Missing, Float64}} = ones(Missing, 0) # contribution of the position to fitness
	ϕgenotype :: Vector{Union{Missing, Float64}} = ones(Missing, 0) # average fitness of genotypes carying this mutation
end


