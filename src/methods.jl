##########################################################################################
################################## Computing trajectories ################################
##########################################################################################

"""
	function get_trajectories(
		X;
		state = 0,
		fitness = missing,
		threshold = 0.05,
		time_scaling = 1,
	)

Assume rows of `X` are frequencies and time goes along columns, *i.e.* `X[:,1]` is the frequency series of the state of position `1`.
Also assume that `X[2*i-1,:]` (resp. `X[2*i,:]`) is the frequency of state `1` (resp. `-1`) at position `i`.
Same for `fitness` if provided.
"""
function get_trajectories(
	X;
	fitness = missing,
	threshold = 0.05,
	time_scaling = 1,
)
	@assert ismissing(fitness) || size(fitness) == size(X) "fitness and frequencies have inconsistent sizes"
	trajectories = Vector{Trajectory}(undef, 0)
	for i in 1:size(X, 2)
		pos = Int(ceil(i/2))
		state = iseven(i) ? -1 : 1
		append!(
			trajectories,
			Trajectory(
				X[:,i];
				pos,
				state,
				fitness = ismissing(fitness) ? missing : fitness[:, i],
				threshold,
				time_scaling,
			)
		)
	end
	return trajectories
end

"""
	function Trajectory(
		X::Vector{Float64};
		pos = 0,
		state = 0,
		fitness = nothing,
		threshold = 0.05,
		time_scaling = 1,
	)

## Arguments
- `X`: frequency vector

"""
function Trajectory(
	X::Vector{Float64};
	pos = 0,
	state = 0,
	fitness = missing,
	threshold = 0.05,
	time_scaling = 1,
)
	indices = get_traj_index(X, threshold)
	# Assert non overlapping trajectories
	for i in 1:length(indices)-1
		@assert isempty(intersect(indices[i], indices[i+1])) "Overlapping trajectories. Frequency series $X"
	end

	trajectories = Trajectory[]
	for t in indices
		final_state = if X[t[end]] < threshold
			:lost
		elseif X[t[end]] > 1-threshold
			:fixed
		else
			:missing
		end

		traj = Trajectory(;
			t = t * time_scaling,
			f = copy(X[t]),
			pos = pos,
			state = state,
			final_state = final_state,
			ϕpos = ismissing(fitness) ? ones(Missing, length(t)) : copy(fitness[t])
		)
		push!(trajectories, traj)
	end

	return trajectories
end

function get_traj_index(X, threshold)
	rx = Union{Missing, Float64}[missing, missing, missing] # last 3 values in X
	indices = []
	idx = Int64[]
	running = false
	for (i,x) in enumerate(X)
		shift!(rx, x)
		if is_start(rx, threshold, i)
			@assert !running
			running = true
		end

		running && push!(idx, i-1)

		if running && is_stop(rx, threshold, i, length(X))
			running = false
			push!(indices, idx)
			idx = Int64[]
		end
	end
	return indices
end

function is_start(rx, threshold, i)
	if i == 2
		# Special case: rx = [missing, val1, val2]
		return rx[2] < threshold && rx[3] >= threshold
	elseif i > 2
		return rx[1] < threshold && rx[2] < threshold && rx[3] >= threshold
	end
	return false
end
function is_stop(rx, threshold, i, L)
	if i == L
		return true
	elseif i > 2
		if rx[2] < threshold && rx[3] < threshold
			return true
		elseif rx[2] > 1-threshold && rx[3] > 1-threshold
			return true
		end
	end
	return false
end


function shift!(rx, x)
	rx[1] = rx[2]
	rx[2] = rx[3]
	rx[3] = x
	return nothing
end

##########################################################################################
################################## Fixation probability ##################################
##########################################################################################

"""
	fixation_probability(f::Function, trajectories::Vector{Trajectory})

Compute fraction of trajectories `T` for which `f(T)` is true and that fix.
"""
function fixation_probability(f::Function, trajectories::Vector{Trajectory})
	nlost = count(T -> f(T) && T.final_state == :lost, trajectories)
	nfix = count(T -> f(T) && T.final_state == :fixed, trajectories)
	return nfix / (nfix + nlost)
end

"""
	fixation_probability(trajectories::Vector{Trajectory}, fb::FrequencyBin)
	fixation_probability(trajectories, f, df)
"""
function fixation_probability(trajectories::Vector{Trajectory}, fb::FrequencyBin)
	return fixation_probability(T -> inbin(T, fb), trajectories)
end
fixation_probability(trajectories, f, df) = fixation_probability(trajectories, FrequencyBin(f, df))

"""
	fixation_probability(trajectories, fb::FrequencyBin, ϕmin)
	fixation_probability(trajectories, f, df, ϕmin)
"""
function fixation_probability(trajectories, fb::FrequencyBin, ϕmin)
	return fixation_probability(T -> has_fitness_above(T, ϕmin, fb), trajectories)
end
function fixation_probability(trajectories, f, df, ϕmin)
	return fixation_probability(trajectories, FrequencyBin(f, df), ϕmin)
end
