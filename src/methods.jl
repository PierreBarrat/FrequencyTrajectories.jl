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
function fixation_probability(f, trajectories::AbstractVector{<:Trajectory})
	nlost = count(T -> f(T) && T.final_state == :lost, trajectories)
	nfix = count(T -> f(T) && T.final_state == :fixed, trajectories)
	return nfix / (nfix + nlost)
end

"""
	fixation_probability(trajectories::Vector{Trajectory}, fb::FrequencyBin)
	fixation_probability(trajectories, f, df)
"""
function fixation_probability(trajectories, fb::FrequencyBin)
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

##########################################################################################
######################################### Misc. ##########################################
##########################################################################################

function StatsBase.mean(
    trajectories::Vector{<:Trajectory}, fb::FrequencyBin;
    K=min(length(trajectories), 50), kwargs...,
)
	T = filter(x -> inbin!(x, fb; kwargs...), trajectories)
	isempty(T) && return Float64[], Float64[], Int[]

	# Constructing grid
	tmin = findmin(x -> x.t[1] - x.time_at_bin[fb], T)[1]
	tmax = findmax(x -> x.t[end] - x.time_at_bin[fb], T)[1]
	L = findmax(length, T)[1]

    grid_len = floor(K*L) |> Int
	tgrid = collect(range(tmin, tmax, grid_len))
	xgrid = zeros(Float64, grid_len)
	Zs = zeros(Int, grid_len)

	# Interpolating and computing values on grid
	for traj in T
		t0 = traj.time_at_bin[fb]
		itp = LinearInterpolation(traj.t .- t0, traj.f)
		for (i, t) in enumerate(tgrid)
			if traj.t[1] - t0 <= t <= traj.t[end] - t0
				xgrid[i] += itp(t)
			elseif t < traj.t[1] - t0
				xgrid[i] += 0.
			elseif t > traj.t[end] - t0
				xgrid[i] += if traj.final_state == :fixed
					1.
				elseif traj.final_state == :lost
					0.
				else
					traj.f[end]
				end
			end
			Zs[i] += 1
		end
	end

	return xgrid ./ Zs, tgrid, Zs
end

function StatsBase.mean(T::Vector{<:Trajectory})
	# Constructing grid
	tmin = findmin(x -> x.t[1], T)[1]
	tmax = findmax(x -> x.t[end], T)[1]
	L = findmax(length, T)[1]
	tgrid = range(tmin, tmax, L)
	xgrid = zeros(Float64, L)
	Zs = zeros(Int, L)

	# Interpolating and computing values on grid
	for traj in T
		itp = LinearInterpolation(traj.t, traj.f)
		for (i, t) in enumerate(tgrid)
			if traj.t[1] <= t <= traj.t[end]
				xgrid[i] += itp(t)
			elseif t < traj.t[1]
				xgrid[i] += 0.
			elseif t > traj.t[end]
				xgrid[i] += if traj.final_state == :fixed
					1.
				elseif traj.final_state == :lost
					0.
				else
					traj.f[end]
				end
			end
			Zs[i] += 1
		end
	end

	return xgrid ./ Zs, tgrid, Zs
end

function mean_freq_at_bin!(trajectories::Vector{Trajectory}, fb::FrequencyBin)
	f = 0.
	Z = 0
	for T in trajectories
		if inbin!(T, fb)
			f += T.freq_at_bin[fb]
			Z += 1
		end
	end
	return f/Z
end



duration(T::Trajectory) = T.t[end] - T.t[1]
