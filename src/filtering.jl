"""
	inbin!(traj::Trajectory, fb::FrequencyBin; always_below=false)

Determine whether `traj` goes through the bin `fb` at any time point.
Set the `time_at_bin` field of `traj` accordingly.
"""
function inbin!(traj::Trajectory, fb::FrequencyBin; always_below=false)
	for (t,f) in traj
		if always_below && f > fb.f + fb.df
			return false
		end

		if in(f, fb)
			traj.freq_at_bin[fb] = f
			traj.time_at_bin[fb] = t
			return true
		end
	end
	return false
end
"""
	inbin(traj::Trajectory, fb::FrequencyBin; always_below=false)

Determine whether `traj` goes through the bin `fb` at any time point.
"""
function inbin(traj, fb; always_below=false)
	for (t,f) in traj
		if always_below && f > fb.f + fb.df
			return false
		end

		if in(f, fb)
			return true
		end
	end
	return false
end

"""
	filter!(trajectories::AbstractVector{<:Trajectory}, fb::FrequencyBin; kwargs...)

`kwargs...` are the ones of `FT.inbin!`.
"""
function filter!(
	trajectories::AbstractVector{<:Trajectory}, fb::FrequencyBin;
	center_in_bin = true, kwargs...
)
	idx = if center_in_bin
		findall(traj -> inbin!(traj, fb; kwargs...), trajectories)
	else
		findall(traj -> inbin(traj, fb; kwargs...), trajectories)
	end
	return keepat!(trajectories, idx)
end
filter!(trajectories, f, df; kwargs...) = filter!(trajectories, FrequencyBin(f, df); kwargs...)


"""
	filter(trajectories::AbstractVector{Trajectory}, fb::FrequencyBin; kwargs...)
"""
function filter(trajectories::AbstractVector{<:Trajectory}, fb::FrequencyBin; kwargs...)
	filter(traj -> inbin(traj, fb; kwargs...), trajectories)
end

"""
	has_fitness_above(traj, ϕmin)

Does `traj` have a fitness *strictly* above `ϕmin` for all its duration?
"""
function has_fitness_above(traj::Trajectory, ϕmin)
	for ϕ in traj.ϕpos
		if ismissing(ϕ) || ϕ <= ϕmin
			return false
		end
	end
	return true
end
"""
	has_fitness_above(traj, ϕmin, fb::FrequencyBin)

Does `traj` have a fitness *strictly* above `ϕmin` when it enters `fb`?
Return `false` if `traj` never enters `fb`.
"""
function has_fitness_above(traj::Trajectory, ϕmin::Float64, fb::FrequencyBin)
	for (i, (t,f)) in enumerate(traj)
		if in(f, fb)
			return traj.ϕpos[i] > ϕmin
		end
	end
	return false
end


