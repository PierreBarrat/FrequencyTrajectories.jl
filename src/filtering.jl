"""
	inbin!(traj::Trajectory, fb::FrequencyBin)

Determine whether `traj` goes through the bin `fb` at any time point.
Set the `timetofreq` field of `traj` accordingly.
"""
function inbin!(traj::Trajectory, fb::FrequencyBin)
	for (t,f) in traj
		if in(f, fb)
			traj.timetofreq[fb] = t
			return true
		end
	end
	return false
end
"""
	inbin(traj::Trajectory, fb::FrequencyBin)

Determine whether `traj` goes through the bin `fb` at any time point.
"""
function inbin(traj, fb)
	for (t,f) in traj
		if in(f, fb)
			return true
		end
	end
	return false
end

function filter!(trajectories::Vector{Trajectory}, fb::FrequencyBin)
	idx = findall(traj -> inbin!(traj, fb), trajectories)
	return trajectories[idx]
end
filter!(trajectories, f, df) = filter!(trajectories, FrequencyBin(f, df))

"""
	filter(f::Function, trajectories::Vector{Trajectory})
"""
function filter(f::Function, trajectories::Vector{Trajectory})
	idx = findall(f, trajectories)
	return trajectories[idx]
end
"""
	filter(trajectories::Vector{Trajectory}, fb::FrequencyBin)
"""
function filter(trajectories::Vector{Trajectory}, fb::FrequencyBin)
	return filter(traj -> inbin(traj, fb), trajectories)
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


