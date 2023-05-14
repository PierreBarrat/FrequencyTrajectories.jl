function trajectory_to_dataframe(traj::Trajectory, id)
	return DataFrame(
		t = traj.t,
		f = traj.f,
		pos = repeat([traj.pos], length(traj)),
		state = traj.state,
		final_state = traj.final_state,
		id = repeat([id], length(traj)),
	)
end

function trajectories_to_dataframe(trajectories)
	X = mapreduce(vcat, enumerate(trajectories)) do (i, T)
		trajectory_to_dataframe(T, i)
	end
	return X
end

"""
	trajectories_from_dataframe(df::DataFrame)

Read a trajectories from a dataframe built using `trajectories_to_dataframe`.
"""
function trajectories_from_dataframe(df::DataFrame)
	@assert names(df) == ["t", "f", "pos", "state", "final_state", "id"]
	return [_trajectory_from_df(X) for X in groupby(df, [:id])]
end

function _trajectory_from_df(df)
	@assert allequal(df.id)
	return Trajectory(;
		t = collect(df.t),
		f = collect(df.f),
		pos = first(df.pos),
		state = first(df.state),
		final_state = Symbol(first(df.final_state)),
	)
end
