function Base.iterate(traj::Trajectory, state=1)
	state > length(traj.f) ? nothing : ((traj.t[state], traj.f[state]), state + 1)
end

Base.eltype(::Type{Trajectory{T}}) where T = Tuple{T, Float64}
Base.length(T::Trajectory) = length(T.f)

Base.in(f, fb::FrequencyBin) = fb.f - fb.df < f < fb.f + fb.df
Base.isless(x::FrequencyBin, y::FrequencyBin) = isless(x.f, y.f)

@recipe function f(traj::Trajectory)
	label --> ""
	return traj.t, traj.f
end
@recipe function f(traj::Trajectory, fb::FrequencyBin)
	if haskey(traj.timetofreq, fb)
		if traj.final_state == :lost
			linecolor --> :blue
		elseif traj.final_state == :fixed
			linecolor --> :red
		else
			linecolor --> :black
		end
		label --> ""
		t0 = traj.timetofreq[fb]
		return traj.t .- t0, traj.f
	else
		label --> ""
		return ()
	end
end
