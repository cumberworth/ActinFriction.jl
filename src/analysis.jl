"""
$(TYPEDSIGNATURES)

Use liner interpolation to create mean of all columns in dataframes.
"""
function meanvar_dfs(dfs, interval=1)
    cols = length(names(dfs[1])) - 1
    ivs = [[] for _ in 1:cols]
    times = 1:interval:dfs[1].t[end]
    df_means = DataFrame()
    df_vars = DataFrame()
    df_means.t = times
    df_vars.t = times
    for df in dfs
        unique!(df)
        for (i, n) in enumerate(names(df)[2:end])
            Interpolations.deduplicate_knots!(df.t, move_knots=true)
            interp = LinearInterpolation(df.t, df[!, n])
            push!(ivs[i], interp(times))
        end
    end
    for (i, n) in enumerate(names(dfs[1])[2:end])
        df_means[!, n] = mean(ivs[i])
        df_vars[!, n] = var(ivs[i])
    end

    return (df_means, df_vars)
end
