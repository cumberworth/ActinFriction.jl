"""
$(TYPEDSIGNATURES)

Use linear interpolation to create mean of all columns in dataframes.
"""
function meanvar_dfs(dfs, interval=1.0)
    cols = length(names(dfs[1])) - 1
    ivs = [[] for _ in 1:cols]
    times = 0:interval:dfs[1].t[end]
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

function calc_velocity(x, t)
    vs = [0.0]
    for i in 2:length(x)
        v = (x[i] - x[i - 1]) / (t[i] - t[i - 1])
        push!(vs, v)
    end

    return vs
end

"""
$(TYPEDSIGNATURES)

Calculate basic quantities.
"""
function calc_basic_quantities(lambda, times, p::RingParams)
    x = p.deltas * lambda
    ls = 1 .+ p.deltas ./ p.deltad .* lambda
    ltot = ls .* overlaps(p)
    R = p.Nsca * (p.Lf .- x) / (2pi)
    force_L_bending = bending_force(lambda, p)
    force_R_bending = force_L_to_R(force_L_bending, p)
    dR_dt = calc_velocity(R, times)
    df = DataFrame(
        t=times,
        lmbda=lambda,
        ls=ls,
        ltot=ltot,
        x=x,
        R=R,
        force_L_bending=force_L_bending,
        force_R_bending=force_R_bending,
        dR_dt=dR_dt)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-binding quasi-equilibrium.
"""
function calc_cX_quantities(lambda, times, p::RingParams)
    df = calc_basic_quantities(lambda, times, p)
    R_max = p.Nsca * p.Lf / (2pi)
    R_eq = equilibrium_ring_radius(p)
    force_L_condensation = fill(condensation_force(p), size(times))
    force_R_condensation = force_L_to_R(force_L_condensation, p)
    df[!, :R_eq_frac] = (R_max .- df.R) / (R_max - R_eq)
    df[!, :force_L_condensation] = force_L_condensation
    df[!, :force_R_condensation] = force_R_condensation
    df[!, :force_L_total] = df.force_L_bending .+ df.force_L_condensation
    df[!, :force_R_total] = df.force_R_bending .+ df.force_R_condensation
    df[!, :zeta_cX] = friction_coefficient_cX(lambda, p)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium.
"""
function calc_Nd_base_quantities(lambda, Ndtot, times, p::RingParams)
    df = calc_basic_quantities(lambda, times, p)
    df[!, :Ndtot] = Ndtot
    df[!, :force_L_entropy] = entropic_force(lambda, Ndtot, p)
    df[!, :force_R_entropy] = force_L_to_R(df.force_L_entropy, p)
    df[!, :force_L_total] = df.force_L_bending .+ df.force_L_entropy
    df[!, :force_R_total] = df.force_R_bending .+ df.force_R_entropy
    df[!, :zeta_Nd_exp] = friction_coefficient_Nd_exp(Ndtot / overlaps(p), p)
    df[!, :Nd_occupancy] = Ndtot ./ df.ltot
    R_max = p.Nsca * p.Lf / (2pi)
    R_eq = equilibrium_ring_radius(p)
    df[!, :R_eq_frac] = (R_max .- df.R) / (R_max - R_eq)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium.
"""
function calc_Nd_exp_quantities(lambda, Ndtot, times, p::RingParams)
    df = calc_Nd_base_quantities(lambda, Ndtot, times, p)
    zeta = df[!, :zeta_Nd_exp]

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium with discrete Nd.
"""
function calc_Nd_exact_quantities(lambda, Ndtot, times, p::RingParams)
    df = calc_Nd_base_quantities(lambda, Ndtot, times, p)
    zeta_Nd_exact = friction_coefficient_Nd_exact(Nd, p)
    df[!, :zeta_Nd_exact] = zeta_Nd_exact

    return df
end
