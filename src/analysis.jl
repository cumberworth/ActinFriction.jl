"""
$(TYPEDSIGNATURES)

Calculate basic quantities.
"""
function calc_basic_quantities(lambda, times, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    x = p.deltas * lambda
    ls = 1 .+ p.deltas ./ p.deltad .* lambda
    ltot = ls .* overlaps
    R = p.Nsca * (p.Lf .- x) / (2pi)
    force_L_bending = bending_force(lambda, p)
    force_L_condensation = condensation_force(p)
    force_R_bending = force_L_to_R(force_L_bending, p)
    force_R_condensation = force_L_to_R(force_L_condensation, p)
    df = DataFrame(
        t=times,
        lmbda=lambda,
        ls=ls,
        ltot=ltot,
        x=x,
        R=R,
        force_L_bending=force_L_bending,
        force_L_condensation=force_L_condensation,
        force_R_bending=force_R_bending,
        force_R_condensation=force_R_condensation)

    return df
end

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

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-binding quasi-equilibrium.
"""
function calc_cX_quantities(lambda, times, p::RingParams)
    df = calc_basic_quantities(lambda, times, p)
    R_max = p.Nsca * p.Lf / (2pi)
    R_eq = calc_equilibrium_ring_radius(p)
    df[!, :R_eq_frac] = (R_max .- df.R) / (R_max - R_eq)
    df[!, :force_L_total] = df.force_L_bending .+ df.force_L_condensation
    df[!, :force_R_total] = df.force_R_bending .+ df.force_R_condensation
    df[!, :zeta_ring_cX] = friction_coefficient_ring_cX(lambda, p)
    df[!, :zeta_overlap_cX] = friction_coefficient_overlap_cX(lambda, p)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-binding quasi-equilibrium.

Ignore singly bound crosslinkers.
"""
function calc_cX_ignore_Ns_quantities(lambda, times, p::RingParams)
    df = calc_basic_quantities(lambda, times, p)
    R_max = p.Nsca * p.Lf / (2pi)
    R_eq = calc_equilibrium_ring_radius_ignore_Ns(p)
    df[!, :R_eq_frac] = (R_max .- df.R) / (R_max - R_eq)
    df[!, :force_L_total] = df.force_L_bending .+ df.force_L_condensation
    df[!, :force_R_total] = df.force_R_bending .+ df.force_R_condensation
    df[!, :zeta_ring_cX] = friction_coefficient_ring_cX(lambda, p)
    df[!, :zeta_overlap_cX] = friction_coefficient_overlap_cX(lambda, p)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium.
"""
function calc_Nd_quantities(lambda, Ndtot, times, p::RingParams)
    df = calc_basic_quantities(lambda, times, p)
    R_max = p.Nsca * p.Lf / (2pi)
    R_eq = calc_equilibrium_ring_radius_ignore_Ns(p)
    df[!, :R_eq_frac] = (R_max .- df.R) / (R_max - R_eq)
    df[!, :Ndtot] = Ndtot
    df[!, :Nd_occupancy] = Ndtot ./ df.ltot
    df[!, :force_L_entropy] = entropic_force(lambda, Ndtot, p)
    df[!, :force_R_entropy] = force_L_to_R(df.force_L_entropy, p)
    df[!, :force_L_total] = df.force_L_bending .+ df.force_L_entropy
    df[!, :force_R_total] = df.force_R_bending .+ df.force_R_entropy
    df[!, :zeta_ring_Nd_double_exp] = friction_coefficient_ring_Nd(lambda, Ndtot, p)
    overlaps = 2p.Nf - p.Nsca
    df[!, :zeta_overlap_Nd_double_exp] = friction_coefficient_overlap_Nd(lambda, Ndtot/overlaps, p)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium.
"""
function calc_Ns_quantities(lambda, Ndtot, Nstot, times, p::RingParams)
    df = calc_Nd_quantities(lambda, Ndtot, times, p)
    df[!, :Nstot] = Nstot

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium with discrete Nd.
"""
function calc_discrete_Nd_quantities(lambda, Ndtot, Nds, times, p::RingParams)
    df = calc_Nd_quantities(lambda, Ndtot, times, p)
    zeta_continuous_l = friction_coefficient_continuous_l_ring_Nd(lambda, Nds, p)
    zeta_continuous_l_ave = friction_coefficient_continuous_l_ave_Nd(lambda, Nds, p)
    zeta_continuous_l_Ndtot = friction_coefficient_continuous_l_Ndtot_ring_Nd(lambda, Nds, p)
    df[!, :zeta_continuous_l] = zeta_continuous_l
    df[!, :zeta_continuous_l_ave] = zeta_continuous_l_ave
    df[!, :zeta_continuous_l_Ndtot] = zeta_continuous_l_Ndtot

    return df
end
