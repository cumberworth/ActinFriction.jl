function zero_overlap(u, t, integrator)
    u[1]
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker-binding quasi-equilibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_cX!(du, u, p, t)
    lambda = u[1]

    zeta = friction_coefficient_cX(lambda, p)
    forcetot = bending_force(lambda, p) + condensation_force(p)

    du[1] = forcetot / (zeta * p.deltas * (2p.Nf - p.Nsca))

    return nothing
end

function noise_ring_cX!(du, u, p, t)
    lambda = u[1]

    zeta = friction_coefficient_cX(lambda, p)
    du[1] = sqrt(2 * kb * p.T / zeta) / p.deltas
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker-binding quasi-equilibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_cX_filaments!(du, u, p, t)
    # make thing to check if pair has been calculated already?
    lambda = u[1]
    R = lambda_to_R(lambda, p)
    for Nf_i in 1:p.Nf
        directions = []
        overlaps = []
        # force on filament i, but just need to know number of overlaps
        # friction on filament i, need to know total amount of overlap
        for Nf_j in 1:p.Nf
            if filaments_adjacent(Nf_i, Nf_j)
                # direction for 0 and pi, also for Nsca = 2
                phi_ji = u[1 + j] - u[1 + i]
                phi_ji_abs = abs(phi_ji)
                direction = sign(phi_ji)
                if phi_ji_abs > pi
                    direction *= -1
                end
                phi_ji_min = mod(phi_ji_abs, pi)
                L_i = p.Lf - R * phi_ji_min
                overlaps, Lij = filament_overlap(Nf_i, Nf_j)
            end
            # check if overlapping by looking at z index
            # add to number of overlaps
            # add to total overlap
        end
        # use collected values to calculate total force and friction on filament i
        # update filament i position
    end

    # calculate total force on radius, and friction
    # update radius
    zeta = friction_coefficient_cX(u[1], p)
    forcetot = bending_force(u[1], p) + condensation_force(p)

    du[1] = forcetot / (zeta * p.deltas * (2p.Nf - p.Nsca))

    return nothing
end

function equation_of_motion_ring_Nd_update!(du, u, p, zeta)
    lambda = u[1]
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(lambda, p) + entropic_force(lambda, Ndtot, p)
    ltot = lambda_to_l(lambda, p) * overlaps
    A = 2 * p.r10 * p.r12 / ((p.r01 + p.r10)^2 + p.r10 * p.r12)
    B = p.r01 * ltot + p.r21 * Ndtot - p.r01 * Ndtot
    C = 2 * p.r21 * Ndtot
    du[1] = forcetot / (zeta * p.deltas * overlaps)
    du[2] = B * A - C

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker diffusion quasi-equlibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_exp!(du, u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exp(Nd, p)
    equation_of_motion_ring_Nd_update!(du, u, p, zeta)

    return nothing
end

function noise_ring_Nd_exp!(du, u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exp(Nd, p)
    du[1] = sqrt(2 * kb * p.T / zeta) / p.deltas
end

function equation_of_motion_ring_Nd_exact_base!(zeta, Ndtot, du, u, p, t)
    lambda = u[1]
    du .= 0

    # println("Ndtot: $Ndtot, lambda: $lambda")
    # l = lambda_to_l(lambda, p) * overlaps
    # bendingforce = bending_force(lambda, p)
    # entropicforce = entropic_force(lambda, Ndtot, p)
    # println("Total sites: $l, crosslinkers: $Ndtot, lambda: $lambda, Bending force: $bendingforce, Entropic force: $entropicforce")
    forcetot = bending_force(lambda, p) + entropic_force(lambda, Ndtot, p)
    overlaps = 2p.Nf - p.Nsca

    du[1] = forcetot / (zeta * p.deltas * overlaps)

    return nothing
end

function equation_of_motion_ring_Nd_exact_static_base!(zeta, Ndtot, du, u, p, t)
    du .= 0

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This tracks the number of crosslinkers per overlap, and the
number of bound crosslinkers per overlap to be updated with a jump process separately. This
is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_exact_Nds!(du, u, p, t)
    zeta = friction_coefficient_Nd_exact(u[3:end], p)
    Ndtot = u[2]
    equation_of_motion_ring_Nd_exact_base!(zeta, Ndtot, du, u, p, t)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_exact_Ndtot!(du, u, p, t)
    Nd = u[2]

    overlaps = 2p.Nf - p.Nsca
    Ndtot = Nd * overlaps
    zeta = friction_coefficient_Nd_exact(Nd, p)
    equation_of_motion_ring_Nd_exact_base!(zeta, Ndtot, du, u, p, t)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_constant_discrete_Ndtot!(du, u, p, t)
    Nd = u[2]

    overlaps = 2p.Nf - p.Nsca
    Ndtot = Nd * overlaps
    equation_of_motion_ring_Nd_exact_base!(p.zeta, Ndtot, du, u, p, t)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_exact_Ndtot_static!(du, u, p, t)
    Nd = u[2]

    overlaps = 2p.Nf - p.Nsca
    Ndtot = Nd * overlaps
    zeta = friction_coefficient_Nd_exact(Nd, p)
    equation_of_motion_ring_Nd_exact_static_base!(zeta, Ndtot, du, u, p, t)

    return nothing
end

function noise_ring_Nd_exact_Ndtot!(du, u, p, t)
    Nd = u[2]

    zeta = friction_coefficient_Nd_exact(Nd, p)
    du[1] = sqrt(2 * kb * p.T / zeta) / p.deltas
end

"""
$(TYPEDSIGNATURES)

Generate crosslinker binding rate function.
"""
function binding_rate_generator(i::Integer)
    function binding_rate(u, p, t)
        lambda = u[1]

        # println("Binding rate generator lambda = $lambda")
        l = lambda_to_l(u[1], p)
        l_discrete = lambda_to_l_discrete(u[1], p)
        # l = lambda_to_l(u[1], p)
        if l_discrete == u[i]
        # if l <= u[i]
            # println("No binding attempt due to max occupancy")
            return 0.0
        else
            A = 2 * p.r10 * p.r12 / ((p.r01 + p.r10)^2 + p.r10 * p.r12)
            B = p.r01 * l + p.r21 * u[i] - p.r01 * u[i]
            return A * B 
        end
    end

    return binding_rate
end

"""
$(TYPEDSIGNATURES)

Generate crosslinker unbinding rate function.
"""
function unbinding_rate_generator(i::Integer)
    function unbinding_rate(u, p, t)
        # lambda = u[1]
        # l = lambda_to_l_discrete(u[1], p)
        # println("Unbinding rate generator lambda = $lambda")
        if u[i] == 1.0
            # println("No unbinding attempt due to min occupancy")
            return 0.0
        else
            return 2 * p.r21 * u[i]
        end
    end

    return unbinding_rate
end

"""
$(TYPEDSIGNATURES)

Generate reaction event function.
"""
function reaction_generator(i::Integer, event::Integer)
    function react!(integrator)
        integrator.u[2] += event
        integrator.u[i] += event
        DiffEqCallbacks.set_proposed_dt!(integrator, 1e-12)

        return nothing
    end

    return react!
end

function reaction_generator_Ndtot(event::Integer)
    function react!(integrator)
        integrator.u[2] += event
        DiffEqCallbacks.set_proposed_dt!(integrator, 1e-12)

        return nothing
    end

    return react!
end

function excess_Nd(u, t, integrator)
    lambda = u[1]

    dt = get_proposed_dt(integrator)
    # println("Test for excess Nd, dt = $dt")
    dlambda = get_du(integrator)[1] * dt
    next_lambda = lambda + dlambda
    # println("Excess Nd next lambda = $next_lambda")
    l = lambda_to_l_discrete(next_lambda, integrator.p)
    # l = lambda_to_l(next_lambda, integrator.p)
    if any(i -> i > l, u.u[3:end])
        return true
    end

    return false
end

function unbind_excess_Nd!(integrator)
    lambda = integrator.u[1]

    dlambda = get_du(integrator)[1] * get_proposed_dt(integrator)
    next_lambda = lambda + dlambda
    # println("Unbind excess Nd next lambda = $next_lambda")
    l = trunc(lambda_to_l(next_lambda, integrator.p))
    Ndtot_diff = 0
    for i in 3:length(integrator.u.u)
        if integrator.u.u[i] > l
            diff = l - integrator.u.u[i]
            integrator.u.u[i] += diff
            Ndtot_diff += diff
        end
    end
    integrator.u.u[2] += Ndtot_diff
    set_proposed_dt!(integrator, 1e-12)

    return nothing
end

function excess_Ndtot(u, t, integrator)
    lambda = u[1]
    Nd = u[2]

    dt = get_proposed_dt(integrator)
    dlambda = get_du(integrator)[1] * dt
    next_lambda = lambda + dlambda
    l = lambda_to_l_discrete(next_lambda, integrator.p)
    # l = lambda_to_l(next_lambda, integrator.p)
    if Nd > l
        return true
    end

    return false
end

function unbind_excess_Ndtot!(integrator)
    lambda = integrator.u[1]
    Nd = integrator.u[2]

    dlambda = get_du(integrator)[1] * get_proposed_dt(integrator)
    next_lambda = lambda + dlambda
    l = lambda_to_l_discrete(next_lambda, integrator.p)
    # l = lambda_to_l(next_lambda, integrator.p)
    if Nd > l
        println("Unbinding excess Nd")
        diff = l - Nd
        integrator.u.u[2] += diff
    end
    set_proposed_dt!(integrator, 1e-12)

    return nothing
end

function occupancy_break!(integrator)
    lambda = integrator.u[1]
    Nd = integrator.u[2]

    dlambda = get_du(integrator)[1] * get_proposed_dt(integrator)
    next_lambda = lambda + dlambda
    l = lambda_to_l_discrete(next_lambda, integrator.p)
    if Nd > l
        println("Breaking due to max occupancy")
    end
    set_proposed_dt!(integrator, 1e-12)
        integrator.u.u[1] += 1

    return nothing
end

function noninteger_Nd(u, t, integrator)
    if any(i -> !isinteger(i), u.u[2:end])
        return true
    end

    return false
end

function round_noninteger_Nd!(integrator)
    for i in 2:length(integrator.u.u)
        if !isinteger(integrator.u.u[i])
            Nd = integrator.u.u[i]
            println("Applying round crosslinker $Nd")
            integrator.u.u[i] = round(integrator.u.u[i])
        end
    end

    return nothing
end

function create_jumps(overlaps)
    jumps = []
    for i in 1:overlaps
        index = i + 2
        binding_rate = binding_rate_generator(index)
        bind! = reaction_generator(index, 1)
        binding_jump = VariableRateJump(binding_rate, bind!)
        push!(jumps, binding_jump)
        unbinding_rate = unbinding_rate_generator(index)
        unbind! = reaction_generator(index, -1)
        unbinding_jump = VariableRateJump(unbinding_rate, unbind!)
        push!(jumps, unbinding_jump)
    end

    return jumps
end

function create_jumps_Ndtot(overlaps)
    jumps = []
    binding_rate = binding_rate_generator(2)
    bind! = reaction_generator_Ndtot(1)
    binding_jump = VariableRateJump(binding_rate, bind!)
    push!(jumps, binding_jump)
    unbinding_rate = unbinding_rate_generator(2)
    unbind! = reaction_generator_Ndtot(-1)
    unbinding_jump = VariableRateJump(unbinding_rate, unbind!)
    push!(jumps, unbinding_jump)
end

function create_callbacks()
    excess_Nd_cb = DiscreteCallback(excess_Nd, unbind_excess_Nd!)
    noninteger_Nd_cb = DiscreteCallback(noninteger_Nd, round_noninteger_Nd!)

    return CallbackSet(excess_Nd_cb, noninteger_Nd_cb)
end

function create_callbacks_Ndtot()
    excess_Nd_cb = DiscreteCallback(excess_Ndtot, unbind_excess_Ndtot!)
    noninteger_Nd_cb = DiscreteCallback(noninteger_Nd, round_noninteger_Nd!)

    return CallbackSet(excess_Nd_cb, noninteger_Nd_cb)
end

function create_callbacks_Ndtot_noise()
    excess_Nd_cb = DiscreteCallback(excess_Ndtot, occupancy_break!)
    noninteger_Nd_cb = DiscreteCallback(noninteger_Nd, round_noninteger_Nd!)
    zero_overlap_cb = ContinuousCallback(zero_overlap, terminate!)

    return CallbackSet(excess_Nd_cb, noninteger_Nd_cb, zero_overlap_cb)
end

function save_and_write_continuous_Nd(df, filebase, params, ifields)
    filename = savename(filebase, params, suffix=".dat", ignored_fields=ifields)
    CSV.write(filename, df, delim=" ")

    return nothing
end

"""
$(TYPEDSIGNATURES)

Solve with crosslinker-binding quasi-equilibrium and write values to file.

Use double exponential friction expression.
"""
function solve_and_write_ring_cX(u0, tspan, params, ifields, filebase)
    prob = ODEProblem(equation_of_motion_ring_cX!, u0, tspan, params)
    sol = solve(prob, Rosenbrock23())
    lambda = [u[1] for u in sol.u]

    df = calc_cX_quantities(lambda, sol.t, params)
    save_and_write_continuous_Nd(df, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_cX_noise(u0, tspan, params, ifields, filebase)
    prob = SDEProblem(equation_of_motion_ring_cX!, noise_ring_cX!, u0, tspan, params)
    cb = ContinuousCallback(zero_overlap, terminate!)
    sol = solve(prob, SKenCarp(), callback=cb)
    lambda = [u[1] for u in sol.u]

    df = calc_cX_quantities(lambda, sol.t, params)
    save_and_write_continuous_Nd(df, filebase, params, ifields)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Solve with crosslinker-diffusion quasi-equilibrium and write values to file.

Use double exponential friction expression.
"""
function solve_and_write_ring_Nd_exp(u0, tspan, params, ifields, filebase)
    prob = ODEProblem(equation_of_motion_ring_Nd_exp!, u0, tspan, params)
    sol = solve(prob, Rosenbrock23())
    lambda = [u[1] for u in sol.u]
    Ndtot = [u[2] for u in sol.u]

    df = calc_Nd_quantities(lambda, Ndtot, sol.t, params)
    save_and_write_continuous_Nd(df, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_exp_noise(u0, tspan, params, ifields, filebase)
    prob = SDEProblem(equation_of_motion_ring_Nd_exp!, noise_ring_Nd_exp!, u0, tspan, params)
    cb = ContinuousCallback(zero_overlap, terminate!)
    sol = solve(prob, SKenCarp(), callback=cb)
    lambda = [u[1] for u in sol.u]
    Ndtot = [u[2] for u in sol.u]

    df = calc_Nd_quantities(lambda, Ndtot, sol.t, params)
    save_and_write_continuous_Nd(df, filebase, params, ifields)

    return nothing
end

function save_and_write_discrete_Nd(dfs, filebase, params, ifields)
    (df_means, df_vars) = meanvar_dfs(dfs, params.interval)
    filename = savename(filebase, params, suffix="_means.dat", ignored_fields=ifields)
    CSV.write(filename, df_means, delim=" ")
    filename = savename(filebase, params, suffix="_vars.dat", ignored_fields=ifields)
    CSV.write(filename, df_vars, delim=" ")

    return nothing
end

function solve_and_write_ring_Nd_exact_Ndtot_base(oprob, trajs, params, ifields, filebase)
    overlaps = 2params.Nf - params.Nsca
    jumps = create_jumps_Ndtot(overlaps)
    jprob = JumpProblem(oprob, Direct(), jumps...)
    eprob = EnsembleProblem(jprob)
    cb = create_callbacks_Ndtot()
    solarray = solve(eprob, Rosenbrock23(), EnsembleThreads(), callback=cb, trajectories=trajs)

    dfs = []
    for (i, sol) in enumerate(solarray)
        lambda = [uti[1] for uti in sol.u]
        Ndtot = [uti[2] * overlaps for uti in sol.u]
        Ndi = [uti[2] for uti in sol.u]

        df = calc_discrete_Nd_quantities(lambda, Ndtot, Ndi, sol.t, params)
        push!(dfs, df)

        filename = savename(filebase, params, suffix="_$i.dat", ignored_fields=ifields)
        CSV.write(filename, df, delim=" ")
    end
    save_and_write_discrete_Nd(dfs, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_exact_Ndtot_noise_base(oprob, trajs, params, ifields, filebase)
    overlaps = 2params.Nf - params.Nsca
    jumps = create_jumps_Ndtot(overlaps)
    jprob = JumpProblem(oprob, Direct(), jumps...)
    eprob = EnsembleProblem(jprob)
    cb = create_callbacks_Ndtot_noise()
    solarray = solve(eprob, SOSRA(), EnsembleThreads(), callback=cb, trajectories=trajs)

    dfs = []
    for (i, sol) in enumerate(solarray)
        lambda = [uti[1] for uti in sol.u]
        Ndtot = [uti[2] * overlaps for uti in sol.u]
        Ndi = [uti[2] for uti in sol.u]

        df = calc_discrete_Nd_quantities(lambda, Ndtot, Ndi, sol.t, params)
        push!(dfs, df)

        filename = savename(filebase, params, suffix="_$i.dat", ignored_fields=ifields)
        CSV.write(filename, df, delim=" ")
    end
    #save_and_write_discrete_Nd(dfs, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_exact_Nds_base(oprob, trajs, params, ifields, filebase)
    overlaps = 2params.Nf - params.Nsca
    jumps = create_jumps(overlaps)
    jprob = JumpProblem(oprob, Direct(), jumps...)
    eprob = EnsembleProblem(jprob)
    cb = create_callbacks()
    solarray = solve(eprob, Rosenbrock23(), EnsembleThreads(), callback=cb, trajectories=trajs)

    dfs = []
    for (i, sol) in enumerate(solarray)
        lambda = [uti[1] for uti in sol.u]
        Ndtot = [uti[2] for uti in sol.u]
        Nds = [uti.u[3:end] for uti in sol.u]

        df = calc_discrete_Nd_quantities(lambda, Ndtot, Nds, sol.t, params)
        push!(dfs, df)

        filename = savename(filebase, params, suffix="_$i.dat", ignored_fields=ifields)
        CSV.write(filename, df, delim=" ")
    end
    save_and_write_discrete_Nd(dfs, filebase, params, ifields)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Sample trajectories with discrete Nd and write values to file.

Use exact expression for friction coefficient but with continuous l. The total number of
crosslinkers varies by the number of overlaps in order to keep the number per overlap an
integer.
"""
function solve_and_write_ring_Nd_exact_Ndtot(u0, tspan, trajs, params, ifields, filebase)
    oprob = ODEProblem(equation_of_motion_ring_Nd_exact_Ndtot!, u0, tspan, params)
    solve_and_write_ring_Nd_exact_Ndtot_base(oprob, trajs, params, ifields, filebase)

    return nothing
end

function solve_and_write_ring_Nd_exact_Ndtot_noise(u0, tspan, trajs, params, ifields, filebase)
    prob = SDEProblem(equation_of_motion_ring_Nd_exact_Ndtot!, noise_ring_Nd_exact_Ndtot!,
                      u0, tspan, params)
    solve_and_write_ring_Nd_exact_Ndtot_noise_base(prob, trajs, params, ifields, filebase)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Sample trajectories with discrete Nd and write values to file.

Use exact expression for friction coefficient but with continuous l. The number of
crosslinkers in each overlap is tracked, and the friction coefficient is calculated by
averaging over the values calculated for each overlap.
"""
function solve_and_write_ring_Nd_exact_Nds(u0, tspan, trajs, params, ifields, filebase)
    oprob = ODEProblem(equation_of_motion_ring_Nd_exact_Nds!, u0, tspan, params)
    solve_and_write_ring_Nd_exact_Nds_base(oprob, trajs, params, ifields, filebase)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Sample trajectories with discrete Nd but no sliding and write values to file.

Use exact expression for friction coefficient but with continuous l. The total number of
crosslinkers varies by the number of overlaps in order to keep the number per overlap an
integer.
"""
function solve_and_write_ring_Nd_exact_Ndtot_static(u0, tspan, trajs, params, ifields, filebase)
    oprob = ODEProblem(equation_of_motion_ring_Nd_exact_Ndtot_static!, u0, tspan, params)
    solve_and_write_ring_Nd_exact_Ndtot_base(oprob, trajs, params, ifields, filebase)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Sample trajectories with discrete Nd but a constant friction coefficient.

The total number of crosslinkers varies by the number of overlaps in order to keep the
number per overlap an integer.
"""
function solve_and_write_ring_Nd_constant_discrete_Ndtot(u0, tspan, trajs, params, ifields,
        filebase)
    oprob = ODEProblem(equation_of_motion_ring_Nd_constant_discrete_Ndtot!, u0, tspan,
                       params)
    solve_and_write_ring_Nd_exact_Ndtot_base(oprob, trajs, params, ifields, filebase)

    return nothing
end
