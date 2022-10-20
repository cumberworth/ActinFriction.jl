"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker-binding quasi-equilibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_cX!(du, u, p, t)
    zeta = friction_coefficient_cX(u[1], p)
    forcetot = bending_force(u[1], p) + condensation_force(p)

    du[1] = forcetot / (zeta * p.deltas * (2p.Nf - p.Nsca))

    return nothing
end

function equation_of_motion_ring_Nd_update!(du, u, p, zeta)
    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(u[1], p) + entropic_force(u[1], u[2], p)
    ltot = lambda_to_l(u[1], p) * overlaps
    A = 2 * p.r10 * p.r12 / ((p.r01 + p.r10)^2 + p.r10 * p.r12)
    B = p.r01 * ltot + p.r21 * u[2] - p.r01 * u[2]
    C = 2 * p.r21 * u[2]
    du[1] = forcetot / (zeta * p.deltas * overlaps)
    du[2] =   B * A - C

    return nothing
end

function equation_of_motion_ring_Nd_Ns_update!(du, u, p, zeta)
    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(u[1], p) + entropic_force(u[1], u[2], p)
    ltot = lambda_to_l(u[1], p) * overlaps
    du[1] = forcetot / (zeta * p.deltas * overlaps)
    du[2] = p.r12 * u[3] * (1 - u[3] / 2 / ltot) - 2 * p.r21 * u[2]
    du[3] = p.cX * p.k01 * (2ltot - u[3] - 2u[2]) - p.r10 * u[3]

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker diffusion quasi-equlibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_exp!(du, u, p, t)
    overlaps = 2p.Nf - p.Nsca
    Nd = u[2] / overlaps
    zeta = friction_coefficient_Nd_exp(u[1], Nd, p)
    equation_of_motion_ring_Nd_update!(du, u, p, zeta)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker diffusion quasi-equlibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_exp_Ns!(du, u, p, t)
    overlaps = 2p.Nf - p.Nsca
    Nd = u[2] / overlaps
    zeta = friction_coefficient_Nd_exp(u[1], Nd, p)
    equation_of_motion_ring_Nd_Ns_update!(du, u, p, zeta)

    return nothing
end

function equation_of_motion_ring_Nd_exact_base!(zeta, Ndtot, du, u, p, t)
    du .= 0
    lambda = u[1]
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

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This tracks the number of crosslinkers per overlap, and the
number of bound crosslinkers per overlap to be updated with a jump process separately. This
is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_exact_Nds!(du, u, p, t)
    zeta = friction_coefficient_Nd_exact(u[1], u[3:end], p)
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
    lambda = u[1]
    Nd = u[2]
    overlaps = 2p.Nf - p.Nsca
    Ndtot = Nd * overlaps
    zeta = friction_coefficient_Nd_exact(lambda, [Nd], p)
    equation_of_motion_ring_Nd_exact_base!(zeta, Ndtot, du, u, p, t)

    return nothing
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
    l = tunc(lambda_to_l(next_lambda, integrator.p))
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
        diff = l - Nd
        integrator.u.u[2] += diff
    end
    set_proposed_dt!(integrator, 1e-12)

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

"""
$(TYPEDSIGNATURES)

Solve with crosslinker-diffusion quasi-equilibrium and write values to file.

Use double exponential friction expression.
"""
function solve_and_write_ring_Nd_exp_Ns(u0, tspan, params, ifields, filebase)
    prob = ODEProblem(equation_of_motion_ring_Nd_exp_Ns!, u0, tspan, params)
    sol = solve(prob, Rosenbrock23())
    lambda = [u[1] for u in sol.u]
    Ndtot = [u[2] for u in sol.u]
    Nstot = [u[3] for u in sol.u]

    df = calc_Ns_quantities(lambda, Ndtot, Nstot, sol.t, params)
    save_and_write_continuous_Nd(df, filebase, params, ifields)

    return nothing
end

function save_and_write_discrete_Nd(dfs, filebase, params, ifields)
    (df_means, df_vars) = meanvar_dfs(dfs)
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
        Nds = [[uti.u[2]] for uti in sol.u]

        df = calc_discrete_Nd_quantities(lambda, Ndtot, Nds, sol.t, params)
        push!(dfs, df)

        filename = savename(filebase, params, suffix="_$i.dat", ignored_fields=ifields)
        CSV.write(filename, df, delim=" ")
    end
    save_and_write_discrete_Nd(dfs, filebase, params, ifields)

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
