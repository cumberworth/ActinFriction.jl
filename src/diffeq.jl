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

function equation_of_motion_ring_cX_noise!(du, u, p, t)
    lambda = u[1]

    overlaps = 2p.Nf - p.Nsca
    zeta = friction_coefficient_cX(lambda, p)
    C = friction_coefficient_cX_C(p)
    forcetot = bending_force(lambda, p) + condensation_force(p)
    sdrift = -kb * p.T * log(C) / (p.deltas * p.deltad * overlaps^2 * zeta)

    du[1] = forcetot / (zeta * p.deltas * overlaps) + sdrift

    return nothing
end

function noise_ring_cX!(du, u, p, t)
    lambda = u[1]

    overlaps = 2p.Nf - p.Nsca
    zeta = friction_coefficient_cX(lambda, p)
    du[1] = sqrt(2 * kb * p.T / zeta) / p.deltas / overlaps
end

function equation_of_motion_ring_Nd_contin_update!(du, u, p, drift)
    lambda = u[1]
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(lambda, p) + entropic_force(lambda, Ndtot, p)
    ltot = lambda_to_l(lambda, p) * overlaps
    A = 2 * p.r10 * p.r12 / ((p.r01 + p.r10)^2 + p.r10 * p.r12)
    B = p.r01 * ltot + p.r21 * Ndtot - p.r01 * Ndtot
    C = 2 * p.r21 * Ndtot

    du[1] = drift
    du[2] = B * A - C

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker diffusion quasi-equlibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_contin_exp!(du, u, p, t)
    lambda = u[1]
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exp(Nd, p)
    equation_of_motion_ring_Nd_contin_update!(du, u, p, zeta)

    return nothing
end

function noise_ring_Nd_exp!(du, u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exp(Nd, p)
    du[1] = sqrt(2 * kb * p.T / zeta) / p.deltas / overlaps

    return nothing
end

function equation_of_motion_ring_Nd_discrete_base!(zeta, Ndtot, du, u, p, t)
    lambda = u[1]
    du .= 0

    forcetot = bending_force(lambda, p) + entropic_force(lambda, Ndtot, p)
    overlaps = 2p.Nf - p.Nsca

    du[1] = forcetot / (zeta * p.deltas * overlaps)

    return nothing
end

function equation_of_motion_ring_Nd_discrete_static_base!(zeta, Ndtot, du, u, p, t)
    du .= 0

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

"""
function equation_of_motion_ring_Nd_discrete_exp!(du, u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exp(Nd, p)
    equation_of_motion_ring_Nd_discrete_base!(zeta, Ndtot, du, u, p, t)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_discrete_exact!(du, u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exact(Nd, p)
    equation_of_motion_ring_Nd_discrete_base!(zeta, Ndtot, du, u, p, t)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_discrete_constant!(du, u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    equation_of_motion_ring_Nd_discrete_base!(p.zeta, Ndtot, du, u, p, t)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd_discrete_exact_static!(du, u, p, t)
    Nd = u[2]
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exact(Nd, p)
    equation_of_motion_ring_Nd_discrete_static_base!(zeta, Ndtot, du, u, p, t)

    return nothing
end

function noise_ring_Nd_exact!(du, u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    zeta = friction_coefficient_Nd_exact(Ndtot, p)
    du[1] = sqrt(2 * kb * p.T / zeta) / p.deltas / overlaps

    return nothing
end

function binding_rate(Nd, l, p)
    A = 2 * p.r10 * p.r12 / ((p.r01 + p.r10)^2 + p.r10 * p.r12)
    B = p.r01 * l + p.r21 * Nd - p.r01 * Nd

    return A * B 
end

function Nd_binding_rate(u, p, t)
    lambda = u[1]
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    l = lambda_to_l(lambda, p)
    l_discrete = lambda_to_l_discrete(u[1], p)
    if l_discrete == Nd
        return 0.0
    else
        return binding_rate(Nd, l, p)
    end
end

function Ndtot_binding_rate(u, p, t)
    lambda = u[1]
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    lambda_tot = lambda * overlaps
    ltot = lambda_to_l(lambda_tot, p)
    ltot_discrete = lambda_to_l_discrete(u[1], p)
    if ltot_discrete == Ndtot
        return 0.0
    else
        return binding_rate(Ndtot, ltot, p)
    end
end

function Nd_unbinding_rate(u, p, t)
    Ndtot = u[2]

    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    if Nd == 1.0
        return 0.0
    else
        return 2 * p.r21 * Nd
    end
end

function Ndtot_unbinding_rate(u, p, t)
    Ndtot = u[2]

    if Ndtot == 1.0
        return 0.0
    else
        return 2 * p.r21 * Ndtot
    end
end

function reaction_generator(event::Integer)
    function react!(integrator)
        integrator.u[2] += event
        DiffEqCallbacks.set_proposed_dt!(integrator, 1e-12)

        return nothing
    end

    return react!
end

# finish converting to u[2] is Ndtot
function excess_Nd(u, t, integrator)
    lambda = u[1]
    Ndtot = u[2]

    dt = get_proposed_dt(integrator)
    dlambda = get_du(integrator)[1] * dt
    next_lambda = lambda + dlambda
    overlaps = 2integrator.p.Nf - integrator.p.Nsca
    Nd = Ndtot / overlaps
    l = lambda_to_l_discrete(next_lambda, integrator.p)
    if Nd > l
        return true
    end

    return false
end

function unbind_excess_Nd!(integrator)
    lambda = integrator.u[1]
    Ndtot = integrator.u[2]

    overlaps = 2integrator.p.Nf - integrator.p.Nsca
    Nd = Ndtot / overlaps
    dlambda = get_du(integrator)[1] * get_proposed_dt(integrator)
    next_lambda = lambda + dlambda
    l = lambda_to_l_discrete(next_lambda, integrator.p)
    println("Unbinding excess Nd")
    diff = l - Nd
    integrator.u.u[2] += diff
    set_proposed_dt!(integrator, 1e-12)

    return nothing
end


function unbind_excess_Ndtot!(integrator)
    lambda = integrator.u[1]
    Ndtot = integrator.u[2]

    overlaps = 2integrator.p.Nf - integrator.p.Nsca
    dlambda = get_du(integrator)[1] * get_proposed_dt(integrator)
    next_lambda = lambda + dlambda
    l = lambda_to_l_discrete(next_lambda, integrator.p)
    ltot = l*overlaps
    println("Unbinding excess Ndtot")
    diff = ltot - Ndtot
    integrator.u.u[2] += diff
    set_proposed_dt!(integrator, 1e-12)

    return nothing
end

# Instead of unbinding excess, the overlap bounces off a barrier
function occupancy_break!(integrator)
    println("Breaking due to max occupancy")
    set_proposed_dt!(integrator, 1e-12)

    # This is an arbitrary selection
    integrator.u.u[1] += 1

    return nothing
end

function noninteger_Ndtot(u, t, integrator)
    return !isinteger(u.u[2])
end

function round_noninteger_Ndtot!(integrator)
    Ndtot = integrator.u.u[2]
    println("Applying round crosslinker $Ndtot")
    integrator.u.u[2] = round(integrator.u.u[2])

    return nothing
end

function create_jumps_Nd(overlaps)
    jumps = []
    bind! = reaction_generator(overlaps)
    binding_jump = VariableRateJump(Nd_binding_rate, bind!)
    push!(jumps, binding_jump)
    unbind! = reaction_generator(-overlaps)
    unbinding_jump = VariableRateJump(Nd_unbinding_rate, unbind!)
    push!(jumps, unbinding_jump)

    return jumps
end

function create_jumps_Ndtot()
    jumps = []
    bind! = reaction_generator(1)
    binding_jump = VariableRateJump(Ndtot_binding_rate, bind!)
    push!(jumps, binding_jump)
    unbind! = reaction_generator(-1)
    unbinding_jump = VariableRateJump(Ndtot_unbinding_rate, unbind!)
    push!(jumps, unbinding_jump)

    return jumps
end

function create_callbacks_Nd()
    excess_Nd_cb = DiscreteCallback(excess_Nd, unbind_excess_Nd!)
    noninteger_Ndtot_cb = DiscreteCallback(noninteger_Ndtot, round_noninteger_Ndtot!)

    return CallbackSet(excess_Nd_cb, noninteger_Ndtot_cb)
end

function create_callbacks_Ndtot()
    excess_Nd_cb = DiscreteCallback(excess_Nd, unbind_excess_Ndtot!)
    noninteger_Ndtot_cb = DiscreteCallback(noninteger_Ndtot, round_noninteger_Ndtot!)

    return CallbackSet(excess_Nd_cb, noninteger_Ndtot_cb)
end

function create_callbacks_Nd_noise()
    excess_Nd_cb = DiscreteCallback(excess_Nd, occupancy_break!)
    noninteger_Ndtot_cb = DiscreteCallback(noninteger_Ndtot, round_noninteger_Ndtot!)
    zero_overlap_cb = ContinuousCallback(zero_overlap, terminate!)

    return CallbackSet(excess_Nd_cb, noninteger_Ndtot_cb, zero_overlap_cb)
end

function create_callbacks_Ndtot_noise()
    excess_Nd_cb = DiscreteCallback(excess_Nd, occupancy_break!)
    noninteger_Ndtot_cb = DiscreteCallback(noninteger_Ndtot, round_noninteger_Ndtot!)
    zero_overlap_cb = ContinuousCallback(zero_overlap, terminate!)

    return CallbackSet(excess_Nd_cb, noninteger_Ndtot_cb, zero_overlap_cb)
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
    prob = SDEProblem(equation_of_motion_ring_cX_noise!, noise_ring_cX!, u0, tspan, params)
    cb = ContinuousCallback(zero_overlap, terminate!)
    sol = solve(prob, SOSRI(), callback=cb)
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
function solve_and_write_ring_Nd_contin_exp(u0, tspan, params, ifields, filebase)
    prob = ODEProblem(equation_of_motion_ring_Nd_contin_exp!, u0, tspan, params)
    sol = solve(prob, Rosenbrock23())
    lambda = [u[1] for u in sol.u]
    Ndtot = [u[2] for u in sol.u]

    df = calc_Nd_quantities(lambda, Ndtot, sol.t, params)
    save_and_write_continuous_Nd(df, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_contin_exp_noise(u0, tspan, params, ifields, filebase)
    prob = SDEProblem(equation_of_motion_ring_Nd_contin_exp!, noise_ring_Nd_exp!, u0, tspan, params)
    cb = ContinuousCallback(zero_overlap, terminate!)
    sol = solve(prob, SOSRI(), callback=cb)
    lambda = [u[1] for u in sol.u]
    Ndtot = [u[2] for u in sol.u]

    df = calc_Nd_quantities(lambda, Ndtot, sol.t, params)
    save_and_write_continuous_Nd(df, filebase, params, ifields)

    return nothing
end

function save_and_write_discrete_Nd(solarray, filebase, params, ifields)
    dfs = []
    for (i, sol) in enumerate(solarray)
        lambda = [uti[1] for uti in sol.u]
        Ndtot = [uti[2] for uti in sol.u]

        df = calc_discrete_Nd_quantities(lambda, Ndtot, sol.t, params)
        push!(dfs, df)

        filename = savename(filebase, params, suffix="_$i.dat", ignored_fields=ifields)
        CSV.write(filename, df, delim=" ")
    end

    #(df_means, df_vars) = meanvar_dfs(dfs, params.interval)
    #filename = savename(filebase, params, suffix="_means.dat", ignored_fields=ifields)
    #CSV.write(filename, df_means, delim=" ")
    #filename = savename(filebase, params, suffix="_vars.dat", ignored_fields=ifields)
    #CSV.write(filename, df_vars, delim=" ")

    return nothing
end

function solve_and_write_ring_Nd_discrete_exp_base(oprob, trajs, params, ifields, filebase)
    jumps = create_jumps_Ndtot()
    jprob = JumpProblem(oprob, Direct(), jumps...)
    eprob = EnsembleProblem(jprob)
    cb = create_callbacks_Ndtot()
    solarray = solve(eprob, Rosenbrock23(), EnsembleThreads(), callback=cb, trajectories=trajs)
    save_and_write_discrete_Nd(solarray, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_discrete_exp_noise_base(sprob, trajs, params, ifields, filebase)
    jumps = create_jumps_Ndtot()
    jprob = JumpProblem(sprob, Direct(), jumps...)
    eprob = EnsembleProblem(jprob)
    cb = create_callbacks_Ndtot_noise()
    solarray = solve(eprob, SOSRI(), EnsembleThreads(), callback=cb, trajectories=trajs,
                     abstol=1, reltol=1)
    save_and_write_discrete_Nd(solarray, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_discrete_exact_base(oprob, trajs, params, ifields, filebase)
    overlaps = 2params.Nf - params.Nsca
    jumps = create_jumps_Nd(overlaps)
    jprob = JumpProblem(oprob, Direct(), jumps...)
    eprob = EnsembleProblem(jprob)
    cb = create_callbacks_Ndtot()
    solarray = solve(eprob, Rosenbrock23(), EnsembleThreads(), callback=cb, trajectories=trajs)
    save_and_write_discrete_Nd(solarray, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_discrete_exact_noise_base(sprob, trajs, params, ifields, filebase)
    overlaps = 2params.Nf - params.Nsca
    jumps = create_jumps_Ndtot()
    jprob = JumpProblem(sprob, Direct(), jumps...)
    eprob = EnsembleProblem(jprob)
    cb = create_callbacks_Ndtot_noise()
    solarray = solve(eprob, SOSRI(), EnsembleThreads(), callback=cb, trajectories=trajs,
                     abstol=1, reltol=1)
    save_and_write_discrete_Nd(solarray, filebase, params, ifields)

    return nothing
end

function solve_and_write_ring_Nd_discrete_exp(u0, tspan, trajs, params, ifields, filebase)
    oprob = ODEProblem(equation_of_motion_ring_Nd_discrete_exp!, u0, tspan, params)
    solve_and_write_ring_Nd_discrete_exp_base(oprob, trajs, params, ifields, filebase)

    return nothing
end

function solve_and_write_ring_Nd_discrete_exp_noise(u0, tspan, trajs, params, ifields, filebase)
    sprob = SDEProblem(equation_of_motion_ring_Nd_discrete_exp!, noise_ring_Nd_exp!, u0,
                       tspan, params)
    solve_and_write_ring_Nd_discrete_exp_noise_base(sprob, trajs, params, ifields, filebase)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Sample trajectories with discrete Nd and write values to file.

Use exact expression for friction coefficient but with continuous l. The total number of
crosslinkers varies by the number of overlaps in order to keep the number per overlap an
integer.
"""
function solve_and_write_ring_Nd_discrete_exact(u0, tspan, trajs, params, ifields, filebase)
    oprob = ODEProblem(equation_of_motion_ring_Nd_discrete_exact!, u0, tspan, params)
    solve_and_write_ring_Nd_discrete_exact_base(oprob, trajs, params, ifields, filebase)

    return nothing
end

function solve_and_write_ring_Nd_discrete_exact_noise(u0, tspan, trajs, params, ifields, filebase)
    sprob = SDEProblem(equation_of_motion_ring_Nd_discrete_exact!, noise_ring_Nd_exact!,
                      u0, tspan, params)
    solve_and_write_ring_Nd_discrete_exact_noise_base(sprob, trajs, params, ifields, filebase)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Sample trajectories with discrete Nd but no sliding and write values to file.

Use exact expression for friction coefficient but with continuous l. The total number of
crosslinkers varies by the number of overlaps in order to keep the number per overlap an
integer.
"""
function solve_and_write_ring_Nd_discrete_exact_static(u0, tspan, trajs, params, ifields, filebase)
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
function solve_and_write_ring_Nd_discrete_constant(u0, tspan, trajs, params, ifields,
        filebase)
    oprob = ODEProblem(equation_of_motion_ring_Nd_discrete_constant!, u0, tspan,
                       params)
    solve_and_write_ring_Nd_discrete_exact_base(oprob, trajs, params, ifields, filebase)

    return nothing
end
