"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker-binding quasi-equilibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_cX!(du, u, p, t)
    zeta = friction_coefficient_ring_cX(u[1], p)
    forcetot = bending_force(u[1], p) + condensation_force(p)

    du[1] = -forcetot / (zeta * p.deltas * (2p.Nf - p.Nsca))

    return nothing
end

"""
$(TYPEDSIGNATURES)

"""
function equation_of_motion_ring_Nd_update(du, u, p, zeta)
    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(u[1], p) + entropic_force(u[1], u[2], p)
    ltot = (1 + p.deltas / p.deltad * u[1]) * overlaps
    du[1] = -forcetot / (zeta * p.deltas * overlaps)
    du[2] = p.cX * p.k01 * p.r12 * ltot - (p.cX * p.k01 * p.r12 + p.r21 * p.r10) * u[2]

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker diffusion quasi-equlibrium.

This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_ring_Nd!(du, u, p, t)
    zeta = friction_coefficient_ring_Nd(u[1], u[2], p)
    equation_of_motion_ring_Nd_update(du, u, p, zeta)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

Does not assume that the spring constant of the crosslinker is larger relative to thermal
fluctations. This is compatible with the DifferentialEquations package.
"""
function equation_of_motion_single_exp_ring_Nd!(du, u, p, t)
    zeta = friction_coefficient_single_exp_ring_Nd(u[1], u[2], p)
    equation_of_motion_ring_Nd_update(du, u, p, zeta)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This requires the number of crosslinkers per overlap to be tracked, and the
number of bound crosslinkers per overlap to be updated with a jump process separately. This
is compatible with the DifferentialEquations package.
"""
function equation_of_motion_continuous_l_ring_Nd!(du, u, p, t)
    du .= 0
    l = trunc(1 + p.deltas / p.deltad * u[1])

    # This should never happen
    if any(isnan, u)
        println(u)
        println("NaN in solution")
        du .= NaN
        return nothing
    end

    # Callbacks are not called before every call of this, so I have to check the domain here
    if any(i -> (i < 0), u)
        println("Setting du to inf from negative value, $u")
        du .= Inf
        return nothing
    end

    # Same issue as above
    if any(i -> (i > l), u[3:end])
        Nds = u[3:end]
        println("Setting du to inf from excess Nd, lambda = $u[1], l = $l, Nds = $Nds")
        du .= Inf
        return nothing
    end

    zeta = friction_coefficient_continuous_l_ring_Nd(u[1], u[3:end], p)
    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(u[1], p) + entropic_force(u[1], u[2], p)
    ltot = (1 + p.deltas / p.deltad * u[1]) * overlaps

    du[1] = -forcetot / (zeta * p.deltas * overlaps)
    lambda = u[1] + du[1]
    l = trunc(1 + p.deltas / p.deltad * lambda)
    for i in 2:length(du)
        du[i] = 0
    end

#    println("Total force = $forcetot")
#    println("zeta = $zeta")
#    dlambda = -forcetot / (zeta * p.deltas * overlaps)
#    println("Delta lambda = $dlambda")
#    println("u = $u")
#    println("du = $du")
    if any(isnan, du)
#        println("Hmm2")
        du .= Inf
        return nothing
    end
    #if any(i -> (i < 0), u + du)
    #    println("Hmm3")
    #    du .= Inf
    #    return nothing
    #end

    # Same issue as above
    #if any(i -> (i > l), (u + du)[3:end])
    #    println("Hmm4")
    #    du .= Inf
    #    return nothing
    #end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Generate crosslinker binding rate function.
"""
function binding_rate_generator(i::Integer)
    function binding_rate(u, p, t)
        l = 1 + p.deltas / p.deltad * u[1]
        if trunc(l) == u[i]
            return 0.0
        else
            return p.cX * p.k01 * p.r12 * (l - u[i])
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
        l = 1 + p.deltas / p.deltad * u[1]
        if u[i] == 1.0
            return 0.0
        else
            return p.r21 * p.r10 * u[i]
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
        DiffEqCallbacks.set_proposed_dt!(integrator, 1e-9)

        return nothing
    end

    return react!
end

function excess_Nd(u, t, integrator)
    lambda = u[1]
    next_lambda = u[1] + get_du(integrator)[1]*get_proposed_dt(integrator)
    l = trunc(1 + integrator.p.deltas / integrator.p.deltad * next_lambda)
    overlaps = u.u[3:end]
#    println("Testing if excess, lambda = $lambda, lambda' = $next_lambda, l = $l, Nds = $overlaps")
    if any(i -> i > l, u.u[3:end])
        return true
    end

    return false
end

function unbind_excess_Nd!(integrator)
    next_lambda = integrator.u[1] + get_du(integrator)[1]*get_proposed_dt(integrator)
    l = trunc(1 + integrator.p.deltas / integrator.p.deltad * next_lambda)
    Ndtot_diff = 0
    for i in 3:length(integrator.u.u)
        if integrator.u.u[i] > l
            diff = l - integrator.u.u[i]
            println("Applying unbind excess")
            integrator.u.u[i] += diff
            Ndtot_diff += diff
        end
    end
    integrator.u.u[2] += Ndtot_diff
    set_proposed_dt!(integrator, 1e-9)

    return nothing
end

function noninteger_Nd(u, t, integrator)
    Nds = u.u[2:end]
#    println("Testing if non-integer, $Nds")
#    println(get_proposed_dt(integrator))
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

function create_callbacks()
    excess_Nd_cb = DiscreteCallback(excess_Nd, unbind_excess_Nd!)
    noninteger_Nd_cb = DiscreteCallback(noninteger_Nd, round_noninteger_Nd!)

    return CallbackSet(excess_Nd_cb, noninteger_Nd_cb)
end

"""
$(TYPEDSIGNATURES)

Solve with crosslinker-binding quasi-equilibrium and write values to file.

Use double exponential friction expression.
"""
function solve_and_write_cX(u0, tspan, params, ifields, filebase)
    prob = ODEProblem(equation_of_motion_ring_cX!, u0, tspan, params)
    sol = solve(prob, Rosenbrock23())
    lambda = [u[1] for u in sol.u]

    df = calc_cX_quantities(lambda, sol.t, params)

    filename = savename(filebase, params, suffix=".dat", ignored_fields=ifields)
    CSV.write(filename, df, delim=" ")
end

"""
$(TYPEDSIGNATURES)

Solve with crosslinker-diffusion quasi-equilibrium and write values to file.

Use double exponential friction expression.
"""
function solve_and_write_double_exp(u0, tspan, params, ifields, filebase)
    prob = ODEProblem(equation_of_motion_ring_Nd!, u0, tspan, params)
    sol = solve(prob, Rosenbrock23())
    lambda = [u[1] for u in sol.u]
    Ndtot = [u[2] for u in sol.u]

    df = calc_Nd_quantities(lambda, Ndtot, sol.t, params)

    filename = savename(filebase, params, suffix=".dat", ignored_fields=ifields)
    CSV.write(filename, df, delim=" ")
end

"""
$(TYPEDSIGNATURES)

Solve with crosslinker-diffusion quasi-equilibrium and write values to file.

Use single exponential friction expression.
"""
function solve_and_write_single_exp(u0, tspan, params, ifields, filebase)
    prob = ODEProblem(equation_of_motion_single_exp_ring_Nd!, u0, tspan, params)
    sol = solve(prob, Rosenbrock23())
    lambda = [u[1] for u in sol.u]
    Ndtot = [u[2] for u in sol.u]

    df = calc_Nd_quantities(lambda, Ndtot, sol.t, params)

    filename = savename(filebase, params, suffix=".dat", ignored_fields=ifields)
    CSV.write(filename, df, delim=" ")
end

"""
$(TYPEDSIGNATURES)

Sample trajectories with discrete Nd and write values to file.

Use exact expression for friction coefficient but with continuous l.
"""
function solve_and_write_continuous_l(u0, tspan, trajs, params, ifields, filebase)
    overlaps = 2params.Nf - params.Nsca
    oprob = ODEProblem(equation_of_motion_continuous_l_ring_Nd!, u0, tspan, params)
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

    (df_means, df_vars) = meanvar_dfs(dfs)
    filename = savename(filebase, params, suffix="_means.dat", ignored_fields=ifields)
    CSV.write(filename, df_means, delim=" ")
    filename = savename(filebase, params, suffix="_vars.dat", ignored_fields=ifields)
    CSV.write(filename, df_vars, delim=" ")
end
