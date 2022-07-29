"""
$(TYPEDSIGNATURES)

Calculate continuous bionomial coefficient.
"""
function binomialc(n, k)
    return gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))
end

"""
$(TYPEDFIELDS)

Parameters for actin-anillin ring system.
"""
Base.@kwdef struct RingParams
    """Per site rate constant of initial crosslinker binding (s^-1)"""
    k01::Float64
    """Per site rate of initial crosslinker binding (M^-1 s^-1)"""
    r01::Float64
    """Per site rate (constant) of singly bound crosslinker unbinding (s^-1)"""
    r10::Float64
    """Per site rate (constant) of singly bound crosslinker doubly binding (s^-1)"""
    r12::Float64
    """Per site rate (constant) of doubly bound crosslinker unbinding one head (s^-1)"""
    r21::Float64
    """Spacing between binding sites on a filament (m)"""
    deltas::Float64
    """Spacing for doubly bound crosslinkers (m)"""
    deltad::Float64
    """Spring constant of crosslinker (N m^-1)"""
    k::Float64
    """Temperature (K)"""
    T::Float64
    """Number of filaments"""
    Nf::Int
    """Number of scaffold filaments"""
    Nsca::Int
    """ Bending rigidity (N m^2)"""
    EI::Float64
    """Filament length (m)"""
    Lf::Float64
    """Friction coefficient of free actin filament"""
    zeta0::Float64
    """Dissociation constant for single crosslinker binding (M)"""
    KsD::Float64
    """Dissociation constant for double crosslinker binding from solution (M)"""
    KdD::Float64
    """Crosslinker concentration (M)"""
    cX::Float64
end

"""
$(TYPEDSIGNATURES)

Generate file name from a set of parameters.

This is a replacement for the DrWatson function that has better handling of rounding and
trailing zeros.
"""
function savename(prefix, params; digits=2, suffix=nothing)
    fields = fieldnames(typeof(params))
    fieldstrings = [string(f) for f in fields]
    sortedindices = sortperm(fieldstrings)
    filename = [prefix]
    for i in sortedindices
        fieldstring = fieldstrings[i]
        value  = getproperty(params, fields[i])
        if isinteger(value)
            fstring = Printf.Format("%.i")
            vstring = Printf.format(fstring, value)
            push!(filename, "$(fieldstring)=$(vstring)")
        else
            fstring = Printf.Format("%.$(digits)e")
            vstring = Printf.format(fstring, round(value, sigdigits=digits + 1))
            push!(filename, "$(fieldstring)=$(vstring)")
        end
    end
    filename = join(filename, "_")

    if suffix !== nothing
        filename *= suffix
    end

    return filename
end

"""
$(TYPEDSIGNATURES)

Calculate current bending force in a ring.
"""
function bending_force(lambda, p::RingParams)
    F = 8 * pi^3 * p.EI * p.Lf * p.Nf / p.Nsca^3
    G = -(p.deltas^3)
    H = 3 * p.Lf * p.deltas^2
    J = -3 * p.Lf^2 * p.deltas
    K = p.Lf^3

    return F ./ (G * lambda.^3 + H * lambda.^2 + J * lambda .+ K)
end

"""
$(TYPEDSIGNATURES)

Calculate condenstation force in a ring.

Since the condenstation force only depends on the total number of overlaps, it is constant
as the ring constricts.
"""
function condensation_force(p::RingParams)
    logarg = 1 + p.KsD^2 * p.cX / (p.KdD * (p.KsD + p.cX)^2)

    return -2pi * kb * p.T * (2p.Nf - p.Nsca) / (p.Nsca * p.deltad) * log.(logarg)
end


"""
$(TYPEDSIGNATURES)

Calculate current entropic force for a ring.
"""
function entropic_force(lambda, Nd, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    l = 1 .+ p.deltas / p.deltad * lambda

    # This is to prevent it from crashing, not a great solution
    if l*overlaps < Nd
        Nd = l*overlaps
    end
    logarg = 1 .- Nd ./ (l * overlaps)

    return overlaps * kb * p.T * log.(logarg) / p.deltad
end


"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker binding quasi-equilibrium.
"""
function friction_coefficient_ring_cX(lambda, p::RingParams)
    zs = p.r01 / p.r10
    zd = p.r01 * p.r12 / (p.r10 * p.r21)
    z = zd / (1 + zs)^2
    rhos = (zs + zs^2) / ((1 + zs)^2 + zd)
    rhod = z / (1 + z)
    B = p.k * p.deltas^2 / (8kb * p.T) - log(2)
    C = (z + 1) / (z * exp.(-B * exp.((rhod + rhos) / (4B))) + 1)

    return p.zeta0 * C.^((1 .+ p.deltas / p.deltad * lambda) * (2p.Nf - p.Nsca))
end


"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.
"""
function friction_coefficient_ring_Nd(lambda, Ndtot, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot ./ overlaps
    B = p.k * p.deltas^2 / (8kb * p.T) - log(2)
    innerexp = Nd ./ ((1 .+ p.deltas / p.deltad * lambda) * 4B)

    return p.zeta0 * exp.(Ndtot * B .* exp.(innerexp))
end

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Does not assume that the spring constant of the crosslinker is larger relative to thermal
fluctations.
"""
function friction_coefficient_single_exp_ring_Nd(lambda, Ndtot, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot ./ overlaps
    B = p.k * p.deltas^2 / (8kb * p.T) - log(2)
    l = 1 .+ p.deltas / p.deltad * lambda

    return p.zeta0 * exp.(Ndtot .* (Nd ./ (4l) .+ B))
end

function friction_coefficient_exact_Nd_integrand(Nd, l, p)
    function parameterized_integrand(NR)
        bt = binomialc(l - NR, Nd - NR) * binomialc(l - Nd + NR, NR) / binomialc(l, Nd)
        et = exp(-p.k * p.deltas^2 * Nd / (2kb * p.T) * (3 / Nd^2 * (NR - Nd / 2)^2 + 1 / 4))
        return bt * et
    end
    return parameterized_integrand
end

function integrate(integrand, a, b)
    n = 10000
    dx = (b - a) / n
    integral = 0
    x = a
    for _ in 0:n
        integral += integrand(x) * dx
        x += dx
    end

    return integral
end

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Uses the exact expression for the friction coefficient, but with continuous overlap
lenghts and continuous number of bound crosslinkers.
"""
function friction_coefficient_exact_ring_Nd(lambda, Ndtot, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    Nd = Ndtot / overlaps
    l = 1 + p.deltas / p.deltad * lambda
#    result = quadgk(friction_coefficient_exact_Nd_integrand(Nd, l, p), 0, Nd - 1)
    result = integrate(friction_coefficient_exact_Nd_integrand(Nd, l, p), 0, Nd)
#    println(result)
    z_ratio = result[1]
    r0 = kb * p.T / (p.deltas^2 * p.zeta0) * sqrt(1 + 3p.k * p.deltas^2 / (4kb * p.T))

    return (kb * p.T / (p.deltas^2 * r0 * z_ratio))^overlaps
end

function friction_coefficient_exact_Nd_summand(NR, Nd, l, p)
	bt = binomial(l - NR, Nd - NR) * binomial(l - Nd + NR, NR) / binomial(l, Nd)
    et = exp(-p.k * p.deltas^2 * Nd / (2kb * p.T) * (3 / Nd^2 * (NR - Nd / 2)^2 + 1 / 4))
	return bt * et
end;

function sum_NR(Nd, l, p::RingParams)
    z_ratio = 0
    for NR in 0:(Nd - 1)
        z_ratio += friction_coefficient_exact_Nd_summand(NR, Nd, l, p)
    end

    return z_ratio
end

function friction_coefficient_exact_ring_Nd(lambdas::Vector, Ndtots::Vector, p::RingParams)
    zetas = []
    for (lambda, Ndtot) in zip(lambdas, Ndtots)
        push!(zetas, friction_coefficient_exact_ring_Nd(lambda, Ndtot, p))
    end

    return zetas
end

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use discrete N and l.
"""
function friction_coefficient_exact_ring_Nd_discrete(Nd, l, overlaps, p::RingParams)
    z_ratio = sum_NR(Nd, l, p)
    r0 = kb * p.T / (p.deltas^2 * p.zeta0) * sqrt(1 + 3p.k * p.deltas^2 / (4kb * p.T))

    return (kb * p.T / (p.deltas^2 * r0 * z_ratio))^overlaps
end

function friction_coefficient_exact_Nd_continuous_l_summand(NR, Nd, l, p)
    bt = binomialc(l - NR, Nd - NR) * binomialc(l - Nd + NR, NR) / binomialc(l, Nd)
    et = exp(-p.k * p.deltas^2 * Nd / (2kb * p.T) * (3 / Nd^2 * (NR - Nd / 2)^2 + 1 / 4))
    return bt * et
end;

function sum_NR_continuous_l(Nd, l, p::RingParams)
    z_ratio = 0
    for NR in 0:(Nd - 1)
        z_ratio += friction_coefficient_exact_Nd_continuous_l_summand(NR, Nd, l, p)
    end

    return z_ratio
end

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use discrete N and continous l.
"""
function friction_coefficient_exact_discrete_N_ring_Nd(lambda, Nd, p::RingParams)
    l = 1 + p.deltas / p.deltad * lambda

    # This is to prevent it from crashing, not a great solution
    if l < Nd
        Nd = l
    end
    z_ratio = sum_NR_continuous_l(Nd, l, p)
    r0 = kb * p.T / (p.deltas^2 * p.zeta0) * sqrt(1 + 3p.k * p.deltas^2 / (4kb * p.T))

    return kb * p.T / (p.deltas^2 * r0 * z_ratio)
end

function friction_coefficient_exact_discrete_N_ring_Nd(lambdas::Vector, Nds::Vector, p::RingParams)
    zetas = []
    for (lambda, Nd) in zip(lambdas, Nds)
        zeta = 1
        for Ndi in Nd
            zeta *= friction_coefficient_exact_discrete_N_ring_Nd(lambda, Ndi, p)
        end
        push!(zetas, zeta)
    end

    return zetas
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

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

Uses the exact expression for the friction coefficient, but with continuous overlap
lenghts and continuous number of bound crosslinkers. This is compatible with the
DifferentialEquations package.
"""
function equation_of_motion_exact_ring_Nd!(du, u, p, t)
    zeta = friction_coefficient_exact_ring_Nd(u[1], u[2], p)
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
function equation_of_motion_exact_discrete_N_ring_Nd!(du, u, p, t)
    zeta = 1
    for Nd in u[3:end]
        zeta *= friction_coefficient_exact_discrete_N_ring_Nd(u[1], Nd, p)
    end
    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(u[1], p) + entropic_force(u[1], u[2], p)
    ltot = (1 + p.deltas / p.deltad * u[1]) * overlaps

    du[1] = -forcetot / (zeta * p.deltas * overlaps)
    for i in 2:length(du)
        du[i] = 0
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Equation of motion for a ring with crosslinker binding quasi-equlibrium.

This uses the exact expression for the friction coefficient with a discrete number of bound
crosslinkers. This requires the number of crosslinkers per overlap to be tracked, and the
number of bound crosslinkers per overlap to be updated with a jump process separately. This
is compatible with the sample_trajectory method implemented in this package.
"""
function equation_of_motion_exact_discrete_N_ring_Nd(lambda, Nd, p)
    zeta = 1
    for Ndi in Nd
        zeta *= friction_coefficient_exact_discrete_N_ring_Nd(lambda, Ndi, p)
    end
    Ndtot = sum(Nd)
    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(lambda, p) + entropic_force(lambda, Ndtot, p)
    ltot = (1 + p.deltas / p.deltad * lambda) * overlaps
    return -forcetot / (zeta * p.deltas * overlaps)
end

"""
$(TYPEDSIGNATURES)

Generate crosslinker binding rate function.
"""
function binding_rate_generator(i::Integer)
    function binding_rate(u, p, t)
        l = 1 + p.deltas / p.deltad * u[1]
        if trunc(l) == u[i]
#            println("Too many crosslinkers, no binding")
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
#        println(u)
        l = 1 + p.deltas / p.deltad * u[1]
        if u[i] == 1.0
            return 0.0
        elseif trunc(l) < u[i]
            tl = trunc(l)
            Nd = u[i]
#            println("Too many crosslinkers ($l $tl $Nd, high unbinding")
            return 1e20
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
#        println("Event $event")
        integrator.u[2] += event
        integrator.u[i] += event

        return nothing
    end

    return react!
end

"""
$(TYPEDSIGNATURES)

Calculate the equilibrium radius of a ring.
"""
function calc_equilibrium_ring_radius(p::RingParams)
    num = p.EI * p.Nf * p.deltad * p.Lf * p.Nsca
    logarg = 1 + p.KsD^2 * p.cX / (p.KdD * (p.KsD + p.cX)^2)
    denom = (2pi * p.T * kb * log(logarg) * (2 * p.Nf - p.Nsca))

    return (num / denom)^(1 / 3)
end
