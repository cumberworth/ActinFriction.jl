"""
$(TYPEDSIGNATURES)

Calculate continuous binonomial coefficient.
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
    """Duration of dynamics (s)"""
    tend::Float64
    """Initial lambda"""
    lambda0::Float64
    """Initial total doubly-bound crosslinkers"""
    Ndtot0::Float64 = NaN
end

"""
$(TYPEDSIGNATURES)

Generate file name from a set of parameters.

This is a replacement for the DrWatson function that has better handling of rounding and
trailing zeros.
"""
function savename(prefix, params; digits=2, suffix=nothing, ignored_fields=[])
    fields = fieldnames(typeof(params))
    fieldstrings = [string(f) for f in fields]
    sortedindices = sortperm(fieldstrings)
    filename = [prefix]
    for i in sortedindices
        fieldstring = fieldstrings[i]
        value  = getproperty(params, fields[i])
        if fieldstring in ignored_fields
        elseif isinteger(value)
            fstring = Printf.Format("%.i")
            vstring = Printf.format(fstring, value)
            push!(filename, "$(fieldstring)=$(vstring)")
        elseif isnan(value)
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
function entropic_force(lambda, Ndtot, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    l = 1 .+ p.deltas / p.deltad * lambda
    logarg = 1 .- Ndtot ./ (l * overlaps)

    ltot = l * overlaps
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

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use discrete N and continous l.
"""
function friction_coefficient_continuous_l_ring_Nd(lambda, Nds, p::RingParams)
    zeta = 1
    for Nd in Nds
        l = 1 + p.deltas / p.deltad * lambda
        z_ratio = sum_NR_continuous_l_overlap_Nd(Nd, l, p)
        r0 = kb * p.T / (p.deltas^2 * p.zeta0) * sqrt(1 + 3p.k * p.deltas^2 / (4kb * p.T))
        zeta *= kb * p.T / (p.deltas^2 * r0 * z_ratio)
    end

    return zeta
end

function friction_coefficient_continuous_l_ring_Nd(lambdas::Vector, Ndss::Vector, p::RingParams)
    zetas = []
    for (lambda, Nds) in zip(lambdas, Ndss)
        zeta = friction_coefficient_continuous_l_ring_Nd(lambda, Nds, p)
        push!(zetas, zeta)
    end

    return zetas
end

function sum_NR_continuous_l_overlap_Nd(Nd, l, p::RingParams)
    z_ratio = 0
    for NR in 0:(Nd - 1)
        z_ratio += friction_coefficient_continuous_l_overlap_Nd_summand(NR, Nd, l, p)
    end

    return z_ratio
end

function friction_coefficient_continuous_l_overlap_Nd_summand(NR, Nd, l, p)
    bt = binomialc(l - NR, Nd - NR) * binomialc(l - Nd + NR, NR) / binomialc(l, Nd)
    et = exp(-p.k * p.deltas^2 * Nd / (2kb * p.T) * (3 / Nd^2 * (NR - Nd / 2)^2 + 1 / 4))
    return bt * et
end;

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use the exact expression (discrete N and l).
"""
function friction_coefficient_exact_ring_Nd(Nd, l, overlaps, p::RingParams)
    z_ratio = sum_NR(Nd, l, p)
    r0 = kb * p.T / (p.deltas^2 * p.zeta0) * sqrt(1 + 3p.k * p.deltas^2 / (4kb * p.T))

    return (kb * p.T / (p.deltas^2 * r0 * z_ratio))^overlaps
end

function sum_NR(Nd, l, p::RingParams)
    z_ratio = 0
    for NR in 0:(Nd - 1)
        z_ratio += friction_coefficient_exact_overlap_Nd_summand(NR, Nd, l, p)
    end

    return z_ratio
end

function friction_coefficient_exact_overlap_Nd_summand(NR, Nd, l, p)
	bt = binomial(l - NR, Nd - NR) * binomial(l - Nd + NR, NR) / binomial(l, Nd)
    et = exp(-p.k * p.deltas^2 * Nd / (2kb * p.T) * (3 / Nd^2 * (NR - Nd / 2)^2 + 1 / 4))
	return bt * et
end;

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

"""
$(TYPEDSIGNATURES)

Calculate basic quantities.
"""
function calc_basic_quantities(lambda, times, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    x = p.deltas * lambda
    ls = 1 .+ p.deltas ./ p.deltad .* lambda
    ltots = ls .* overlaps
    R = p.Nsca * (p.Lf .- x) / (2pi)
    R_max = p.Nsca * p.Lf / (2pi)
    R_eq = calc_equilibrium_ring_radius(p)
    R_eq_frac = (R_max .- R) / (R_max - R_eq)
    force_bending = bending_force(lambda, p)
    force_condensation = condensation_force(p)
    df = DataFrame(
        t=times,
        lmbda=lambda,
        ls=ls,
        ltots=ltots,
        x=x,
        R=R,
        R_eq_frac=R_eq_frac,
        force_bending=force_bending,
        force_condensation=force_condensation)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-binding quasi-equilibrium.
"""
function calc_cX_quantities(lambda, times, p::RingParams)
    df = calc_basic_quantities(lambda, times, p)
    df[!, :force_total] = df.force_bending .+ df.force_condensation
    df[!, :zeta_cX] = friction_coefficient_ring_cX(lambda, p)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium.
"""
function calc_Nd_quantities(lambda, Ndtot, times, p::RingParams)
    df = calc_basic_quantities(lambda, times, p)
    df[!, :Ndtot] = Ndtot
    df[!, :force_entropy] = entropic_force(lambda, Ndtot, p)
    df[!, :force_total] = df.force_bending .+ df.force_entropy
    df[!, :zeta_Nd_double_exp] = friction_coefficient_ring_Nd(lambda, Ndtot, p)
    df[!, :zeta_Nd_single_exp] = friction_coefficient_single_exp_ring_Nd(lambda, Ndtot, p)

    return df
end

"""
$(TYPEDSIGNATURES)

Calculate quantities for crosslinker-diffusion quasi-equilibrium with discrete Nd.
"""
function calc_discrete_Nd_quantities(lambda, Ndtot, Nds, times, p::RingParams)
    df = calc_Nd_quantities(lambda, Ndtot, times, p)
    zeta_continuous_l = friction_coefficient_continuous_l_ring_Nd(lambda, Nds, p)
    df[!, :zeta_continuous_l] = zeta_continuous_l

    return df
end
