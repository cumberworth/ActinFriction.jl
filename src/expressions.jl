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
    r0::Float64
    """Jump rate prefactor"""
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

Calculate friction coefficient prefactor.
"""
function zeta0(p::RingParams)
    return kb * p.T / p.deltas^2 / p.r0 * sqrt(1 + 3 * p.k * p.deltas^2 / 4 / kb / p.T)
end

"""
$(TYPEDSIGNATURES)

Calculate jump rate prefactor.
"""
function r0(zeta0, p::RingParams)
    return kb * p.T / p.deltas^2 / zeta0 * sqrt(1 + 3 * p.k * p.deltas^2 / 4 / kb / p.T)
end

"""
$(TYPEDSIGNATURES)

Convert from force on L to force on R
"""
function force_L_to_R(force, p::RingParams)
    return force * -2pi ./ (p.Nsca)
end

"""
$(TYPEDSIGNATURES)

Calculate current bending force on L in a ring.
"""
function bending_force(lambda, p::RingParams)
    return -4pi^2 * p.EI * p.Lf * p.Nf ./ (p.Nsca^2 * (p.Lf .- p.deltas * lambda).^3)
end

"""
$(TYPEDSIGNATURES)

Calculate condenstation force on L in a ring.

Since the condenstation force only depends on the total number of overlaps, it is constant
as the ring constricts.
"""
function condensation_force(p::RingParams)
    logarg = 1 + p.KsD^2 * p.cX / (p.KdD * (p.KsD + p.cX)^2)

    return kb * p.T * (2p.Nf - p.Nsca) / p.deltad * log.(logarg)
end

"""
$(TYPEDSIGNATURES)

Calculate condenstation force on L in a ring, ignoring singly bound crosslinkers.

Since the condenstation force only depends on the total number of overlaps, it is constant
as the ring constricts.
"""
function condensation_force_ignore_Ns(p::RingParams)
    logarg = 1 + p.cX / p.KdD

    return kb * p.T * (2p.Nf - p.Nsca) / p.deltad * log.(logarg)
end

"""
$(TYPEDSIGNATURES)

Calculate current entropic force on L for a ring.
"""
function entropic_force(lambda, Ndtot, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    logarg = 1 .- p.deltad * Ndtot ./ ((p.deltas * lambda .+ p.deltad) * overlaps)

    return -overlaps * kb * p.T / p.deltad * log.(logarg)
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

    return zeta0(p) * C.^((1 .+ p.deltas / p.deltad * lambda) * (2p.Nf - p.Nsca))
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

    return zeta0(p) * exp.(Ndtot * B .* exp.(innerexp))
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
        zeta *= kb * p.T / (p.deltas^2 * p.r0 * z_ratio)
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

    return (kb * p.T / (p.deltas^2 * p.r0 * z_ratio))^overlaps
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

Calculate the equilibrium radius of a ring, ignoring singly bound crosslinkers.
"""
function calc_equilibrium_ring_radius_ignore_Ns(p::RingParams)
    num = p.EI * p.Nf * p.deltad * p.Lf * p.Nsca
    logarg = 1 + p.cX / p.KdD
    denom = (2pi * p.T * kb * log(logarg) * (2 * p.Nf - p.Nsca))

    return (num / denom)^(1 / 3)
end
