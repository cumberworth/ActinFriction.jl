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
Base.@kwdef mutable struct RingParams
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
    """Jump rate prefactor"""
    r0::Float64 = NaN
    """Actin filament diameter (m)"""
    Df::Float64
    """Viscosity of fluid (kg m^-1 s^-1)"""
    eta::Float64
    """Diffusion coefficient of singly bound crosslinkers on filament (m^2 s^-1)"""
    Ds::Float64
    """Dissociation constant for single crosslinker binding from solution (M)"""
    KsD::Float64
    """Dissociation constant for double crosslinker binding from solution (M)"""
    KdD::Float64
    """Crosslinker concentration (M)"""
    cX::Float64
    """Number of overlaps moving collectively during constriction"""
    n::Float64
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

Calculate friction coefficient prefactor with preset r0.
"""
function zeta0(p::RingParams)
    return kb * p.T / p.deltas^2 / p.r0 * sqrt(1 + 3 * p.k * p.deltas^2 / 4 / kb / p.T)
end

"""
$(TYPEDSIGNATURES)

Calculate jump rate prefactor from friction coefficient prefactor.
"""
function zeta0_to_r0(zeta0, p::RingParams)
    return kb * p.T / p.deltas^2 / zeta0 * sqrt(1 + 3 * p.k * p.deltas^2 / 4 / kb / p.T)
end


"""
$(TYPEDSIGNATURES)

Calculate parallel translation diffusion coefficient of a rigid rod.
"""
function rod_parallel_diffusion_coefficient(p::RingParams)
    gamma = -0.114
    kb * p.T * (log(p.Lf / p.Df) - gamma) / (2pi * p.eta * p.Lf)
end

function singly_bound_jump_rate_prefactor(p::RingParams)
    return p.Ds / p.deltas^2
end

"""
$(TYPEDSIGNATURES)

Calculate diffusion coefficient of filament at peak of sliding barrier.
"""
function barrier_diffusion_coefficient(l, Nd, p::RingParams)
    Dm = rod_parallel_diffusion_coefficient(p)
    h0 = singly_bound_jump_rate_prefactor(p)

    return (Dm / p.deltas^2 + h0 / Nd * (1 - Nd / l)) / 4
#    return (Dm / p.deltas^2 + h0 / Nd) / 4
end

"""
$(TYPEDSIGNATURES)

Calculate exact free energy barrier to sliding.
"""
function free_energy_barrier_Nd_exact(l, Nd, p::RingParams)
    return -kb * p.T * log(sum_NR_continuous_l(Nd, l, p))
end

"""
$(TYPEDSIGNATURES)

Calculate jump rate prefactor from Kramers' theory.
"""
function kramers_r0(lambda, Nd, p::RingParams)
    l = lambda_to_l(lambda, p)
    beta = 1 / kb / p.T
    DF = free_energy_barrier_Nd_exact(l, Nd, p)
    D = barrier_diffusion_coefficient(l, Nd, p)

    return 8 * (beta * DF)^3 * D * pi / (8 * (beta * DF)^2 + 4beta * DF + 5)
end

"""
$(TYPEDSIGNATURES)

Convert from force on L to force on R.
"""
function force_L_to_R(force, p::RingParams)
    return force * -2pi ./ (p.Nsca)
end

"""
$(TYPEDSIGNATURES)

Convert from continuous number of sites in an overlap to lambda.
"""
function l_to_lambda(l, p::RingParams)
    return (l .- 1) * p.deltad / p.deltas
end

"""
$(TYPEDSIGNATURES)

Convert from lambda to continuous number of sites in an overlap.
"""
function lambda_to_l(lambda, p::RingParams)
    return 1 .+ p.deltas / p.deltad .* lambda
end

"""
$(TYPEDSIGNATURES)

Convert from lambda to discrete number of sites in an overlap.
"""
function lambda_to_l_discrete(lambda, p::RingParams)
    # return floor(p.deltas / p.deltad * (lambda + 1)) + 1
    return floor(p.deltas / p.deltad * lambda) + 1
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

Calculate current entropic force on L for a ring.
"""
function entropic_force(lambda, Ndtot, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    logarg = 1 .- p.deltad * Ndtot ./ ((p.deltas * lambda .+ p.deltad) * overlaps)

    return -overlaps * kb * p.T / p.deltad * log.(logarg)
end

function friction_coefficient_B(p::RingParams)
    return p.k * p.deltas^2 / (8kb * p.T) - log(2)
end

function friction_coefficient_cX_C(p::RingParams)
    zs = p.cX / p.KsD
    zd = p.cX / p.KdD
    z = zd / (1 + zs)^2
    B = friction_coefficient_B(p)
    C = (z + 1) / (z * exp.(-B) + 1)

    return C
end

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker binding quasi-equilibrium.
"""
function friction_coefficient_cX(lambda, p::RingParams)
    C = friction_coefficient_cX_C(p)

    return zeta0(p) * C.^((1 .+ p.deltas / p.deltad * lambda) * p.n)
end

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.
"""
function friction_coefficient_Nd_exp(lambda, Nd, p::RingParams)
    B = friction_coefficient_B(p)
#    innerexp = Nd ./ ((1 .+ p.deltas / p.deltad * lambda) * 4B)

    return zeta0(p) * exp.(p.n * Nd * B)
end

"""
$(TYPEDSIGNATURES)

Calculate mean friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use discrete N and continous l.
"""
function friction_coefficient_Nd_exact(lambda, Nds, p::RingParams)
    zeta = 0
    for Nd in Nds
        # println("Ndi = $Nd")
        l = lambda_to_l(lambda, p)
        z_ratio = sum_NR_continuous_l(Nd, l, p)
        zeta += kb * p.T / (p.deltas^2 * p.r0 * z_ratio)
    end

    return zeta / length(Nds)
end

function friction_coefficient_Nd_exact(lambdas::Vector, Ndss::Vector, p::RingParams)
    zetas = []
    for (lambda, Nds) in zip(lambdas, Ndss)
        zeta = friction_coefficient_Nd_exact(lambda, Nds, p)
        push!(zetas, zeta)
    end

    return zetas
end

function sum_NR_continuous_l(Nd, l, p::RingParams)
    z_ratio = 0
    for NR in 0:(p.n * Nd)
        z_ratio += friction_coefficient_continuous_l_summand(NR, Nd, l, p)
    end

    return z_ratio
end

function friction_coefficient_continuous_l_summand(NR, Nd, l, p)
#    bt = binomialc(l - NR, Nd - NR) * binomialc(l - Nd + NR, NR) / binomialc(l, Nd)
    bt = binomialc(p.n * Nd, NR)
    et = exp(-p.k * p.deltas^2 * p.n * Nd / (2kb * p.T) *
            (3 / (p.n * Nd)^2 * (NR - p.n * Nd / 2)^2 + 1 / 4))
    return bt * et
end

"""
$(TYPEDSIGNATURES)

Calculate friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use the exact expression (discrete N and l).
"""
function friction_coefficient_Nd_exact_discrete(Nd, l, overlaps, p::RingParams)
    z_ratio = sum_NR(Nd, l, p)

    return (kb * p.T / (p.deltas^2 * p.r0 * z_ratio))^overlaps
end

function sum_NR_discrete_l(Nd, l, p::RingParams)
    z_ratio = 0
    for NR in 0:(Nd)
        z_ratio += friction_coefficient_discrete_l_summand(NR, Nd, l, p)
    end

    return z_ratio
end

function friction_coefficient_discrete_l_summand(NR, Nd, l, p)
#	bt = binomial(l - NR, Nd - NR) * binomial(l - Nd + NR, NR) / binomial(l, Nd)
	bt = binomial(p.n * Nd, NR)
    et = exp(-p.k * p.deltas^2 * p.n * Nd / (2kb * p.T) *
            (3 / (p.n * Nd)^2 * (NR - p.n * Nd / 2)^2 + 1 / 4))
	return bt * et
end;

"""
$(TYPEDSIGNATURES)

Calculate the equilibrium radius of a ring.
"""
function equilibrium_ring_radius(p::RingParams)
    num = p.EI * p.Nf * p.deltad * p.Lf * p.Nsca
    logarg = 1 + p.KsD^2 * p.cX / (p.KdD * (p.KsD + p.cX)^2)
    denom = (2pi * p.T * kb * log(logarg) * (2 * p.Nf - p.Nsca))

    return (num / denom)^(1 / 3)
end
