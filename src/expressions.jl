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
    k01::Float64 = NaN
    """Per site rate of initial crosslinker binding (M^-1 s^-1)"""
    r01::Float64 = NaN
    """Per site rate (constant) of singly bound crosslinker unbinding (s^-1)"""
    r10::Float64 = NaN
    """Per site rate (constant) of singly bound crosslinker doubly binding (s^-1)"""
    r12::Float64 = NaN
    """Per site rate (constant) of doubly bound crosslinker unbinding one head (s^-1)"""
    r21::Float64 = NaN
    """Spacing between binding sites on a filament (m)"""
    deltas::Float64 = NaN
    """Spacing for doubly bound crosslinkers (m)"""
    deltad::Float64 = NaN
    """Spring constant of crosslinker (N m^-1)"""
    k::Float64 = NaN
    """Temperature (K)"""
    T::Float64 = NaN
    """Number of filaments"""
    Nf::Int = 0
    """Number of scaffold filaments"""
    Nsca::Int = 0
    """ Bending rigidity (N m^2)"""
    EI::Float64 = NaN
    """Filament length (m)"""
    Lf::Float64 = NaN
    """Jump rate prefactor"""
    r0::Float64 = NaN
    """Actin filament diameter (m)"""
    Df::Float64 = NaN
    """Viscosity of fluid (kg m^-1 s^-1)"""
    eta::Float64 = NaN
    """Diffusion coefficient of singly bound crosslinkers on filament (m^2 s^-1)"""
    Ds::Float64 = NaN
    """Dissociation constant for single crosslinker binding from solution (M)"""
    KsD::Float64 = NaN
    """Dissociation constant for double crosslinker binding from solution (M)"""
    KdD::Float64 = NaN
    """Crosslinker concentration (M)"""
    cX::Float64 = NaN
    """Number of overlaps moving collectively during constriction"""
    n::Float64 = NaN
    """Duration of dynamics (s)"""
    tend::Float64 = NaN
    """Initial lambda"""
    lambda0::Float64 = NaN
    """Initial total doubly-bound crosslinkers"""
    Ndtot0::Float64 = NaN
    """Write inteval for mean and variance of stochastic simulations"""
    interval::Float64 = NaN
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
function barrier_diffusion_coefficient(Nd, p::RingParams)
    Dm = rod_parallel_diffusion_coefficient(p)
    h0 = singly_bound_jump_rate_prefactor(p)

    return (Dm / p.deltas^2 + h0 / Nd) / 4
end

"""
$(TYPEDSIGNATURES)

Calculate exact free energy barrier to sliding in units of kb T.
"""
function free_energy_barrier_cX(l, p::RingParams)
    A = 0.5 * log(1 + 3 * p.k * p.deltas^2 / (4 * kb * p.T))
    C = friction_coefficient_cX_C(p)

    return A + l * log(C)
end

"""
$(TYPEDSIGNATURES)

Calculate exact free energy barrier to sliding in units of kb T.
"""
function free_energy_barrier_Nd_exp(Nd, p::RingParams)
    A = 0.5 * log(1 + 3 * p.k * p.deltas^2 / (4 * kb * p.T))
    B = friction_coefficient_B(p)

    return A + p.n * Nd * B
end

"""
$(TYPEDSIGNATURES)

Calculate exact free energy barrier to sliding in units of kb T.
"""
function free_energy_barrier_Nd_exact(Nd, p::RingParams)
    return -log(sum_NR(Nd, p))
end

"""
$(TYPEDSIGNATURES)

Calculate jump rate prefactor from Kramers' theory.
"""
function kramers_r0(Nd, p::RingParams)
    beta = 1 / kb / p.T
    DF = kb * p.T * free_energy_barrier_Nd_exact(Nd, p)
    D = barrier_diffusion_coefficient(Nd, p)

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
    return floor.(p.deltas / p.deltad * lambda) .+ 1
end

"""
$(TYPEDSIGNATURES)

Convert from ring radius to lambda.
"""
function lambda_to_R(lambda, p::RingParams)
    return p.Nsca / (2pi) * (p.Lf .- p.deltas * lambda)
end

"""
$(TYPEDSIGNATURES)

Convert from ring radius to lambda.
"""
function R_to_lambda(R, p::RingParams)
    return 1 / p.deltas * (p.Lf .- 2pi * R ./ p.Nsca)
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
function friction_coefficient_Nd_exp(Nd, p::RingParams)
    B = friction_coefficient_B(p)

    return zeta0(p) * exp.(p.n * Nd * B)
end

"""
$(TYPEDSIGNATURES)

Calculate mean friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use discrete N.
"""
function friction_coefficient_Nd_exact_ave(Nds, p::RingParams)
    zeta = 0
    for Nd in Nds
        # println("Ndi = $Nd")
        z_ratio = sum_NR(Nd, p)
        zeta += kb * p.T / (p.deltas^2 * p.r0 * z_ratio)
    end

    return zeta / length(Nds)
end

function friction_coefficient_Nd_exact_ave(Ndss::Vector, p::RingParams)
    zetas = []
    for Nds in Ndss
        zeta = friction_coefficient_Nd_exact(Nds, p)
        push!(zetas, zeta)
    end

    return zetas
end

"""
$(TYPEDSIGNATURES)

Calculate mean friction coefficient for a ring with crosslinker diffusion quasi-equilibrium.

Use discrete N.
"""
function friction_coefficient_Nd_exact(Nd, p::RingParams)
    zeta = 0
    # println("Ndi = $Nd")
    z_ratio = sum_NR(Nd, p)

    return kb * p.T / (p.deltas^2 * p.r0 * z_ratio)
end

function friction_coefficient_Nd_exact(Ndss::Vector, p::RingParams)
    zetas = []
    for Nds in Ndss
        zeta = friction_coefficient_Nd_exact(Nds, p)
        push!(zetas, zeta)
    end

    return zetas
end

function sum_NR(Nd, p::RingParams)
    z_ratio = 0
    for NR in 0:(p.n * Nd)
        z_ratio += friction_coefficient_summand(NR, Nd, p)
    end

    return z_ratio
end

function friction_coefficient_summand(NR, Nd, p)
    bt = binomialc(p.n * Nd, NR)
    et = exp(-p.k * p.deltas^2 * p.n * Nd / (2kb * p.T) *
            (3 / (p.n * Nd)^2 * (NR - p.n * Nd / 2)^2 + 1 / 4))
    return bt * et
end

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

"""
$(TYPEDSIGNATURES)

Calculate the equilibrium occupancy.
"""
function equilibrium_occupancy(p::RingParams)
    xi_d = p.cX / p.KdD
    xi_s = p.cX / p.KsD

    return xi_d / ((1 + xi_s)^2 + xi_d)
end
