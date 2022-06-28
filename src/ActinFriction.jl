module ActinFriction

using Printf

const kb = 1.380649e-23

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

function bending_force(lambda, p::RingParams)
    F = 8 * pi^3 * p.EI * p.Lf * p.Nf / p.Nsca^3
    G = -(p.deltas^3)
    H = 3 * p.Lf * p.deltas^2
    J = -3 * p.Lf^2 * p.deltas
    K = p.Lf^3

    return F ./ (G * lambda.^3 + H * lambda.^2 + J * lambda .+ K)
end

function condensation_force(p::RingParams)
    logarg = 1 + p.KsD^2 * p.cX / (p.KdD * (p.KsD + p.cX)^2)

    return -2pi * kb * p.T * (2p.Nf - p.Nsca) / (p.Nsca * p.deltad) * log.(logarg)
end


function entropic_force(lambda, Nd, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    logarg = 1 .- Nd ./ ((1 .+ p.deltas / (p.deltad * overlaps) * lambda) * overlaps)

    return overlaps * kb * p.T * log.(logarg) / p.deltad
end


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


function friction_coefficient_ring_Nd(lambda, Nd, p::RingParams)
    overlaps = 2p.Nf - p.Nsca
    B = p.k * p.deltas^2 / (8kb * p.T) - log(2)
    innerexp = Nd ./ ((1 .+ p.deltas / (p.deltad * overlaps) * lambda) * overlaps * 4B)

    return p.zeta0 * exp.(Nd .* B * exp.(innerexp))
end

function equation_of_motion_ring_cX!(du, u, p, t)
    zeta = friction_coefficient_ring_cX(u[1], p)
    forcetot = bending_force(u[1], p) + condensation_force(p)

    du[1] = -forcetot / (zeta * p.deltas * (2p.Nf - p.Nsca))
end

function equation_of_motion_ring_Nd!(du, u, p, t)
    zeta = friction_coefficient_ring_Nd(u[1], u[2], p)
    overlaps = 2p.Nf - p.Nsca
    forcetot = bending_force(u[1], p) + entropic_force(u[1], u[2], p)
    ltot = (1 + p.deltas / p.deltad * u[1]) * overlaps
    du[1] = -forcetot / (zeta * p.deltas * overlaps)
    du[2] = p.cX * p.k01 * p.r12 * ltot - (p.cX * p.k01 * p.r12 - p.r21 * p.r10) * u[2]
end

function calc_equilibrium_ring_radius(p::RingParams)
    num = p.EI * p.Nf * p.deltad * p.Lf * p.Nsca
    logarg = 1 + p.KsD^2 * p.cX / (p.KdD * (p.KsD + p.cX)^2)
    denom = (2pi * p.T * kb * log(logarg) * (2 * p.Nf - p.Nsca))

    return (num / denom)^(1 / 3)
end

end
