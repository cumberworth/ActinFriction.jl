"""$(README)

$(EXPORTS)
"""
module ActinFriction

export RingParams, savename
export free_energy_barrier_cX, free_energy_barrier_Nd_exp, free_energy_barrier_Nd_exact
export kramers_r0
export force_L_to_R, l_to_lambda, lambda_to_l, lambda_to_R, R_to_lambda
export bending_force, condensation_force, entropic_force
export friction_coefficient_cX, friction_coefficient_Nd_exp, friction_coefficient_Nd_exact
export equilibrium_ring_radius, equilibrium_occupancy
export solve_and_write_ring_cX
export solve_and_write_ring_cX_noise
export solve_and_write_ring_Nd_contin_exp
export solve_and_write_ring_Nd_contin_exp_noise
export solve_and_write_ring_Nd_discrete_exp
export solve_and_write_ring_Nd_discrete_exp_noise
export solve_and_write_ring_Nd_discrete_exact
export solve_and_write_ring_Nd_discrete_exact_noise

using CSV
using DataFrames
using DifferentialEquations
using DocStringExtensions
using Interpolations
using Printf
using SpecialFunctions
using Statistics

const kb = 1.380649e-23

include("expressions.jl")
include("analysis.jl")
include("diffeq.jl")

end
