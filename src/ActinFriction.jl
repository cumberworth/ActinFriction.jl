module ActinFriction

using DataFrames
using DiffEqCallbacks
using DocStringExtensions
using JumpProcesses
using Printf
using QuadGK
using SpecialFunctions

const kb = 1.380649e-23

include("expressions.jl")
include("diffeq.jl")

end
