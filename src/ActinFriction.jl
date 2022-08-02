module ActinFriction

using DataFrames
using DifferentialEquations
using DocStringExtensions
using Interpolations
using Printf
using QuadGK
using SpecialFunctions
using Statistics

const kb = 1.380649e-23

include("expressions.jl")
include("diffeq.jl")
include("analysis.jl")

end
