module ActinFriction

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
