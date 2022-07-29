module ActinFriction

using DocStringExtensions
using Printf
using QuadGK
using SpecialFunctions

const kb = 1.380649e-23

include("expressions.jl")
include("sampler.jl")

end
