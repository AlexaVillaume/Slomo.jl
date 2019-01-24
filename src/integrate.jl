module Integrate

import DifferentialEquations: solve, ODEProblem, Tsit5

const alg = Tsit5()

function solve(f, a::Float64, b::Float64; kwargs...)
    solve(ODEProblem((u, p, t) -> f(t), 0.0, (a, b)), alg; kwargs...)
end

"""
Integrate f from lower bound a to upper bound b.
"""
function integrate(f, a::Float64, b::Float64)::Float64
    return solve(f, a, b; save_everystep=false)[end]
end

function integrate(f, a::Float64, b::Array{Float64, 1})::Array{Float64, 1}
    bmax = maximum(b)
    return solve(f, a, bmax; saveat = b).u
end

end
