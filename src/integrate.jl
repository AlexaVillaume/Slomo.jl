import OrdinaryDiffEq: solve
using OrdinaryDiffEq: ODEProblem, Tsit5

const alg = Tsit5()

function solve(f, a::Float64, b::Float64)
    solve(ODEProblem((u, p, t) -> f(t), f(a), (a, b)), alg)
end

"""
Integrate f from lower bound a to upper bound b.
"""
function integrate(f, a::Float64, b::Float64)::Float64
    return solve(f, a, b)[end]
end

function integrate(f, a::Float64, b::Array{Float64, 1})::Array{Float64, 1}
    bmax = maximum(b)
    return solve(f, a, bmax)(b)
end

function integrate(f, a::Array{Float64, 1}, b::Float64)::Array{Float64, 1}
    # First integrate over the widest range (from the minimum lower bound)
    amin = minimum(a)
    ymin = solve(f, amin, b)[1]
    # Then subtract the integral from the minimum lower bound up
    amax = maximum(a)
    return ymin - solve(f, amin, amax)(a)
end
    

