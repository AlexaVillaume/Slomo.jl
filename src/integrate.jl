using OrdinaryDiffEq: ODEProblem, solve, Tsit5

const alg = Tsit5()

"""
Integrate f from lower bound a to upper bound b.
"""
function integrate(f, a::Float64, b::Float64)::Float64
    prob = ODEProblem((u, p, t) -> f(t), f(a), (a, b))
    sol = solve(prob, alg)
    return sol[end]
end

function integrate(f, a::Float64, b::Array{Float64, 1})::Array{Float64, 1}
    bmax = maximum(b)
    prob = ODEProblem((u, p, t) -> f(t), f(a), (a, bmax))
    sol = solve(prob, alg)
    return sol(b)
end

function integrate(f, a::Array{Float64, 1}, b::Float64)::Array{Float64, 1}
    # First integrate over the widest range (from the minimum lower bound)
    amin = minimum(a)
    prob = ODEProblem((u, p, t) -> f(t), f(amin), (amin, b))
    ymin = solve(prob, alg)[1]
    # Then subtract the integral from the minimum lower bound up
    amax = maximum(a)
    prob = ODEProblem((u, p, t) -> f(t), f(amin), (amin, amax))
    sol = solve(prob, alg)
    return ymin - sol(a)
end
    

