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
    integral = zero(a)
    for (i, ai) in enumerate(a)
        prob = ODEProblem((u, p, t) -> f(t), f(ai), (ai, b))
        sol = solve(prob, alg)
        integral[i] = sol[end]
    end
end
    

