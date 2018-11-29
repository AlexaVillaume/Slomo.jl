using DifferentialEquations: ODEProblem, solve
using Interpolations: LinearInterpolation

using Slomo.Models: update, DensityModel, Tracer
using Slomo.Models: mass, volume_density, surface_density, K_jeans

include("constants.jl")

"""
Calculate the projected momements of the line-of-sight velocity distribution by
numerically integrating the spherical Jeans equations.

    R : projected radii in kpc
"""
function jeans(mass_model::DensityModel, tracer::Tracer, R::Array{Float64, 1};
               n_interp::Int = 20, max_factor::Float64 = 100.0,
               parameters...)
    p = Dict{Symbol, Float64}(key => value for (key, value) in parameters)
    mass_model = update(mass_model, p)
    tracer = update(tracer, p)
    
    M(r) = mass(mass_model, r)
    rho(r) = volume_density(tracer.density_model, r)
    sigma(R) = surface_density(tracer.density_model, R)
    K(r, R) = K_jeans(tracer.anisotropy_model, r, R)

    Rmin = minimum(R)
    Rmax = maximum(R)
    Rgrid = 10 .^ collect(range(log10(Rmin), stop=log10(Rmax), length=n_interp))
    sgrid = zeros(n_interp)
    
    for (i, R) in enumerate(Rgrid)
        integrand(y, p, r) = K(r, R) .* rho(r) .* M(r) ./ r
        y0 = integrand(0.0, 0.0, R)
        prob = ODEProblem(integrand, y0, (R, 100 * Rmax))
        sol = solve(prob)
        sgrid[i] = 2 * G / sigma(R) * sol[end]
    end
    s_itp = LinearInterpolation(Rgrid, sgrid)
    return s_itp(R)
end
