using Interpolations: LinearInterpolation

using Slomo.Models: update, DensityModel, Tracer
using Slomo.Models: mass, density, density2d, K_jeans

include("constants.jl")
include("integrate.jl")

"""
Calculate the projected momements of the line-of-sight velocity distribution by
numerically integrating the spherical Jeans equations.

    R : projected radii in kpc
    n_interp : number of points on the interpolation grid for the moment
    max_factor : upper limit of projection integration, in units of the maximum radius
    parameters : keywords to pass on to mass and tracer models
"""
function jeans(mass_model::DensityModel, tracer::Tracer, R::Array{Float64, 1};
               n_interp::Int = 20, max_factor::Float64 = 100.0,
               parameters...)
    
    p = Dict{Symbol, Float64}(key => value for (key, value) in parameters)
    mass_model = update(mass_model, p)
    tracer = update(tracer, p)
    
    M(r) = mass(mass_model, r)
    ρ(r) = density(tracer.density_model, r)
    Σ(R) = density2d(tracer.density_model, R)
    K(r, R) = K_jeans(tracer.anisotropy_model, r, R)
    
    Rmin = minimum(R)
    Rmax = maximum(R)
    Rgrid = 10 .^ collect(range(log10(Rmin), stop=log10(Rmax), length=n_interp))
    sgrid = zeros(n_interp)
    for (i, Ri) in enumerate(Rgrid)
        integrand(r) = K(r, Ri) * ρ(r) * M(r) / r
        integral = integrate(integrand, Ri, max_factor * Rmax)
        sgrid[i] = √(2 * G / Σ(Ri) * integral)
    end
    s_itp = LinearInterpolation(Rgrid, sgrid)
    return s_itp(R)
    
end
