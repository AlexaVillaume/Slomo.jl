using Interpolations: LinearInterpolation

using Slomo.Models: update, DensityModel, Tracer
using Slomo.Models: mass, density, density2d, K_jeans, g_jeans, beta

include("constants.jl")
include("integrate.jl")

"""
Construct a function to compute the squared radial velocity dispersion.
"""
function func_sigmar2(mass_model::DensityModel, tracer::Tracer,
                      rmin::Float64; rmax::Float64 = 1e3)
    M(r) = mass(mass_model, r)
    ρ(r) = density(tracer.density_model, r)
    g(r) = g_jeans(tracer.anisotropy_model, r)
    integrand(r) = M(r) * ρ(r) * g(r) / r ^ 2
    sol = solve(integrand, rmin, rmax)
    return r -> sol(r) / (ρ(r) * g(r))
end

"""
Calculate the projected momements of the line-of-sight velocity distribution by
numerically integrating the spherical Jeans equations.

    R : projected radii in kpc
    n_interp : number of points on the interpolation grid for the moment
    max_factor : upper limit of projection integration, in units of the maximum radius
    parameters : keywords to pass on to mass and tracer models
"""
function jeans(mass_model::DensityModel, tracer::Tracer, R;
               n_interp::Int = 20, max_factor::Float64 = 100.0,
               parameters...)
    
    p = Dict{Symbol, Float64}(key => value for (key, value) in parameters)
    mass_model = update(mass_model, p)
    tracer = update(tracer, p)

    compute_with_kernel = hasmethod(K_jeans, (typeof(tracer.anisotropy_model), Any, Any))
    
    Rmin = minimum(R)
    Rmax = maximum(R)
    rmax = max_factor * Rmax
    Rgrid = 10 .^ collect(range(log10(Rmin), stop=log10(Rmax), length=n_interp))
    sgrid = zeros(n_interp)
    
    M(r) = mass(mass_model, r)
    ρ(r) = density(tracer.density_model, r)
    Σ(R) = density2d(tracer.density_model, R)
    K(r, R) = K_jeans(tracer.anisotropy_model, r, R)
    β(r) = beta(tracer.anisotropy_model, r)
    σr2 = func_sigmar2(mass_model, tracer, Rmin; rmax = rmax)
    
    for (i, Ri) in enumerate(Rgrid)
        integrand = compute_with_kernel ?
            integrand(r) = K(r, Ri) * ρ(r) * M(r) / r :
            integrand(r) = (1.0 - β(r) * (Ri / r) ^ 2) * ρ(r) * σr2(r) * r / sqrt(r^2 - Ri^2)
        integral = integrate(integrand, Ri, rmax)
        sgrid[i] = √(2 * G / Σ(Ri) * integral)
    end

    s_itp = LinearInterpolation(Rgrid, sgrid)
    return s_itp(R)
end

