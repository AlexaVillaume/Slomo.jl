using Interpolations: LinearInterpolation

using Slomo.Models: update, DensityModel, Tracer
using Slomo.Models: mass, density, density2d, K_jeans, g_jeans, beta

include("constants.jl")
include("integrate.jl")

"""
Calculate the projected momements of the line-of-sight velocity distribution by
numerically integrating the spherical Jeans equations.

    R : projected radii in kpc
    n_interp : number of points on the interpolation grid for the moment
    rmax : upper limit of projection integration, in kpc
    parameters : keywords to pass on to mass and tracer models
"""
function jeans(mass_model::DensityModel, tracer::Tracer, R;
               n_interp::Int = 10, rmax::Float64 = 1e5, fudge=1e-6,
               parameters...)

    # update models with new parameters
    p = Dict{Symbol, Float64}(key => value for (key, value) in parameters)
    mass_model = update(mass_model, p)
    tracer = update(tracer, p)

    # set up interpolation grid
    Rmin = minimum(R)
    Rmax = maximum(R)
    Rgrid = 10 .^ collect(range(log10(Rmin), stop=log10(Rmax), length=n_interp))
    sgrid = zeros(n_interp)
    
    # construct required functions
    M(r) = mass(mass_model, r)
    ρ(r) = density(tracer.density_model, r)
    Σ(R) = density2d(tracer.density_model, R)
    # K(r, R) = K_jeans(tracer.anisotropy_model, r, R)
    β(r) = beta(tracer.anisotropy_model, r)
    g(r) = g_jeans(tracer.anisotropy_model, r)
    # construct radial velocity dispersion function
    ρσr2_integrand(r) = -M(r) * ρ(r) * g(r) / r ^ 2
    ρσr2_sol = solve(ρσr2_integrand, rmax, Rmin)
    # account for any numerical noise that fluctuates the integral to be negative
    σr2(r) = begin
        y = ρσr2_sol(r) ./ (ρ(r) .* g(r))
        if typeof(y) <: Array
            mask = y .< 0
            y[mask] = zeros(count_true(mask)) .+ eps()
            return y        
        elseif y < 0.0
            return 0.0
        end
        return y
    end
    
    # integrate for los velocity dispersion
    for (i, Ri) in enumerate(Rgrid)
        integrand(r) = (1.0 - β(r) * Ri^2 / r^2) * ρ(r) * σr2(r) * r / √(r^2 - Ri^2)
        sol = solve(integrand, rmax, Ri * (1.0 + fudge))
        integral = sol[1] - sol[end]
        sgrid[i] = √(2 * G / Σ(Ri) * integral)
    end

    s_itp = LinearInterpolation(Rgrid, sgrid)
    return s_itp(R)
end

