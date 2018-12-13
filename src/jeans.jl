using Interpolations: LinearInterpolation

using Slomo.Models: update, JeansModel, NotImplemented, has_analytic_profile
using Slomo.Models: mass, density, density2d, K_jeans, g_jeans, beta
using Slomo.Constants: G, rmax
using Slomo.Integrate: solve, integrate

"""
Calculate the line-of-sight velocity dispersion by
numerically integrating the spherical Jeans equations.

    model : potential and tracer models
    R : projected radii in kpc
    n_interp : number of points on the interpolation grid for the moment
    rmax : upper limit of projection integration, in kpc
    parameters : keywords to pass on to mass and tracer models
"""
function sigma_los(model::JeansModel, R;
                   n_interp::Int = 10, fudge = 1e-6,
                   parameters...)

    # update models with new parameters
    p = Dict{Symbol, Float64}(key => value for (key, value) in parameters)
    if length(keys(p)) > 0
        model = update(model, p)
    end
    mass_model = model.mass_model
    density_model = model.density_model
    anisotropy_model = model.anisotropy_model

    # check for analytic functions
    has_kernel = has_analytic_profile(K_jeans, anisotropy_model)
    has_g = has_analytic_profile(g_jeans, anisotropy_model)
    has_mass = has_analytic_profile(mass, mass_model)
    has_sd = has_analytic_profile(density2d, density_model)

    # set up interpolation grid
    Rmin = minimum(R)
    Rmax = maximum(R)
    Rgrid = 10 .^ collect(range(log10(Rmin), stop=log10(Rmax), length=n_interp))

    # construct required functions
    M(r) = mass(mass_model, r)
    ρ(r) = density(density_model, r)
    Σ(R) = density2d(density_model, R)
    K(r, R) = K_jeans(anisotropy_model, r, R)
    β(r) = beta(anisotropy_model, r)
    g(r) = g_jeans(anisotropy_model, r)
    
    # if we have an analytic Jeans kernel, use it!
    if has_kernel
        integral = integral_from_kernel(Rgrid, M, ρ, K; fudge = fudge)
    else
        integral = integral_from_srad(Rgrid, M, ρ, β, g; fudge = fudge)
    end
 
    sgrid = @. √(2 * G / Σ(Rgrid) * integral)
    s_itp = LinearInterpolation(Rgrid, sgrid)
    return s_itp(R)        
end

function integral_from_kernel(Rgrid::Array{Float64, 1},
                              M, ρ, K; fudge = 1e-6)::Array{Float64, 1}
    integral = zeros(size(Rgrid))
    for (i, Ri) in enumerate(Rgrid)
        integrand(r) = K(r, Ri) * M(r) * ρ(r) / r
        sol = solve(integrand, rmax, Ri * (1.0 + fudge))
        integral[i] = sol[1] - sol[end]
    end
    return integral    
end

function integral_from_srad(Rgrid::Array{Float64, 1},
                            M, ρ, β, g; fudge = 1e-6)::Array{Float64, 1}
    # construct radial velocity dispersion function
    ρσr2_integrand(r) = -M(r) * ρ(r) * g(r) / r ^ 2
    ρσr2_sol = solve(ρσr2_integrand, rmax, Rgrid[1])
    # account for any numerical noise that fluctuates the integral to be negative
    σr2(r) = max(ρσr2_sol(r) / (ρ(r) * g(r)), 0.0)
    integral = zeros(size(Rgrid))
    for (i, Ri) in enumerate(Rgrid)
        integrand(r) = (1.0 - β(r) * Ri^2 / r^2) * ρ(r) * σr2(r) * r / √(r^2 - Ri^2)
        sol = solve(integrand, rmax, Ri * (1.0 + fudge))
        integral[i] = sol[1] - sol[end]
    end
    return integral
end
