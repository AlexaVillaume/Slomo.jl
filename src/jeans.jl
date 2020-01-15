using Interpolations: LinearInterpolation
using DifferentialEquations

using Slomo.Models: update, JeansModel, NotImplemented, has_analytic_profile
using Slomo.Models: mass, density, density2d, K_jeans, g_jeans, beta
using Slomo.Constants: G
using Slomo.Integrate: solve, integrate

"""
Calculate the line-of-sight velocity dispersion by
numerically integrating the spherical Jeans equations.

    model : potential and tracer models
    R : projected radii in kpc
    n_interp : number of points on the interpolation grid for the moment
    
    parameters : keywords to pass on to mass and tracer models
"""
function sigma_los(model::JeansModel, R;
                   n_interp::Int = 10, fudge = 1e-6, interp = true, rmax = 1e2,
                   parameters...)
    
    # update models with new parameters
    if length(keys(parameters)) > 0
        model = update(model; parameters...)
    end
    mass_model = model.mass_model
    density_model = model.density_model
    anisotropy_model = model.anisotropy_model

    # check for analytic functions
    has_kernel = has_analytic_profile(K_jeans, anisotropy_model)
    # has_g = has_analytic_profile(g_jeans, anisotropy_model)
    # has_mass = has_analytic_profile(mass, mass_model)
    # has_sd = has_analytic_profile(density2d, density_model)

    # set up interpolation grid
    Rmin = minimum(R)
    Rmax = maximum(R)
    if (length(R) <= n_interp) || ~interp
        Rgrid = R
        interp = false
    else
        Rgrid = exp10.(collect(range(log10(Rmin), stop=log10(Rmax), length=n_interp)))
    end

    # construct required functions
    M(r) = mass(mass_model, r)
    ρ(r) = density(density_model, r)
    Σ(R) = density2d(density_model, R)
    K(r, R) = K_jeans(anisotropy_model, r, R)
    β(r) = beta(anisotropy_model, r)
    g(r) = g_jeans(anisotropy_model, r)
    
    # if we have an analytic Jeans kernel, use it!
    if has_kernel
        integral = integral_from_kernel(Rgrid, M, ρ, K; fudge = fudge,
                                        rmax = rmax)
    else
        integral = integral_from_vr2(Rgrid, M, ρ, β, g; fudge = fudge,
                                     rmax = rmax)
    end

    # restrict to non-negative quantities
    integral[integral .< eps()] .= 0
    sgrid = @. √(2 / Σ(Rgrid) * integral)
    if interp
        return LinearInterpolation(Rgrid, sgrid)(R)
    else
        return sgrid
    end
end

"""
Use multi-threading to map out the integration to different threads.

parameter_sets should be an array of dictionaries containing parameter values
to use for the integration
"""
function sigma_los_parallel(model::JeansModel, R,
                            parameter_sets::Array{Dict{Symbol,T}} where T<:Number;
                            n_interp::Int = 10, fudge = 1e-6, interp = true, rmax = 1e2)
    sigma_profiles = zeros(length(parameter_sets), length(R))
    Threads.@threads for i = 1:length(parameter_sets)
        sigma_profiles[i, :] = sigma_los(model, R;
                                         n_interp = n_interp, fudge = fudge,
                                         interp = interp, rmax = rmax,
                                         parameter_sets[i]...)
    end
    return sigma_profiles
end
                            
"""
Calculate the line-of-sight excess kurtosis by numerically integrating the 
spherical Jeans equations.

    model : potential and tracer models
    R : projected radii in kpc
    parameters : keywords to pass on to mass and tracer models
    return_sigma : if true, will return both excess kurtosis and dispersion
"""
function kappa_los(model::JeansModel, R;
                   n_interp::Int = 10, fudge = 1e-6, interp = true,
                   return_sigma = false, rmax = 1e2,
                   parameters...)
    # update models with new parameters
    if length(keys(parameters)) > 0
        model = update(model; parameters...)
    end
    mass_model = model.mass_model
    density_model = model.density_model
    anisotropy_model = model.anisotropy_model

    # set up interpolation grid
    Rmin = minimum(R)
    Rmax = maximum(R)
    if (length(R) <= n_interp) || ~interp
        Rgrid = R
    else
        Rgrid = exp10.(collect(range(log10(Rmin), stop=log10(Rmax), length=n_interp)))
    end

    # construct required functions
    M(r) = mass(mass_model, r)
    ρ(r) = density(density_model, r)
    Σ(R) = density2d(density_model, R)
    β(r) = beta(anisotropy_model, r)
    g(r) = g_jeans(anisotropy_model, r)

    rmin = minimum(Rgrid)
    
    # integrand for ρ * g * v_r^2
    ρgvr2_integrand(r) = -G * M(r) * ρ(r) * g(r) / r ^ 2
    prob = ODEProblem((u, p, t) -> ρgvr2_integrand(t), 0.0, (rmax, rmin * 1e-2))
    ρgvr2_sol = solve(prob, Tsit5() , isoutofdomain=(u, p, t) -> any(u .<= 0))

    # capture bad inputs, restrict output to non-negative values
    ρgvr2(r) = begin
        if r > rmax
            return 0.0
        end
        return max(ρgvr2_sol(r), 0.0)
    end
    
    # calculate the second moment of the losvd
    integral = zeros(size(Rgrid))
    for (i, R) in enumerate(Rgrid)
        integrand(r) = (1.0 - β(r) * (R / r)^2) * ρgvr2(r) / g(r) * r / √(r^2 - R^2)
        sol = solve(integrand, rmax, R * (1.0 + fudge))
        integral[i] = sol[1] - sol[end]
    end
    sgrid = @. √(2 / Σ(Rgrid) * integral)

    # integrand for ρ * g * v_r^4 / G
    ρgvr4_integrand(r) = -3 * ρgvr2(r) * G * M(r) / r ^ 2
    prob = ODEProblem((u, p, t) -> ρgvr4_integrand(t), 0.0, (rmax, rmin * 1e-1))
    ρgvr4_sol = solve(prob, Tsit5() , isoutofdomain=(u, p, t) -> any(u .<= 0))
   
    # capture bad inputs, restrict output to non-negative values
    ρvr4(r) = begin
        if r > rmax
            return 0.0
        end
        y = ρgvr4_sol(r) / g(r)
        return max(y, 0.0) 
    end
    
    # calculate fourth moment of the losvd
    integral = zeros(size(Rgrid))
    for (i, R) in enumerate(Rgrid)
        integrand(r) = begin
            factor1 = (1 - 2 * β(r) * (R / r)^2 + β(r) * (1 + β(r)) / 2 * (R / r)^4)
            factor2 = ρvr4(r) * r / √(r^2 - R^2)
            factor1 * factor2
        end
        sol = solve(integrand, rmax, R * (1.0 + fudge))
        integral[i] = sol[1] - sol[end]
    end
    v4los = @. 2 / Σ(Rgrid) * integral
    kgrid = @. v4los / sgrid ^ 4 - 3
    
    if return_sigma
        return kgrid, sgrid
    else
        return kgrid
    end
end
       
function integral_from_kernel(Rgrid, M, ρ, K; fudge = 1e-6, rmax = 1e2)
    integral = zeros(size(Rgrid))
    for (i, Ri) in enumerate(Rgrid)
        integrand(r) = K(r, Ri) * G * M(r) * ρ(r) / r
        sol = solve(integrand, rmax, Ri * (1.0 + fudge))
        integral[i] = sol[1] - sol[end]
    end
    return integral    
end

function integral_from_vr2(Rgrid::Array{Float64, 1}, M, ρ, β, g; fudge = 1e-6, rmax = 1e2)
    
    rmin = minimum(Rgrid)
    
    # integrand for ρ * g * v_r^2
    ρgvr2_integrand(r) = -G * M(r) * ρ(r) * g(r) / r ^ 2
    prob = ODEProblem((u, p, t) -> ρgvr2_integrand(t), 0.0, (rmax, rmin * 1e-2))
    ρgvr2_sol = solve(prob, Tsit5() , isoutofdomain=(u, p, t) -> any(u .<= 0))

    ρgvr2(r) = begin
        if r <= rmax
            return max(ρgvr2_sol(r), 0.0)
        else
            return 0.0
        end
    end
    
    integral = zeros(size(Rgrid))
    for (i, R) in enumerate(Rgrid)
        integrand(r) = (1.0 - β(r) * (R / r)^2) * ρgvr2(r) / g(r) * r / √(r^2 - R^2)
        sol = solve(integrand, rmax, R * (1.0 + fudge))
        integral[i] = sol[1] - sol[end]
    end
    return integral
end

