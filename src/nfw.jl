import Slomo.Models: DensityModel, mass, density

"""
Enclosed mass for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function M_NFW(r, rs::Float64, rhos::Float64)
    x = r / rs
    mu = @. log(1.0 + x) - x / (1.0 + x)
    return 4pi * rs^3 * rhos * mu
end

"""
Local volume density for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function rho_NFW(r, rs::Float64, rhos::Float64)
    x = r / rs
    return rhos * (@. x ^ -1 * (1 + x) ^ -2)
end

struct NFWModel <: DensityModel
    rs::Float64
    rhos::Float64
end

NFWModel() = NFWModel(50.0, 1e7)

function density(model::NFWModel, r)
    return rho_NFW(r, model.rs, model.rhos)
end

function mass(model::NFWModel, r)
    return M_NFW(r, model.rs, model.rhos)
end
