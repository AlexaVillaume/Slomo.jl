import Slomo.Models: DensityModel, mass, volume_density

"""
Enclosed mass for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function M_NFW(r::Array{Float64, 1}, rs, rhos)::Array{Float64, 1}
    x = r / rs
    factor = 4.0 * pi * rs^3 * rhos
    mu = log.(1.0 .+ x) .- x ./ (1.0 .+ x)
    return factor .* mu
end

"""
Local volume density for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function rho_NFW(r::Array{Float64, 1}, rs, rhos)::Array{Float64, 1}
    x = r / rs
    return rhos * x .^ -1 .* (1 .+ x) .^ -2
end

struct NFWModel <: DensityModel
    rs::Float64
    rhos::Float64
    NFWModel(rs, rhos) = begin
        @assert rs > 0 "rs must be positive"
        @assert rhos > 0 "rhos must be positive"
        new(rs, rhos)
    end
end

function volume_density(model::NFWModel, r::Array{Float64, 1})
    return rho_NFW(r, model.rs, model.rhos)
end

function mass(model::NFWModel, r::Array{Float64, 1})
    return M_NFW(r, model.rs, model.rhos)
end
