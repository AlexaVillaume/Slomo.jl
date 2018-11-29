using GSL: sf_gamma

import Slomo.Models: DensityModel, mass, surface_density, volume_density

"""
Surface density for Sersic model.

    R : radii in kpc
    Re : effective radius in kpc
    n : index
"""
function s_sersic(R::Array{Float64, 1}, Re, n)::Array{Float64, 1}
    bn = -1. / 3 + 2. * n + 4 / (405. * n) + 46 / (25515. * n^2)
    return exp.(-bn * (R / Re) .^ (1.0 / n))
end


"""
Volume density for Sersic model.

    r : radii in kpc
    Re : effective radius in kpc
    n : index
"""
function rho_sersic(r::Array{Float64, 1}, Re, n)::Array{Float64, 1}
    bn = -1. / 3 + 2. * n + 4 / (405. * n) + 46 / (25515. * n^2)
    pn = 1. - 0.6097 / n + 0.05463 / n^2
    x = r / Re
    nu0 = bn ^ n * sf_gamma(2 * n) / (2 * Re * sf_gamma((3 - pn) * n))
    return nu0 * (bn ^ n * x) .^ (-pn) .* exp.(-bn * x .^ (1.0 / n))
end

struct SersicModel <: DensityModel
    Re::Float64
    n::Float64
    SersicModel(Re, n) = begin
        @assert Re > 0 "Re must be positive"
        @assert n > 0 "n must be positive"
        new(Re, n)
    end
end

function surface_density(model::SersicModel, R::Array{Float64, 1})
    return s_sersic(R, model.Re, model.n)
end

function volume_density(model::SersicModel, r::Array{Float64, 1})
    return rho_sersic(r, model.Re, model.n)
end


    
