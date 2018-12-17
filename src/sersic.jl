using SpecialFunctions: gamma
using StatsFuns: gammacdf

import Slomo.Models: DensityModel, mass, density, density2d

"""
'b' parameter in the Sersic function, from the Ciotti & Bertin (1999) approximation.

n is the Sersic index
"""
function b_cb(n::Float64)::Float64
    -1. / 3 + 2. * n + 4 / (405. * n) + 46 / (25515. * n^2)
end

"""
'p' parameter in fitting deprojected Sersic function, from Lima Neto et al. (1999).

n is the Sersic index
"""
function p_ln(n::Float64)::Float64
    1. - 0.6097 / n + 0.05463 / n^2
end

"""
Surface density for Sersic model.

    R : projected (2d) radii in kpc
    Re : effective radius in kpc
    n : index
    Mtot : total mass
"""
function s_sersic(R, Re::Float64, n::Float64, Mtot::Float64)
    bn = b_cb(n)
    s0 = Mtot * bn ^ (2 * n) / (2π * n * Re ^ 2 * gamma(2 * n))
    return s0 * @. exp(-bn * (R / Re) ^ (1.0 / n))
end

"""
Volume density for Sersic model.

    r : deprojected (3d) radii in kpc
    Re : effective radius in kpc
    n : index
    Mtot : total mass
"""
function rho_sersic(r, Re::Float64, n::Float64, Mtot::Float64)
    bn = b_cb(n)
    pn = p_ln(n)
    x = r / Re
    rho0 = Mtot * bn ^ (3 * n) / (4π * n * Re ^ 3 * gamma((3 - pn) * n))
    return rho0 * (@. (bn ^ n * x) ^ (-pn) * exp(-bn * x ^ (1.0 / n)))
end

"""
Mass enclosed within deprojected radius for the Sersic model.

    r : deprojected (3d) radii in kpc
    Re : effective radius in kpc
    n : index
    Mtot : total mass
"""
function M_sersic(r, Re::Float64, n::Float64, Mtot::Float64)
    pn = p_ln(n)
    bn = b_cb(n)
    alpha = (3 - pn) * n
    x = @. bn * (r / Re) ^ (1.0 / n)
    return @. Mtot * gammacdf(alpha, 1.0, x)
end

struct SersicModel <: DensityModel
    Re::Float64
    n::Float64
    Mtot::Float64
end

SersicModel(Re, n) = SersicModel(Re, n, 1.0)
SersicModel() = SersicModel(10.0, 4.0)

function density2d(model::SersicModel, R)
    return s_sersic(R, model.Re, model.n, model.Mtot)
end

function density(model::SersicModel, r)
    return rho_sersic(r, model.Re, model.n, model.Mtot)
end

function mass(model::SersicModel, r)
    return M_sersic(r, model.Re, model.n, model.Mtot)
end


    
