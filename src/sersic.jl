using GSL: sf_gamma

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

    R : radii in kpc
    Re : effective radius in kpc
    n : index
"""
function s_sersic(R, Re::Float64, n::Float64)
    bn = b_cb(n)
    s0 = bn ^ (2 * n) / (2pi * n * Re ^ 2 * sf_gamma(2 * n))
    return s0 * @. exp(-bn * (R / Re) ^ (1.0 / n))
end

"""
Volume density for Sersic model.

    r : radii in kpc
    Re : effective radius in kpc
    n : index
"""
function rho_sersic(r, Re::Float64, n::Float64)
    bn = b_cb(n)
    pn = p_ln(n)
    x = r / Re
    rho0 = bn ^ (3 * n) / (4pi * n * Re ^ 3 * sf_gamma((3 - pn) * n))
    return rho0 * (@. (bn ^ n * x) ^ (-pn) * exp(-bn * x ^ (1.0 / n)))
end

"""
Mass enclosed within deprojected radius for the Sersic model.

    r : radii in kpc
    Re : effective radius in kpc
    n : index
    M : total mass
"""
function M_sersic(r, Re, n::Float64, Mtot::Float64)
    pn = p_ln(n)
    bn = b_cb(n)
    return Mtot * @. sf_gamma_inc((3 - pn) * n, bn * (r / Re) ^ (1.0 / n))
end

struct SersicModel <: DensityModel
    Re::Float64
    n::Float64
    Mtot::Float64
end

SersicModel(Re, n) = SersicModel(Re, n, 1.0)
SersicModel() = SersicModel(10.0, 4.0)

function density2d(model::SersicModel, R)
    return s_sersic(R, model.Re, model.n)
end

function density(model::SersicModel, r)
    return rho_sersic(r, model.Re, model.n)
end

function mass(model::SersicModel, r)
    return M_sersic(r, model.Re, model.n, model.Mtot)
end


    
