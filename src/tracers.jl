""" Collection of tracer density models.

[`SersicModel`](@ref) : Representation of a Sersic surface density profile
"""
module Tracers

import Slomo.Models: DensityModel, mass, density, density2d
import Slomo.Utils: gamma, gamma_inc

export SersicModel

"""
    b_cb(n)

Compute the "b" parameter in the Sersic function, from the Ciotti & Bertin (1999) 
approximation.  `n` is the Sersic index.
"""
function b_cb(n)
    -1. / 3 + 2. * n + 4 / (405. * n) + 46 / (25515. * n^2)
end

"""
    p_ln(n)

Compute the "p" parameter in fitting deprojected Sersic function, from 
Lima Neto et al. (1999).  `n` is the Sersic index
"""
function p_ln(n)
    1. - 0.6097 / n + 0.05463 / n^2
end

"""
    s_sersic(R, Re, n, Mtot)

Compute the surface density profile for Sersic model.

* `R`: projected (2d) radii in kpc
* `Re`: effective radius in kpc
* `n`: index
* `Mtot`: total mass
"""
function s_sersic(R, Re, n, Mtot)
    bn = b_cb(n)
    s0 = Mtot * bn ^ (2 * n) / (2π * n * Re ^ 2 * gamma(2 * n))
    return s0 * @. exp(-bn * (R / Re) ^ (1.0 / n))
end

"""
    rho_sersic(r, Re, n, Mtot)

Compute the volume density profile for Sersic model.

* `r`: deprojected (3d) radii in kpc
* `Re`: effective radius in kpc
* `n`: index
* `Mtot`: total mass
"""
function rho_sersic(r, Re, n, Mtot)
    bn = b_cb(n)
    pn = p_ln(n)
    x = r / Re
    rho0 = Mtot * bn ^ (3 * n) / (4π * n * Re ^ 3 * gamma((3 - pn) * n))
    return rho0 * (@. (bn ^ n * x) ^ (-pn) * exp(-bn * x ^ (1.0 / n)))
end

"""
    M_sersic(r, Re, n, Mtot)

Compute the enclosed mass within deprojected radius for the Sersic model.

* `r`: deprojected (3d) radii in kpc
* `Re`: effective radius in kpc
* `n`: index
* `Mtot`: total mass
"""
function M_sersic(r, Re, n, Mtot)
    pn = p_ln(n)
    bn = b_cb(n)
    alpha = (3 - pn) * n
    x = @. bn * (r / Re) ^ (1.0 / n)
    return @. Mtot * gamma_inc(alpha, x)
end

"""
    SersicModel(Re, n, Mtot)

Representation of a Sersic surface density profile.

* `Re`: the effective radius in kpc.
* `n`: the Sersic index
* `Mtot`: the total mass associated to the profile.

For the purposes of using the `SersicModel` to track the tracer population's
density distribution, `Mtot` is somewhat arbitrary.  For using it as a component
of the total mass profile, this should take units of solar masses.
"""
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

end
