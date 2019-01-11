using Slomo.Utils: gamma, gamma_inc

#================================================
Einasto model, Einasto 1965
http://adsabs.harvard.edu/abs/1965TrAlm...5...87E
================================================#

"""
Enclosed mass for Einasto model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    alpha : Einasto shape index, lower is more cored
"""
function M_Einasto(r, rs, rhos, alpha)
    x = r / rs
    M = (4π * rhos * rs *
         exp(2.0 / alpha) * gamma(3.0 / alpha) * (2.0 / alpha) ^ (-3.0 / alpha))
    return @. M * gamma_inc(3.0 / alpha, 2.0 / alpha * x ^ alpha)
end

"""
Local volume density for Einasto model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    alpha : Einasto shape index, lower is more cored
"""
function rho_Einasto(r, rs, rhos, alpha)
    x = r / rs
    return @. rhos * exp(- 2.0 / alpha * (x^alpha - 1.0))
end

"""
Einasto halo model

    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    alpha : Einasto shape index, lower is more cored
"""
struct EinastoModel <: HaloModel
    rs::Float64
    rhos::Float64
    alpha::Float64
end

# default to 1e12 Mvir / Msun halo
EinastoModel() = EinastoModel(25.3, 8.4e5, 0.16)
density(halo::EinastoModel, r) = rho_Einasto(r, halo.rs, halo.rhos, halo.alpha)
mass(halo::EinastoModel, r) = M_Einasto(r, halo.rs, halo.rhos, halo.alpha)
scale_radius(halo::EinastoModel) = halo.rs

function Einasto_from_virial(Mvir, cvir, alpha;
                             mdef = default_mdef,
                             cosmo = default_cosmo,
                             z = 0.0)
    Rvir = Rvir_from_Mvir(Mvir; mdef = mdef, cosmo = cosmo, z = z)
    rs = Rvir / cvir
    rhos = Mvir / M_Einasto(Rvir, rs, 1.0, alpha)
    return EinastoModel(rs, rhos, alpha)
end
