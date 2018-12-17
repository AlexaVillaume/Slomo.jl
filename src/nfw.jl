"""
Enclosed mass for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function M_NFW(r, rs::Float64, rhos::Float64)
    x = r / rs
    y = 4π * rhos * rs^3
    return @. y * (log(1.0 + x) - x / (1.0 + x))
end

"""
Local volume density for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function rho_NFW(r, rs::Float64, rhos::Float64)
    x = r / rs
    return @. rhos * (x ^ -1 * (1 + x) ^ -2)
end

struct NFWModel <: HaloModel
    rs::Float64
    rhos::Float64
end

NFWModel() = NFWModel(50.0, 1e7)
density(halo::NFWModel, r) = rho_NFW(r, halo.rs, halo.rhos)
mass(halo::NFWModel, r) = M_NFW(r, halo.rs, halo.rhos)
scale_radius(halo::NFWModel) = halo.rs

function NFW_from_virial(Mvir, cvir;
                         mdef = "200c",
                         cosmo = default_cosmo,
                         z = 0.0)
    Rvir = Rvir_from_Mvir(Mvir; mdef = mdef, cosmo = cosmo, z = z)
    rs = Rvir / cvir
    rhos = Mvir / M_NFW(Rvir, rs, 1.0)
    return NFWModel(rs, rhos)
end
