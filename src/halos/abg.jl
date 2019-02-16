using Slomo.Utils: hyp2f1

#================================================
alpha-beta-gamma model

A double power law where -γ is the inner log 
slope, -β is the outer log slope, and α is the 
sharpness of the transition.

See e.g., Merritt et al. 2006.
http://adsabs.harvard.edu/abs/2006AJ....132.2685M
================================================#

"""
Enclosed mass for αβγ model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    alpha : transition sharpness
    beta : negative outer log slope
    gamma : negative inner log slope
"""
function M_ABG(r, rs, rhos, alpha, beta, gamma)
    x = r / rs
    y = 4π * rhos * rs^3
    omega = 3.0 - gamma
    hyp = hyp2f1.(omega / alpha,
                 (beta - gamma) / alpha,
                 1 + omega / alpha,
                 -(x .^ alpha))
    return @. y * x ^ omega * hyp / omega
end

"""
Local volume density for αβγ model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    alpha : transition sharpness
    beta : negative outer log slope
    gamma : negative inner log slope
"""
function rho_ABG(r, rs, rhos, alpha, beta, gamma)
    x = r / rs
    return @. rhos * (x ^ -gamma * (1 + x ^ alpha) ^ ((-beta - gamma) / alpha))
end

"""
Linear density derivative for αβγ model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    alpha : transition sharpness
    beta : negative outer log slope
    gamma : negative inner log slope
"""
function drhodr_ABG(r, rs, rhos, alpha, beta, gamma)
    x = r / rs
    α = alpha
    β = beta
    γ = gamma
    term1 = @. -γ * x ^ (-γ - 1) * (1 + x ^ α)^(-(β - γ) / α)
    term2 = @. -(β - γ) * x ^ (α - γ - 1) * (1 + x ^ α)^(-(β - γ) / α - 1)
    return @. rhos / rs * (term1 + term2)
end

"""
ABG halo density model.

    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    alpha : transition sharpness
    beta : negative outer log slope
    gamma : negative inner log slope
"""
struct ABGModel <: HaloModel
    rs::Float64
    rhos::Float64
    alpha::Float64
    beta::Float64
    gamma::Float64
end

# default to 1e12 Mvir / Msun NFW halo
ABGModel() = ABGModel(25.3, 3.7e6, 1, 3, 1)
density(halo::ABGModel, r) = rho_ABG(r, halo.rs, halo.rhos,
                                     halo.alpha, halo.beta, halo.gamma)
mass(halo::ABGModel, r) = M_ABG(r, halo.rs, halo.rhos,
                                halo.alpha, halo.beta, halo.gamma)
scale_radius(halo::ABGModel) = begin
    ((2 - halo.gamma) / (halo.beta - 2)) ^ (1 / halo.alpha) * halo.rs
end


function ABG_from_virial(Mvir, cvir, alpha, beta, gamma;
                         mdef = default_mdef,
                         cosmo = default_cosmo,
                         z = 0.0)
    Rvir = Rvir_from_Mvir(Mvir; mdef = mdef, cosmo = cosmo, z = z)
    r2 = Rvir / cvir
    rs = ((2 - gamma) / (beta - 2))^(-1 / alpha) * r2
    rhos = Mvir / M_ABG(Rvir, rs, 1.0, alpha, beta, gamma)
    return ABGModel(rs, rhos, alpha, beta, gamma)
end
