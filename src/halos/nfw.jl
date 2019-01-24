using Slomo.Utils: hyp2f1

#================================================
NFW model, Navarro, Frenk, & White 1997
http://adsabs.harvard.edu/abs/1997ApJ...490..493N
================================================#

"""
Enclosed mass for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function M_NFW(r, rs, rhos)
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
function rho_NFW(r, rs, rhos)
    x = r / rs
    return @. rhos * (x ^ -1 * (1 + x) ^ -2)
end

"""
Linear density derivative for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
function drhodr_NFW(r, rs, rhos)
    x = r / rs
    return @. -rhos / rs * ((1.0 + x)^-2 * x^-2 + 2.0 * (1.0 + x)^-3 * x^-1)
end

"""
NFW halo density model.

    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
"""
struct NFWModel <: HaloModel
    rs::Float64
    rhos::Float64
end

# default to 1e12 Mvir / Msun halo
NFWModel() = NFWModel(25.3, 3.7e6)
density(halo::NFWModel, r) = rho_NFW(r, halo.rs, halo.rhos)
mass(halo::NFWModel, r) = M_NFW(r, halo.rs, halo.rhos)
scale_radius(halo::NFWModel) = halo.rs

function NFW_from_virial(Mvir, cvir;
                         mdef = default_mdef,
                         cosmo = default_cosmo,
                         z = 0.0)
    Rvir = Rvir_from_Mvir(Mvir; mdef = mdef, cosmo = cosmo, z = z)
    rs = Rvir / cvir
    rhos = Mvir / M_NFW(Rvir, rs, 1.0)
    return NFWModel(rs, rhos)
end

#================================================
Generalized NFW model (gNFW)

Instance of a generalized Hernquist model with
the transition parameter α = 1, the outer slope 
β = 3, and the inner slope, γ, left as a free 
parameter.  See Hernquist 1990 and Zhao 1996.

http://adsabs.harvard.edu/abs/1990ApJ...356..359H
http://adsabs.harvard.edu/abs/1996MNRAS.278..488Z
================================================#

"""
Enclosed mass for gNFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    gamma : negative of inner density log slope
"""
function M_GNFW(r, rs, rhos, gamma)
    omega = 3.0 - gamma
    x = r / rs
    y = 4π * rhos * rs^3 / omega
    return @. y * x ^ omega * hyp2f1(omega, omega, omega + 1.0, -x)
end

"""
Local volume density for NFW model.

    r : radii in kpc
    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    gamma : negative of inner density log slope
"""
function rho_GNFW(r, rs, rhos, gamma)
    x = r / rs
    return @. rhos * (x ^ -gamma) * (1 + x) ^ (gamma - 3.0)
end

"""
gNFW halo density model.

    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    gamma : negative of inner density log slope    
"""
struct GNFWModel <: HaloModel
    rs::Float64
    rhos::Float64
    gamma::Float64
end

# default to 1e12 Mvir cusped halo
GNFWModel() = GNFWModel(25.3, 3.7e6, 1.0)
density(halo::GNFWModel, r) = rho_GNFW(r, halo.rs, halo.rhos, halo.gamma)
mass(halo::GNFWModel, r) = M_GNFW(r, halo.rs, halo.rhos, halo.gamma)
scale_radius(halo::GNFWModel) = halo.rs

function GNFW_from_virial(Mvir, cvir, gamma;
                         mdef = default_mdef,
                         cosmo = default_cosmo,
                         z = 0.0)
    Rvir = Rvir_from_Mvir(Mvir; mdef = mdef, cosmo = cosmo, z = z)
    rs = Rvir / cvir
    rhos = Mvir / M_GNFW(Rvir, rs, 1.0, gamma)
    return GNFWModel(rs, rhos, gamma)
end

#================================================
coreNFW model, Read, Agertz, & Collins 2016
http://adsabs.harvard.edu/abs/2016MNRAS.459.2573R
================================================#

"""
coreNFW halo density model

    rs : scale radius in kpc
    rhos : scale density in Msun / kpc3
    rc : core radius in kpc
    nc : core index (0 -> no core, 1 -> full core)
"""
struct CoreNFWModel <: HaloModel
    rs::Float64
    rhos::Float64
    rc::Float64
    nc::Float64
end

CoreNFWModel() = CoreNFWModel(25.3, 3.7e6, 12.5, 0.5)

scale_radius(halo::CoreNFWModel) = halo.rs

density(halo::CoreNFWModel, r) = begin
    f = tanh.(r / halo.rc)
    fn = f .^ halo.nc
    fnm1 = f .^ (halo.nc - 1.0)
    ρnfw = rho_NFW.(r, halo.rs, halo.rhos)
    Mnfw = M_NFW.(r, halo.rs, halo.rhos)
    return @. fn * ρnfw + halo.nc * fnm1 * (1.0 - f ^ 2) * Mnfw / (4π * r^2 * halo.rc)
end

mass(halo::CoreNFWModel, r) = begin
    @. M_NFW(r, halo.rs, halo.rhos) * tanh(r / halo.rc) ^ halo.nc
end

"""
Construct a CoreNFWModel from the scaling relations in Read et al. 2016.

    Mvir : virial mass in Msun
    cvir : halo concentration
    Re : effective radius of stellar distribution, in kpc
    t_sf : time since start of star formation, in Gyr
"""
function CoreNFW_from_virial(Mvir, cvir, Re, t_sf;
                             mdef = default_mdef,
                             cosmo = default_cosmo,
                             z = 0.0,
                             κ = 0.04,
                             η = 1.75)
    nfw_halo = NFW_from_virial(Mvir, cvir; mdef = mdef, cosmo = cosmo, z = z)
    G_alt = 4.49850215e-06 # G in Msun^-1 kpc^3 Gyr^-2
    rs = nfw_halo.rs
    rhos = nfw_halo.rhos
    rc = κ * Re * 4.0 / 3.0
    t_dyn = 2π * sqrt(rs ^ 3 / (G_alt * M_NFW(rs, rs, rhos)))
    nc = tanh(κ * t_sf / t_dyn)
    return CoreNFWModel(rs, rhos, rc, nc)
end    
    
