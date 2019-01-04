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

GNFWModel() = GNFWModel(50.0, 1e7, 1.0)
density(halo::GNFWModel, r) = rho_GNFW(r, halo.rs, halo.rhos, halo.gamma)
mass(halo::GNFWModel, r) = M_GNFW(r, halo.rs, halo.rhos, halo.gamma)
scale_radius(halo::GNFWModel) = halo.rs

function GNFW_from_virial(Mvir, cvir, gamma;
                         mdef = "200c",
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

CoreNFWModel(halo::NFWModel, rc::Float64, nc::Float64) = CoreNFWModel(halo.rs, halo.rhos, rc, nc)

scale_radius(halo::CoreNFWModel) = halo.rs

density(halo::CoreNFWModel, r) = begin
    f = tanh(r / halo.rc)
    fn = f ^ halo.nc
    fnm1 = f ^ (halo.nc - 1.0)
    ρnfw = rho_NFW(r, halo.rs, halo.rhos)
    Mnfw = M_NFW(r, halo.rs, halo.rhos)
    return fn * ρnfw + halo.nc * fnm1 * (1.0 - f ^ 2) * Mnfw / (4π * r^2 * halo.rc)
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
                             mdef = "200c",
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
    
#==========================================================
Soliton NFW model, see Marsh & Pop 2015, Robles et al. 2019
https://ui.adsabs.harvard.edu/#abs/2015MNRAS.451.2479M
https://ui.adsabs.harvard.edu/#abs/arXiv:1807.06018
==========================================================#

function rho_sol(r, rsol, rhosol)
    @. rhosol * (1.0 + (r / rsol)^2) ^ -8
end

# it's analytic!
function M_sol(r, rsol, rhosol)
    x = r / rsol
    t = atan.(x, 1.0)
    term0 = 27720. * t
    term1 = 17325. * sin.(2 * t)
    term2 = -1155. * sin.(4 * t)
    term3 = -4235. * sin.(6 * t)
    term4 = -2625. * sin.(8 * t)
    term5 = -903. * sin.(10 * t)
    term6 = -175. * sin.(12 * t)
    term7 = -15. * sin.(14 * t)
    factor = 4π * rhosol * rsol^3 / 1720320.
    return @. factor * (term0 + term1 + term2 + term3 + term4 + term5 + term6 + term7)
end

function drhodr_sol(r, rsol, rhosol)
    x = r / rsol
    return @. -16 * rhosol * x * (1.0 + x^2) ^ -9
end
    
"""
Compute the matching radius, repsilon, from the specified NFW and SolNFW parameters.
"""
function matching_radius(rs, rhos, rsol, rhosol; xstart = 2.0, rtol = 1e-9)
    ρ1(r) = rho_NFW(r, rs, rhos)
    dρ1(r) = drhodr_NFW(r, rs, rhos)
    ρ2(r) = rho_sol(r, rsol, rhosol)
    dρ2(r) = drhodr_sol(r, rsol, rhosol)
    f(r) = (1.0 - ρ1(r) / ρ2(r)) ^ 2
    fp(r) = -2.0 * (1.0 - ρ1(r) / ρ2(r)) * (dρ1(r) / ρ2(r) - ρ1(r) * dρ2(r) / ρ2(r)^2)
    return fzero(f, fp, xstart * rsol; rtol = rtol)
end

"""
Soliton NFW model.  For consistency, the NFW density must match the soliton 
density at repsilon.  If that is not the case, the constructor will update
repsilon.

    rs : NFW scale radius in kpc
    rhos : NFW scale density in Msun / kpc3
    rsol : soliton scale radius in kpc
    rhosol : soliton scale density in kpc
    repsilon : transition radius
"""
struct SolNFWModel <: HaloModel
    rs::Float64
    rhos::Float64
    rsol::Float64
    rhosol::Float64
    repsilon::Float64
    # enforce density matching
    SolNFWModel(rs, rhos, rsol, rhosol, repsilon) = begin
        ρnfw = rho_NFW(repsilon, rs, rhos)
        ρsol = rho_sol(repsilon, rsol, rhosol)
        if ((ρnfw - ρsol) / ρsol) ^ 2 < 1e-9
            return new(rs, rhos, rsol, rhosol, repsilon)
        else
            @warn("recalculating matching radius", maxlog=1)
            repsilon = matching_radius(rs, rhos, rsol, rhosol)
            return new(rs, rhos, rsol, rhosol, repsilon)
        end
    end
end

density(halo::SolNFWModel, r::T where T <: Real) = begin
    if r < halo.repsilon
        return rho_sol(r, halo.rsol, halo.rhosol)
    end
    return rho_NFW(r, halo.rs, halo.rhos)
end

density(halo::SolNFWModel, r::Array{T, 1} where T <: Real) = begin
    ρ = zero(r)
    idx_sol = r .< halo.repsilon
    idx_nfw = .~ idx_sol
    ρ[idx_sol] = rho_sol(r[idx_sol], halo.rsol, halo.rhosol)
    ρ[idx_nfw] = rho_NFW(r[idx_nfw], halo.rs, halo.rhos)
    return ρ
end

mass(halo::SolNFWModel, r::T where T <: Real) = begin
    if r < halo.repsilon
        return M_sol(r, halo.rsol, halo.rhosol)
    end
    dM_epsilon = (M_sol(halo.repsilon, halo.rsol, halo.rhosol) -
                  M_NFW(halo.repsilon, halo.rs, halo.rhos))
    return M_NFW(r, halo.rs, halo.rhos) + dM_epsilon
end

mass(halo::SolNFWModel, r::Array{T, 1} where T <: Real) = begin
    M = zero(r)
    dM_epsilon = (M_sol(halo.repsilon, halo.rsol, halo.rhosol) -
                  M_NFW(halo.repsilon, halo.rs, halo.rhos))
    idx_sol = r .< halo.repsilon
    idx_nfw = .~ idx_sol
    M[idx_sol] = M_sol(r[idx_sol], halo.rsol, halo.rhosol)
    M[idx_nfw] = M_NFW(r[idx_nfw], halo.rs, halo.rhos) .+ dM_epsilon
    return M
end

scale_radius(halo::SolNFWModel) = halo.rs
