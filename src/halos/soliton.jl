"""
    rho_sol(r, rsol, rhosol)

Compute the local volume density for the soliton model.

* `rsol`: soliton scale radius in kpc
* `rhosol`: soliton scale density in kpc
"""
function rho_sol(r, rsol, rhosol)
    @. rhosol * (1.0 + (r / rsol)^2) ^ -8
end

"""
    M_sol(r, rsol, rhosol)

Compute the enclosed mass for the soliton model.  It's analytic!

* `rsol`: soliton scale radius in kpc
* `rhosol`: soliton scale density in kpc
"""
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

"""
    drhodr_sol(r, rsol, rhosol)

Compute the linear density derivative for the soliton model.

* `rsol`: soliton scale radius in kpc
* `rhosol`: soliton scale density in kpc
"""
function drhodr_sol(r, rsol, rhosol)
    x = r / rsol
    return @. -16 * rhosol * x * (1.0 + x^2) ^ -9
end

"""
    m22_from_sol(rsol, rhosol; <keyword arguments>)

Calculate the axion mass from the soliton scale parameters. See Marsh & Pop 2015, 
equation 8.  Axion mass has units of 1e-22 eV.  

* `rsol`: soliton scale radius in kpc
* `rhosol`: soliton scale density in kpc

# Arguments
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `alpha_mp::Real = 0.23`: the alpha fitting parameter from Marsh & Pop 2015
"""
function m22_from_sol(rsol, rhosol;
                      cosmo = default_cosmo, z = 0.0, alpha_mp = 0.23)
    Δ_sol = rhosol / ρcrit(z; cosmo = cosmo)
    m22 = √(Δ_sol ^ -1 * rsol ^ -4 * (cosmo.h / 0.7) ^ 2 * 5e4 * alpha_mp ^ -4)
    return m22
end

"""
    rhosol_from_rsol(m22, rsol; <keyword arguments>)

Calculate the scale density of the soliton core from m22 and the scale radius.

* `m22`: ultralight axion mass in units of 1e-22 eV
* `rsol`: soliton scale radius in kpc

# Arguments
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `alpha_mp::Real = 0.23`: the alpha fitting parameter from Marsh & Pop 2015
"""
function rhosol_from_rsol(m22, rsol;
                         cosmo = default_cosmo, z = 0.0, alpha_mp = 0.23)
    Δ_sol = (5e4 / alpha_mp^4) * (cosmo.h / 0.7)^-2 * m22^-2 * rsol^-4
    return Δ_sol * ρcrit(z; cosmo = cosmo)
end

"""
    rsol_from_Mvir(m22, Mvir)

Calculate the scale radius of the soliton core from the halo mass scaling relation 
(see Robles+2019).
"""
function rsol_from_Mvir(m22, Mvir)
    rsol = 3.315 * 1.6 * (Mvir / 1e9)^(-1/3.) * m22^-1
end

"""
    matching_radius(rs, rhos, rsol, rhosol; <keyword arguments>)

Compute the matching radius, repsilon, from the specified NFW and Sol parameters.

# Arguments
- `xstart::Real = 2.0`: factor of the soliton radius to start the search for the root
- `rtol::Real = 1e-9`: relative error tolerance for the root finding
- `maxevals:Int = 100`: maximum number of function evalutions for the root finding
"""
function matching_radius(rs, rhos, rsol, rhosol;
                         xstart = 2.0, rtol = 1e-9, maxevals = 100)
    ρ1(r) = rho_NFW(r, rs, rhos)
    dρ1(r) = drhodr_NFW(r, rs, rhos)
    ρ2(r) = rho_sol(r, rsol, rhosol)
    dρ2(r) = drhodr_sol(r, rsol, rhosol)
    f(r) = (1.0 - ρ1(r) / ρ2(r)) ^ 2
    fp(r) = -2.0 * (1.0 - ρ1(r) / ρ2(r)) * (dρ1(r) / ρ2(r) - ρ1(r) * dρ2(r) / ρ2(r)^2)
    return fzero(f, fp, xstart * rsol; rtol = rtol, maxevals = maxevals)
end

"""
    SolNFWModel(rs, rhos, rsol, rhosol, repsilon)

Soliton NFW model.  For consistency, the NFW density must match the soliton 
density at repsilon.  If that is not the case, the constructor will update
repsilon.

- [Marsh & Pop 2015](https://ui.adsabs.harvard.edu/#abs/2015MNRAS.451.2479M)
- [Robles, Bullock, and Boylan-Kolchin 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483..289R/)

* `rs`: NFW scale radius in kpc
* `rhos`: NFW scale density in Msun / kpc3
* `rsol`: soliton scale radius in kpc
* `rhosol`: soliton scale density in kpc
* `repsilon`: transition radius
"""
struct SolNFWModel <: HaloModel
    rs::Float64
    rhos::Float64
    rsol::Float64
    rhosol::Float64
    repsilon::Float64
    # enforce density matching
    SolNFWModel(rs, rhos, rsol, rhosol, repsilon;
                rtol = 1e-9, maxevals = 100) = begin
        ρnfw = rho_NFW(repsilon, rs, rhos)
        ρsol = rho_sol(repsilon, rsol, rhosol)
        if ((ρnfw - ρsol) / ρsol) ^ 2 < 1e-9
            return new(rs, rhos, rsol, rhosol, repsilon)
        else
            @warn("recalculating matching radius", maxlog=1)
            try
                repsilon = matching_radius(rs, rhos, rsol, rhosol;
                                           rtol = rtol, maxevals = maxevals)
            catch err
                @warn("failed to find matching radius for " *
                      "\trs = $rs \n\trhos = $rhos \n\trsol = $rsol \n\trhosol = $rhosol")
                throw(err)
            end
            return new(rs, rhos, rsol, rhosol, repsilon)
        end
    end
end

SolNFWModel() = SolNFWModel(21.1, 5.6e6, 0.53, 3.1e10, 0.48)

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

"""
    SolNFW_from_virial(Mvir, cvir, m22; <keyword arguments>)

Construct a Solition-NFW halo from the halo mass, concentration, and axion mass.

# Arguments
- `rsol = nothing`: for default, use the soliton size-halo mass scaling relationship
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `xstart::Real = 10.0`: factor of the scale radius to start the search for the root
- `rtol::Real = 1e-2`: relative error tolerance for the root finding
- `maxevals:Int = 100`: maximum number of function evalutions for the root finding
"""
function SolNFW_from_virial(Mvir, cvir, m22;
                            rsol = nothing,
                            mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                            rtol = 1e-9, maxevals=100)
    # soliton parameters, calculate rsol from Mvir scaling if not passed in
    if rsol == nothing
        rsol = rsol_from_Mvir(m22, Mvir)
    end
    rhosol = rhosol_from_rsol(m22, rsol, cosmo = cosmo, z = z)
    # NFW parameters
    Rvir = Rvir_from_Mvir(Mvir; mdef = mdef, cosmo = cosmo, z = z)
    rs = Rvir / cvir

    # first guess, rhos from normal NFW profile
    rhos = Mvir / M_NFW(Rvir, rs, 1.0)
    repsilon = matching_radius(rs, rhos, rsol, rhosol;
                               rtol = rtol, maxevals = maxevals)

    return SolNFWModel(rs, rhos, rsol, rhosol, repsilon)
end

function matching_radius(rs, rhos, alpha, beta, gamma, rsol, rhosol;
                         xstart = 2.0, rtol = 1e-6, maxevals = 100)
    ρ1(r) = rho_ABG(r, rs, rhos, alpha, beta, gamma)
    dρ1(r) = drhodr_ABG(r, rs, rhos, alpha, beta, gamma)
    ρ2(r) = rho_sol(r, rsol, rhosol)
    dρ2(r) = drhodr_sol(r, rsol, rhosol)
    f(r) = (1.0 - ρ1(r) / ρ2(r)) ^ 2
    fp(r) = -2.0 * (1.0 - ρ1(r) / ρ2(r)) * (dρ1(r) / ρ2(r) - ρ1(r) * dρ2(r) / ρ2(r)^2)
    return fzero(f, fp, xstart * rsol; rtol = rtol, maxevals = 100)
end

"""
    SolABGModel(rs, rhos, alpha, beta, gamma, rsol, rhosol, repsilon)

Soliton + αβγ model.  For consistency, the αβγ density profile must match the 
soliton density profile at repsilon.  If that is not the case, the constructor 
will update repsilon.

* `rs`: NFW scale radius in kpc
* `rhos`: NFW scale density in Msun / kpc3
* `alpha`: transition sharpness
* `beta`: negative outer log slope
* `gamma`: negative inner log slope
* `rsol`: soliton scale radius in kpc
* `rhosol`: soliton scale density in kpc
* `repsilon`: transition radius
"""
struct SolABGModel <: HaloModel
    rs::Float64
    rhos::Float64
    alpha::Float64
    beta::Float64
    gamma::Float64
    rsol::Float64
    rhosol::Float64
    repsilon::Float64
    # enforce density matching
    SolABGModel(rs, rhos, alpha, beta, gamma, rsol, rhosol, repsilon;
                rtol = 1e-6, maxevals = 100) = begin
        ρabg = rho_ABG(repsilon, rs, rhos, alpha, beta, gamma)
        ρsol = rho_sol(repsilon, rsol, rhosol)
        if ((ρabg - ρsol) / ρsol) ^ 2 < rtol
            return new(rs, rhos, alpha, beta, gamma, rsol, rhosol, repsilon)
        else
            @warn("recalculating matching radius", maxlog=1)
            try
                repsilon = matching_radius(rs, rhos, alpha, beta, gamma,
                                           rsol, rhosol;
                                           rtol = rtol, maxevals = maxevals)
            catch err
                @warn("failed to find matching radius for " *
                      "\trs = $rs \n\trhos = $rhos \n\trsol = $rsol \n\trhosol = $rhosol")
                throw(err)
            end
            return new(rs, rhos, alpha, beta, gamma, rsol, rhosol, repsilon)
        end
    end
end

SolABGModel() = SolABGModel(21.1, 5.6e6, 1, 3, 1, 0.53, 3.1e10, 0.48)

density(halo::SolABGModel, r::T where T <: Real) = begin
    if r < halo.repsilon
        return rho_sol(r, halo.rsol, halo.rhosol)
    end
    return rho_ABG(r, halo.rs, halo.rhos, halo.alpha, halo.beta, halo.gamma)
end

density(halo::SolABGModel, r::Array{T, 1} where T <: Real) = begin
    ρ = zero(r)
    idx_sol = r .< halo.repsilon
    idx_abg = .~ idx_sol
    ρ[idx_sol] = rho_sol(r[idx_sol], halo.rsol, halo.rhosol)
    ρ[idx_abg] = rho_ABG(r[idx_abg], halo.rs, halo.rhos, halo.alpha, halo.beta, halo.gamma)
    return ρ
end

mass(halo::SolABGModel, r::T where T <: Real) = begin
    if r < halo.repsilon
        return M_sol(r, halo.rsol, halo.rhosol)
    end
    dM_epsilon = (M_sol(halo.repsilon, halo.rsol, halo.rhosol) -
                  M_ABG(halo.repsilon, halo.rs, halo.rhos, halo.alpha, halo.beta, halo.gamma))
    return M_ABG(r, halo.rs, halo.rhos, halo.alpha, halo.beta, halo.gamma) + dM_epsilon
end

mass(halo::SolABGModel, r::Array{T, 1} where T <: Real) = begin
    M = zero(r)
    dM_epsilon = (M_sol(halo.repsilon, halo.rsol, halo.rhosol) -
                  M_ABG(halo.repsilon, halo.rs, halo.rhos, halo.alpha, halo.beta, halo.gamma))
    idx_sol = r .< halo.repsilon
    idx_nfw = .~ idx_sol
    M[idx_sol] = M_sol(r[idx_sol], halo.rsol, halo.rhosol)
    M[idx_nfw] = M_ABG(r[idx_nfw], halo.rs, halo.rhos, halo.alpha, halo.beta, halo.gamma) .+ dM_epsilon
    return M
end

scale_radius(halo::SolABGModel) = begin
    ((2 - halo.gamma) / (halo.beta - 2)) ^ (1 / halo.alpha) * halo.rs
end

"""
    SolABG_from_virial(Mvir, cvir, alpha, beta, gamma, m22; <keyword arguments>)

Construct a Soliton + αβγ halo from the halo mass, concentration, shape parameters, and 
axion mass.

# Arguments
- `rsol = nothing`: for default, use the soliton size-halo mass scaling relationship
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `xstart::Real = 10.0`: factor of the scale radius to start the search for the root
- `rtol::Real = 1e-2`: relative error tolerance for the root finding
- `maxevals:Int = 100`: maximum number of function evalutions for the root finding
"""
function SolABG_from_virial(Mvir, cvir, alpha, beta, gamma, m22;
                            rsol = nothing,
                            mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                            rtol = 1e-9, maxevals=100)
    # soliton parameters, calculate rsol from Mvir scaling if not passed in
    if rsol == nothing
        rsol = rsol_from_Mvir(m22, Mvir)
    end
    rhosol = rhosol_from_rsol(m22, rsol, cosmo = cosmo, z = z)
    # ABG parameters
    Rvir = Rvir_from_Mvir(Mvir; mdef = mdef, cosmo = cosmo, z = z)
    r2 = Rvir / cvir
    rs = ((2 - gamma) / (beta - 2))^(-1 / alpha) * r2

    # first guess, rhos from normal ABG profile
    rhos = Mvir / M_ABG(Rvir, rs, 1.0, alpha, beta, gamma)
    repsilon = matching_radius(rs, rhos, alpha, beta, gamma, rsol, rhosol;
                               rtol = rtol, maxevals = maxevals)

    return SolABGModel(rs, rhos, alpha, beta, gamma, rsol, rhosol, repsilon)
end
