"""
Dark matter halo models and relations.

Significant credit to Benedikt Diemer for developing [Colossus](https://bitbucket.org/bdiemer/colossus).

While much of this is inspired by Colossus; any implementation faults are my own.

The default halo mass definition is "200c", i.e. the virial radius of the halo is defined
to be the radius that encloses a density that is 200 times the critical density of the
universe.

The default cosmology (from `Slomo.CosmologyTools`) is currently set to the values
from Planck 2018 (h = 0.6766, Neff = 3.046, OmegaM = 0.3111).

# Common HaloModel methods
* [`scale_radius`](@ref): scale radius of the halo
* [`virial_radius`](@ref): virial radius of the halo
* [`virial_mass`](@ref): virial mass of the halo
* [`concentration`](@ref): halo concentration (defined as rvir / rs)

# HaloModel subtypes
* [`NFWModel`](@ref): Navarro, Frenk, and White (1997) profile
* [`CoreNFWModel`](@ref): Core-NFW (Read et al. 2016) profile
* [`GNFWModel`](@ref): generalized NFW model (i.e., αβγ where α = 1, β = 3)
* [`ABGModel`](@ref): α-β-γ profile (double power law)
* [`EinastoModel`](@ref): Einasto (1965) profile
* [`SolNFWModel`](@ref): Soliton + NFW model (Marsh & Pop 2015)
* [`SolABGModel`](@ref): Soliton + αβγ model (Wasserman et al. 2019)

# Halo relations
* [`hmcr`](@ref): halo mass-concentration relation from Dutton & Maccio (2014)
* [`shmr`](@ref): stellar-halo mass relation from Rodriguez-Puebla et al. (2017)
* [`abg_from_logshm`](@ref): αβγ scaling with log(Mstar / Mhalo) from DiCintio et al. (2014)
"""
module Halos

using Roots
using ForwardDiff

using Slomo.Utils: log_gauss
using Slomo.CosmologyTools: Ωm, ρm, ρcrit, default_cosmo
import Slomo.Models: DensityModel, mass, density, NotImplemented
using Slomo.Utils: hyp2f1, gamma, gamma_inc

export virial_radius, virial_mass, scale_radius, concentration
export NFWModel, CoreNFWModel, ABGModel, GNFWModel, EinastoModel, SolNFWModel, SolABGModel
export NFW_from_virial, CoreNFW_from_virial, ABG_from_virial, GNFW_from_virial, Einasto_from_virial
export SolNFW_from_virial, SolABG_from_virial
export hmcr, hmcr_prior, shmr, shmr_prior, abg_from_logshm
    
const default_mdef = "200c"

abstract type HaloModel <: DensityModel end

"""
    scale_radius(halo::HaloModel)

Return the scale radius of the halo.  The convention adopted here is that the scale radius
occurs where the logarithmic slope of the halo density profile is equal to -2.  For the 
standard NFW model this is equal to `rs`, but other halo models may compute some other
value to stay consistent with this convention (e.g., the `ABGModel`).
"""
function scale_radius(halo::HaloModel)::Float64
    throw(NotImplemented("needs to be implemented for subtypes of HaloModel"))
end

"""
    Δρ_from_mdef(mdef; <keyword arguments>)
    
Parse the halo mass definition to get the spherical overdensity in Msun / kpc3.
If mdef is "vir", then use the spherical overdensity relation of Bryan & Norman 1998.

`mdef` should be either "vir" or of the form "[0-9]*[c|m]" where "c" refers to a multiple
of the critical density and "m" refers to a multiple of the mass density, e.g.,

    mdef in ["vir", "200c", "500c", "200m", "500m", ...]

# Arguments
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function Δρ_from_mdef(mdef; cosmo = default_cosmo, z = 0.0)
    mdef = lowercase(strip(mdef))
    if mdef == "vir"
        x = Ωm(z; cosmo = cosmo) .- 1.0
        Δ = @. 18π^2 + 82.0 * x - 39.0 * x^2
        ρ = ρcrit(z; cosmo = cosmo)
    else
        mdef_type = mdef[end]
        Δ = parse(Float64, mdef[1:end - 1])
        if mdef_type == 'm'
            ρ = ρm(z; cosmo = cosmo)
        elseif mdef_type == 'c'
            ρ = ρcrit(z; cosmo = cosmo)
        else
            throw("$mdef not recognized as a valid halo mass definition")
        end
    end
    return @. Δ * ρ
end

"""
    Rvir_from_Mvir(Mvir; <keyword arguments>)

Compute the virial radius from the virial mass, using the definition

```math
M_\\mathrm{vir} = 4\\pi \\Delta\\rho / (3 R_\\mathrm{vir}^3)
```

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function Rvir_from_Mvir(Mvir; mdef = default_mdef, cosmo = default_cosmo, z = 0.0)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    return @. (Mvir * 3.0 / 4π / Δρ) ^ (1 / 3)
end

"""
    Mvir_from_Rvir(Rvir; <keyword arguments>)

Compute the virial mass from the virial radius, using the definition

```math
M_\\mathrm{vir} = 4\\pi \\Delta\\rho / (3 R_\\mathrm{vir}^3)
```

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function Mvir_from_Rvir(Rvir; mdef = default_mdef, cosmo = default_cosmo, z = 0.0)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    return @. 4π / 3.0 * Rvir ^ 3 * Δρ
end

"""
    virial_radius(halo::HaloModel, <keyword arguments>)

Compute the virial radius of the halo via root finding.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `xstart::Real = 10.0`: factor of the scale radius to start the search for the root
- `rtol::Real = 1e-2`: relative error tolerance for the root finding
- `maxevals:Int = 100`: maximum number of function evalutions for the root finding
"""
function virial_radius(halo::HaloModel;
                       mdef = default_mdef,
                       cosmo = default_cosmo,
                       z = 0.0,
                       xstart = 10.0,
                       rtol = 1e-2,
                       maxevals = 100)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    # function to find roots
    f(r) = 3.0 / 4π * mass(halo, r) / r^3 - Δρ
    # derivative
    fp(r) = 3.0 * (r^-1 * density(halo, r) - r^-4 * mass(halo, r))
    rstart = xstart * scale_radius(halo)
    return fzero(f, fp, rstart; rtol = rtol, maxevals = maxevals)
end

"""
    virial_mass(halo::HaloModel; <keyword arguments>)

Compute the virial mass via root finding.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `xstart::Real = 10.0`: factor of the scale radius to start the search for the root
- `rtol::Real = 1e-2`: relative error tolerance for the root finding
- `maxevals:Int = 100`: maximum number of function evalutions for the root finding
"""
function virial_mass(halo::HaloModel;
                     mdef = default_mdef,
                     cosmo = default_cosmo,
                     z = 0.0,
                     xstart = 10.0,
                     rtol = 1e-2,
                     maxevals = 100)
    Rvir = virial_radius(halo;
                         mdef = mdef, cosmo = cosmo, z = z,
                         xstart = xstart, rtol = rtol, maxevals = maxevals)
    return Mvir_from_Rvir(Rvir, mdef = mdef, cosmo = cosmo, z = z)
end

"""
    virial_mass(Rvir::Real; <keyword arguments>)

Compute the virial mass from the virial radius.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
"""
function virial_mass(Rvir::Real;
                     mdef = default_mdef,
                     cosmo = default_cosmo,
                     z = 0.0)
    return Mvir_from_Rvir(Rvir, mdef = mdef, cosmo = cosmo, z = z)
end

"""
    concentration(halo::HaloModel; <keyword arguments>)

Compute the halo concentration.

# Arguments
- `mdef::AbstractString = default_mdef`: halo mass definition (e.g., "200c", "vir")
- `cosmo::AbstractCosmology = default_cosmo`: cosmology under which to evaluate the overdensity
- `z::Real = 0.0`: redshift at which to evaluate the overdensity
- `xstart::Real = 10.0`: factor of the scale radius to start the search for the root
- `rtol::Real = 1e-2`: relative error tolerance for the root finding
- `maxevals:Int = 100`: maximum number of function evalutions for the root finding
"""
function concentration(halo::HaloModel;
                       mdef = default_mdef,
                       cosmo = default_cosmo,
                       z = 0.0,
                       xstart = 10.0,
                       rtol = 1e-2,
                       maxevals = 100)
    Rvir = virial_radius(halo;
                         mdef = mdef, cosmo = cosmo, z = z,
                         xstart = xstart, rtol = rtol, maxevals = maxevals)
    rs = scale_radius(halo)
    return Rvir / rs
end

include("nfw.jl")
include("relations.jl")
include("abg.jl")
include("einasto.jl")
include("soliton.jl")

end
