"""
Dark matter halo models and relations.

Significant credit to Benedikt Diemer for developing Colossus:

https://bitbucket.org/bdiemer/colossus

Much of this is inspired by Colossus; any implementation faults are my own.
"""
module Halos

using Roots: fzero

using Slomo.Constants: default_cosmo
using Slomo.CosmologyTools: Ωm, ρm, ρcrit
import Slomo.Models: DensityModel, mass, density, NotImplemented

abstract type HaloModel <: DensityModel end

function scale_radius(halo::HaloModel)::Float64
    throw(NotImplemented("needs to be implemented for subtypes of HaloModel"))
end

"""
Parse the halo mass definition to get the spherical overdensity in Msun / kpc3.
If mdef is "vir", then use the spherical overdensity relation of Bryan & Norman 1998.

    mdef ∈ ("vir", "200c", "500c", "200m", "500m", ...)
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
Compute the virial radius from the virial mass, using the definition

    Mvir = 4π / 3 Rvir^3 * Δρ
"""
function Rvir_from_Mvir(Mvir; mdef = "200c", cosmo = default_cosmo, z = 0.0)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    return @. (Mvir * 3.0 / 4π / Δρ) ^ (1 / 3)
end

"""
Compute the virial mass from the virial radius, using the definition

    Mvir = 4π / 3 Rvir^3 * Δρ
"""
function Mvir_from_Rvir(Rvir; mdef = "200c", cosmo = default_cosmo, z = 0.0)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    return @. 4π / 3.0 * Rvir ^ 3 * Δρ
end

"""
Compute the virial radius.
"""
function virial_radius(halo::HaloModel;
                       mdef = "200c",
                       cosmo = default_cosmo,
                       z = 0.0,
                       xstart = 10.0,
                       rtol = 1e-3)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    # function to find roots
    f(r) = 3.0 / 4π * mass(halo, r) / r^3 - Δρ
    # derivative
    fp(r) = 3.0 * (r^-1 * density(halo, r) - r^-4 * mass(halo, r))
    rstart = xstart * scale_radius(halo)
    return fzero(f, fp, rstart; rtol = rtol)
end

"""
Compute the virial mass.
"""
function virial_mass(halo::HaloModel;
                     mdef = "200c",
                     cosmo = default_cosmo,
                     z = 0.0,
                     xstart = 10.0,
                     rtol = 1e-3)
    Rvir = virial_radius(halo;
                         mdef = mdef, cosmo = cosmo, z = z,
                         xstart = xstart, rtol = rtol)
    return Mvir_from_Rvir(Rvir, mdef = mdef, cosmo = cosmo, z = z)
end

"""
Compute the halo concentration from the z = 0 relation of Dutton & Maccio 2014.
"""
function concentration(Mvir; mdef = "200c", cosmo = default_cosmo, z = 0.0)
    mdef == lowercase(strip(mdef))
    @assert(mdef == "200c", "Currently only implemented for 200c definition")
    # convert Mvir to h-scaled units
    h = cosmo.h
    a = 0.520 + (0.905 - 0.520) * exp(-0.617 * z ^ 1.21)
    b = -0.101 + 0.026 * z
    return @. exp10(a + b * (log10(Mvir * h) - 12))
end

include("nfw.jl")

end
