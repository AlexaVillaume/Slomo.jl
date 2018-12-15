"""
Dark matter halo models and relations.

Significant credit to Benedikt Diemer for developing Colossus:

https://bitbucket.org/bdiemer/colossus

Much of this is inspired by Colossus; any implementation faults are my own.
"""
module Halos

using Slomo.Constants: default_cosmo
using Slomo.CosmologyTools: Ωm, ρm, ρcrit
import Slomo.Models: DensityModel, mass, density

include("halo_models/nfw.jl")

"""
Parse the halo mass definition to get the spherical overdensity in Msun / kpc3.
If mdef is "vir", then use the spherical overdensity relation of Bryan & Norman 1998.

    mdef ∈ ("vir", "200c", "500c", "200m", "500m", ...)
"""
function Δρ_from_mdef(mdef; cosmo = default_cosmo, z = 0.0)
    mdef = lowercase(strip(mdef))
    if mdef == "vir"
        x = Ωm(z; cosmo = cosmo) .- 1.0
        Δ = @. 18.0 * pi^2 + 82.0 * x - 39.0 * x^2
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

    Mvir = 4pi / 3 Rvir^3 * delta_rho
"""
function Rvir_from_Mvir(Mvir; mdef = "200c", cosmo = default_cosmo, z = 0.0)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    return @. (Mvir * 3.0 / 4pi / Δρ) ^ (1 / 3)
end

"""
Compute the virial mass from the virial radius, using the definition

    Mvir = 4pi / 3 Rvir^3 * delta_rho
"""
function Mvir_from_Rvir(Rvir; mdef = "200c", cosmo = default_cosmo, z = 0.0)
    Δρ = Δρ_from_mdef(mdef; cosmo = cosmo, z = z)
    return @. 4pi / 3.0 * Rvir ^ 3 * Δρ
end

end
