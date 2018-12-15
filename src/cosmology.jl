"""
Extensions to Cosmology.jl; might try to merge there eventually...

Significant credit to Benedikt Diemer for developing Colossus:

https://bitbucket.org/bdiemer/colossus

Much of this is inspired by Colossus; any implementation faults are my own.
"""
module CosmologyTools

import Cosmology
using Slomo.Constants: default_cosmo, rho_crit_h2

Cosmology.E(c::Cosmology.AbstractCosmology, z::Array{T} where T <: Real) = begin
    [E(c, i) for i in z]
end

"""
Compute the matter density of the universe at redshift z, in units of the critical density.
"""
function Ωm(z; cosmo = default_cosmo)
    return cosmo.Ω_m .* (1.0 .+ z) .^ 3 ./ Cosmology.E(cosmo, z) .^ 2
end

"""
Compute the critical density of the universe at redshift z, in Msun / kpc3
"""
function ρcrit(z; cosmo = default_cosmo)
    return rho_crit_h2 * cosmo.h^2 .* Cosmology.E(cosmo, z) .^ 2
end

"""
Compute the mass density of the universe at redshift z, in Msun / kpc3
"""
function ρm(z; cosmo = default_cosmo)
    return Ωm(z; cosmo = cosmo) .* ρcrit(z; cosmo = cosmo)
end

end
