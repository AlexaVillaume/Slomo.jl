using Slomo: CosmologyTools
using Slomo: Constants

using Test

cosmo = Constants.default_cosmo
rho_crit_h2 = Constants.rho_crit_h2

z = [0.0, 0.1, 0.5]
@test size(E(cosmo, z)) == size(z)
@test size(Ωm(z; cosmo = cosmo)) == size(z)
@test size(ρcrit(z; cosmo = cosmo)) == size(z)
@test size(ρm(z; cosmo = cosmo)) == size(z)

z = 0
@test Ωm(z; cosmo = cosmo) == cosmo.Ω_m
@test ρcrit(z; cosmo = cosmo) == rho_crit_h2 * cosmo.h ^ 2
@test ρm(z; cosmo = cosmo) == cosmo.Ω_m * rho_crit_h2 * cosmo.h ^ 2

z = 0.0
@test Ωm(z; cosmo = cosmo) == cosmo.Ω_m
@test ρcrit(z; cosmo = cosmo) == rho_crit_h2 * cosmo.h ^ 2
@test ρm(z; cosmo = cosmo) == cosmo.Ω_m * rho_crit_h2 * cosmo.h ^ 2
