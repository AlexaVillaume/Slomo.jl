"""
Copy of Cosmology.jl to remove Unitful dependence.

Original license:

The MIT License (MIT)

Copyright (c) 2013 Mike Nolta <mike@nolta.net>

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
"""
module CosmologyTools

using Slomo.Constants: rho_crit_h2
export Ωm, ρcrit, ρm, default_cosmo

abstract type AbstractCosmology end
abstract type AbstractClosedCosmology <: AbstractCosmology end
abstract type AbstractFlatCosmology <: AbstractCosmology end
abstract type AbstractOpenCosmology <: AbstractCosmology end

struct FlatLCDM{T <: Real} <: AbstractFlatCosmology
    h::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    FlatLCDM(promote(float(h), float(Ω_Λ), float(Ω_m), float(Ω_r))...)


a2E(c::FlatLCDM, a) = sqrt(c.Ω_r + c.Ω_m * a + c.Ω_Λ * a^4)

struct ClosedLCDM{T <: Real} <: AbstractClosedCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    ClosedLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                       float(Ω_r))...)


struct OpenLCDM{T <: Real} <: AbstractOpenCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) =
    OpenLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                     float(Ω_r))...)


function a2E(c::Union{ClosedLCDM,OpenLCDM}, a)
    a2 = a * a
    sqrt(c.Ω_r + c.Ω_m * a + (c.Ω_k + c.Ω_Λ * a2) * a2)
end

for c in ("Flat", "Open", "Closed")
    name = Symbol("$(c)WCDM")
    @eval begin
        struct $(name){T <: Real} <: $(Symbol("Abstract$(c)Cosmology"))
            h::T
            Ω_k::T
            Ω_Λ::T
            Ω_m::T
            Ω_r::T
            w0::T
            wa::T
        end
        function $(name)(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real,
                         w0::Real, wa::Real)
            $(name)(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                            float(Ω_r), float(w0), float(wa))...)
        end
    end
end

function WCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real, w0::Real, wa::Real)
    if Ω_k < 0
        ClosedWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    elseif Ω_k > 0
        OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    else
        FlatWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    end
end

function a2E(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    ade = exp((1 - 3 * (c.w0 + c.wa)) * log(a) + 3 * c.wa * (a - 1))
    sqrt(c.Ω_r + (c.Ω_m + c.Ω_k * a) * a + c.Ω_Λ * ade)
end

"""
    cosmology(;h = 0.69,
               Neff = 3.04,
               OmegaK = 0,
               OmegaM = 0.29,
               OmegaR = nothing,
               Tcmb = 2.7255,
               w0 = -1,
               wa = 0)
# Parameters
* `h` - Dimensionless Hubble constant
* `OmegaK` - Curvature density (Ω_k)
* `OmegaM` - Matter density (Ω_m)
* `OmegaR` - Radiation density (Ω_r)
* `Tcmb` - CMB temperature in Kelvin; used to compute Ω_γ
* `Neff` - Effective number of massless neutrino species; used to compute Ω_ν
* `w0` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
* `wa` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
# Examples
```jldoctest
julia> c = cosmology()
Cosmology.FlatLCDM{Float64}(0.69, 0.7099122024007928, 0.29, 8.77975992071536e-5)
julia> c = cosmology(OmegaK=0.1)
Cosmology.OpenLCDM{Float64}(0.69, 0.1, 0.6099122024007929, 0.29, 8.77975992071536e-5)
julia> c = cosmology(w0=-0.9, OmegaK=-0.1)
Cosmology.ClosedWCDM{Float64}(0.69, -0.1, 0.8099122024007929, 0.29, 8.77975992071536e-5, -0.9, 0.0)
```
"""
function cosmology(;h = 0.69,
                   Neff = 3.04,
                   OmegaK = 0,
                   OmegaM = 0.29,
                   OmegaR = nothing,
                   Tcmb = 2.7255,
                   w0 = -1,
                   wa = 0)

    if OmegaR === nothing
        OmegaG = 4.48131e-7 * Tcmb^4 / h^2
        OmegaN = Neff * OmegaG * (7 / 8) * (4 / 11)^(4 / 3)
        OmegaR = OmegaG + OmegaN
    end

    OmegaL = 1 - OmegaK - OmegaM - OmegaR

    if !(w0 == -1 && wa == 0)
        return WCDM(h, OmegaK, OmegaL, OmegaM, OmegaR, w0, wa)
    end

    if OmegaK < 0
        return ClosedLCDM(h, OmegaK, OmegaL, OmegaM, OmegaR)
        elseif OmegaK > 0
        return OpenLCDM(h, OmegaK, OmegaL, OmegaM, OmegaR)
        else
        return FlatLCDM(h, OmegaL, OmegaM, OmegaR)
    end
end

# hubble rate
scale_factor(z) = 1 / (1 + z)
E(c::AbstractCosmology, z) = (a = scale_factor(z); a2E(c, a) / a^2)

const planck18 = cosmology(h = 0.6766,
                           Neff = 3.046,
                           OmegaM = 0.3111)                           

const planck15 = cosmology(h = 0.6774,
                           Neff = 3.046,
                           OmegaM = 0.3089)

const planck13 = cosmology(h = 0.6777,
                           Neff = 3.046,
                           OmegaM = 0.3071)

const wmap9 = cosmology(h = 0.6932,
                        Neff = 3.046,
                        OmegaM = 0.2865)

default_cosmo = planck18

"""
Compute the matter density of the universe at redshift z, in units of the critical density.
"""
function Ωm(z; cosmo = default_cosmo)
    return cosmo.Ω_m * (1.0 + z) ^ 3 / E(cosmo, z) ^ 2
end

"""
Compute the critical density of the universe at redshift z, in Msun / kpc3
"""
function ρcrit(z; cosmo = default_cosmo)
    return rho_crit_h2 * cosmo.h^2 * E(cosmo, z) ^ 2
end

"""
Compute the mass density of the universe at redshift z, in Msun / kpc3
ρm(z) = Ωm(z) * ρcrit(z)
"""
function ρm(z; cosmo = default_cosmo)
    return cosmo.Ω_m * (1 + z)^3 * rho_crit_h2 * cosmo.h^2
end

end
