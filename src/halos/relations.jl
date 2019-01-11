"""
Compute the expected halo concentration from halo mass using the relation of 
Dutton & Maccio 2014 (equations 7, 10-13)
"""
function hmcr(Mvir; mdef = default_mdef, cosmo = default_cosmo, z = 0.0)
    mdef == lowercase(strip(mdef))
    if mdef == "200c"
        a = 0.520 + (0.905 - 0.520) * exp(-0.617 * z ^ 1.21)
        b = -0.101 + 0.026 * z
    elseif mdef == "vir"
        a = 0.537 + (1.025 - 0.537) * exp(-0.718 * z ^ 1.08)
        b = -0.097 + 0.024 * z
    else
        throw("only implemented for 200c and vir halo mass definitions")
    end
    # convert Mvir to h-scaled units
    return @. exp10(a + b * (log10(Mvir * cosmo.h) - 12.0))
end

"""
Compute the probability of drawing a halo with virial parameters (Mvir, cvir),
using the halo mass-concentration relation of Dutton & Maccio 2014.

     sigma_logc : scatter in logc at fixed Mvir (in dex)
"""
function hmcr_prior(Mvir, cvir;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logc = 0.16)
    logc_expected = log10(hmcr(Mvir; mdef = mdef, cosmo = cosmo, z =z))
    logc = log10(cvir)
    return log_gauss(logc, logc_expected, sigma_logc)
end

function hmcr_prior(halo::HaloModel;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logc = 0.16)
    Mvir = virial_mass(halo; mdef = mdef, z = z, cosmo = cosmo)
    logc_expected = log10(hmcr(Mvir; mdef = mdef, cosmo = cosmo, z =z))
    logc_expected = log10(cvir_expected)
    logc = log10(concentration(halo; mdef = mdef, z = z, cosmo = cosmo))
    return log_gauss(logc, logc_expected, sigma_logc)
end


"""
Compute the expected stellar mass from halo mass using the relation of 
Rodriguez-Puebla et al. 2017.

See equations 25-33 for the functional form and 49 - 55 for the parameters.
"""
function shmr(halo::HaloModel; z = 0.0, cosmo = default_cosmo)
    
    Mvir = virial_mass(halo; mdef = "vir", z = z, cosmo = cosmo)
    
    P(x, y, z) = y * z - x * z / (1.0 + z)
    Q(z) = exp(-4.0 / (1.0 + z)^2)
    
    ϵ0     = -1.758
    ϵ1     = +0.110
    ϵ2     = -0.061
    ϵ3     = -0.023
    logM00 = +11.548
    M01    = -1.297
    M02    = -0.026
    α0     = +1.975
    α1     = +0.714
    α2     = +0.042
    δ0     = +3.390
    δ1     = -0.472
    δ2     = -0.931
    γ0     = +0.498
    γ1     = -0.157

    logϵ = ϵ0 + P(ϵ1, ϵ2, z) * Q(z) + P(ϵ3, 0.0, z)
    logM0 = logM00 + P(M01, M02, z) * Q(z)
    α = α0 + P(α1, α2, z) * Q(z)
    δ = δ0 + P(δ1, δ2, z) * Q(z)
    γ = γ0 + P(γ1, 0.0, z) * Q(z)

    x = Mvir / exp10(logM0)
    g(x) = δ * log10(1.0 + exp(x))^γ / (1.0 + exp(10.0^-x)) - log10(10^(-α * x) + 1.0)
    logMstar = logϵ + logM0 + g(x) - g(0)
    
    return exp10(logMstar)
end

function shmr(Mvir; mdef = default_mdef, z = 0.0, cosmo = default_cosmo)
    mdef == lowercase(strip(mdef))
    if mdef != "vir"
        @warn("converting halo mass from $mdef to vir definition", maxlog = 1)
        @warn("assuming an NFW profile", maxlog = 1)
        cvir = hmcr(Mvir; mdef = mdef, z = z, cosmo = cosmo)
        halo = NFW_from_virial(Mvir, cvir; mdef = mdef, z = z, cosmo = cosmo)
        Mvir = virial_mass(halo; mdef = "vir", z = z, cosmo = cosmo)
    end
    return shmr(halo; z = z, cosmo = cosmo)
end

"""
Compute the probability of drawing a halo mass and stellar mass pair (Mvir, Mstar).

    sigma_logMstar : scatter in logMstar (dex)
"""
function shmr_prior(Mvir, Mstar;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logMstar = 0.15)
    logMstar_expected = log10(shmr(Mvir; mdef = mdef, z = z, cosmo = cosmo))
    logMstar = log10(Mstar)
    return log_gauss(logMstar, logMstar_expected, sigma_logMstar)
end

function shmr_prior(halo::HaloModel, Mstar;
                    mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                    sigma_logMstar = 0.15)
    logMstar_expected = log10(shmr(halo; z = z, cosmo = cosmo))
    logMstar = log10(Mstar)
    return log_gauss(logMstar, logMstar_expected, sigma_logMstar)
end
