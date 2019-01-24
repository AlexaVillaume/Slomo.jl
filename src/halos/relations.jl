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
function shmr(Mvir; z = 0.0, cosmo = default_cosmo, mdef = default_mdef)

    if mdef != "vir"
        @warn("converting $mdef to vir halo mass definition", maxlog=1)
        @warn("assuming an NFW profile with standard concentration", maxlog=1)
        cvir = hmcr(Mvir; z = z, cosmo = cosmo, mdef = mdef)
        halo = NFW_from_virial(Mvir, cvir; z = z, cosmo = cosmo, mdef = mdef)
        Mvir = virial_mass(halo; z = z, cosmo = cosmo, mdef = "vir")
    end
    
    P(x, y, z) = y * z - x * z / (1.0 + z)
    Q(z) = exp(-4.0 / (1.0 + z)^2)
    
    ϵ0  =  -1.758
    ϵ1  =  +0.110
    ϵ2  =  -0.061
    ϵ3  =  -0.023
    M00 = +11.548
    M01 =  -1.297
    M02 =  -0.026
    α0  =  +1.975
    α1  =  +0.714
    α2  =  +0.042
    δ0  =  +3.390
    δ1  =  -0.472
    δ2  =  -0.931
    γ0  =  +0.498
    γ1  =  -0.157

    logϵ = ϵ0 + P(ϵ1, ϵ2, z) * Q(z) + P(ϵ3, 0.0, z)
    logM0 = M00 + P(M01, M02, z) * Q(z)
    α = α0 + P(α1, α2, z) * Q(z)
    δ = δ0 + P(δ1, δ2, z) * Q(z)
    γ = γ0 + P(γ1, 0.0, z) * Q(z)

    x = log10(Mvir) - logM0
    g(x) = δ * log10(1.0 + exp(x))^γ / (1.0 + exp(10.0^-x)) - log10(10^(-α * x) + 1.0)
    logMstar = logϵ + logM0 + g(x) - g(0)
    
    return exp10(logMstar)
end

function shmr(halo::HaloModel; mdef = default_mdef, z = 0.0, cosmo = default_cosmo)
    Mvir = virial_mass(halo; mdef = "vir", z = z, cosmo = cosmo)
    return shmr(Mvir; z = z, cosmo = cosmo, mdef = "vir")
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

"""
Computes the dimensionless peak height, ν = δ_crit(z) / σ(M, z), using the
approximation of Dutton & Maccio 2014 for a Planck 2013 cosmology.
"""
function peak_height(Mvir, z; cosmo = default_cosmo, mdef="200c")
    if mdef != "200c"
        throw("only implemented for 200c definition")
    end
    m = @. log10(M / 1e12 * cosmo.h)
    ν0 = @. exp10(-0.11 + 0.146 * m + 0.0138 * m^2 + 0.00123 * m^3)
    return @. ν0 * (0.033 + 0.79 * (1.0 + z) + 0.176 * exp(-1.356 * z))
end

"""
α-ν relation of Gao+2008
"""
function alpha_peak(Mvir; z = 0.0, cosmo = default_cosmo, mdef = "200c")
    if mdef != "200c"
        throw("only implemented for 200c definition")
    end
    ν = peak_height(Mvir, z; cosmo = cosmo, mdef = mdef)
    return @. 0.0095 * ν^2 + 0.155
end

function alpha_peak_prior(Mvir, alpha;
                          mdef = default_mdef, cosmo = default_cosmo, z = 0.0,
                          sigma0 = 0.16, sigmaz = 0.03)
    sigma_logalpha = sigma0 + sigmaz * z
    logalpha_expected = log10(alpha_peak(Mvir; z = z, cosmo = cosmo, mdef = mdef))
    logalpha = log10(alpha)
    return log_gauss(logalpha, logalpha_expected, sigma_logalpha)
end

    
