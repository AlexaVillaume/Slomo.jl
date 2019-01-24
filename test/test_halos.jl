using Slomo: Halos
using Slomo: Constants

using Test

rtol = 1e-5

for T in [NFWModel, GNFWModel, CoreNFWModel, EinastoModel, SolNFWModel]
    for Mvir in exp10.([8.0, 10.0, 12.0, 14.0])
        for mdef in ["vir", "200c", "200m", "500c"]
            for z in [0.0, 0.5, 1.0]
                if mdef == "vir" || mdef = "200c"
                    cvir = Halos.hmcr(Mvir; mdef = mdef, z = z)
                else
                    cvir = 10.0
                end
                if T == NFWModel
                    halo = Halos.NFW_from_virial(Mvir, cvir; mdef = mdef, z = z)
                elseif T == GNFWModel
                    halo = Halos.GNFW_from_virial(Mvir, cvir, 0.3; mdef = mdef, z = z)
                elseif T == CoreNFWModel
                    halo = Halos.CoreNFW_from_virial(Mvir, cvir, 10.0, 10.0; mdef = mdef, z = z)
                elseif T == EinastoModel
                    halo = Halos.Einasto_from_virial(Mvir, cvir, 0.16; mdef = mdef, z = z)
                elseif T == SolNFWModel
                    halo = Halos.SolNFW_from_virial(Mvir, cvir, 1.0; mdef = mdef, z = z)
                end
                rs = scale_radius(halo)
                Rvir = virial_radius(halo; mdef = mdef, z = z)            
                @test isapprox(Mvir, mass(halo, Rvir), rtol = rtol)
                @test isapprox(cvir, Rvir / rs, rtol = rtol)
                @test isapprox(Mvir, virial_mass(halo; mdef = mdef, z = z), rtol = rtol)
                @test isapprox(cvir, concentration(halo; mdef = mdef, z = z), rtol = rtol)
                @test isapprox(Mvir, Mvir_from_Rvir(Rvir; mdef = mdef, z = z), rtol = rtol)
                @test isapprox(Rvir, Rvir_from_Mvir(Mvir; mdef = mdef, z = z), rtol = rtol)
                for rs_fraction in [1, 1.0, [1], [1.0], [0.5, 1.0, 1.0]]
                    r = rs .* rs_fraction
                    @test size(mass(halo, r)) == size(r)
                    @test size(density(halo, r)) == size(r)
                end
            end
        end
    end
end
