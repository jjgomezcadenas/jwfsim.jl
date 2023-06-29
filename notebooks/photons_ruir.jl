### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c72e53d2-a21b-11ed-29d3-c72f8dda6658
using Pkg; Pkg.activate(ENV["JWfSim"])

# ╔═╡ 7926a440-16a6-4e7a-b9a4-ebc1102699e1
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Images
	using Colors
	using Clustering
	using Plots
	using Printf
	using InteractiveUtils
	using Statistics
	using StatsBase
	using StatsPlots
	using Distributions
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
	using PoissonRandom
	using ImageFiltering
	using ImageSegmentation
	using Interpolations
end

# ╔═╡ f9be2504-3c14-4ab2-bd67-90f3859cd96c
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    fs, ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M


# ╔═╡ 842e709a-8c93-4764-844b-e9193d6a2ed5
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 13513708-7ed8-4a6b-80c8-dd0456aea681
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
                :(eval(x) = $(Expr(:core, :eval))($name, x)),
                :(include(x) = $(Expr(:top, :include))($name, x)),
                :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
                :(include($path))))
    m
end

# ╔═╡ 6f0c7d83-04e6-456c-8a65-b6ecccf4ae0b
jwf = ingredients("../src/jwfsim.jl")

# ╔═╡ b0d80e2d-caaa-449d-a509-3487536dc9ae
PlutoUI.TableOfContents(title="Number of photons Ru & Ir", indent=true)

# ╔═╡ 78ff48a5-a903-44a3-af8d-ffc4c832cea5
md"""
# Absorption cross sections for Ru & Ir
"""

# ╔═╡ 6d047015-90c5-4774-9fa6-3e5f49160c19
md"""
## Load dataframes
"""

# ╔═╡ abf3d6d5-ca26-434c-b248-4c2395c79451
begin
	pdata = joinpath(ENV["JWfSim"], "data")
	mecrusl ="VS014_RuSL_1E-5_molar_extinction_coefficient.csv"
	mecirsl ="VS020_IrSL_1E-5_molar_extinction_coefficient.csv"
	ruslex  ="RuSl_emission.csv"
	irslex  ="IrSl_emission.csv"
end

# ╔═╡ 8bc26294-daad-409e-8eca-b03f1c79288e
md"""
#### Absorption spectrum 
"""

# ╔═╡ babadc24-94dd-4c00-80ae-42a0bbc5fcc7
begin
	mecrusldf = sort(jwf.jwfsim.load_df_from_csv(pdata, mecrusl, jwf.jwfsim.spG))
	names(mecrusldf)
	mecrusldf[!, "ϵ"] = mecrusldf[!, "e_M_1_cm_1"] .* cm^-1 .* M^-1
	mecrusldf[!, "λ"] = mecrusldf[!, "λ_nm"] .* nm
	mecrusldf

	mecirsldf = sort(jwf.jwfsim.load_df_from_csv(pdata, mecirsl, jwf.jwfsim.spG))
	mecirsldf[!, "ϵ"] = mecirsldf[!, "e_M_1_cm_1"] .* cm^-1 .* M^-1
	mecirsldf[!, "λ"] = mecirsldf[!, "λ_nm"] .* nm
	mecirsldf
end

# ╔═╡ c5bacabf-70a3-4c29-9562-a1a030c769f2
md"""
#### Emission spectrum 
"""

# ╔═╡ 41577c16-0c4b-4849-a265-3acc017bbf4c
begin
	ruslexdf = sort(jwf.jwfsim.load_df_from_csv(pdata, ruslex, jwf.jwfsim.spG))
	ruslexdf[!, "λ"] = ruslexdf[!, "λ_nm"] .* nm
	ruslexdf
end

# ╔═╡ 72b5d01a-b6f1-49ed-9674-7a51ae9ba1df
ruc = names(ruslexdf)[2:end-1]

# ╔═╡ 1d0779eb-73da-4163-b8bf-b4895ef4484e
begin
	irslexdf = sort(jwf.jwfsim.load_df_from_csv(pdata, irslex, jwf.jwfsim.spG))
	irslexdf[!, "λ"] = irslexdf[!, "λ_nm"] .* nm
	irslexdf
end

# ╔═╡ 4ee2f115-6188-413d-a3fa-04a22773cd8e
irc = names(irslexdf)[2:end-1]

# ╔═╡ e02432fe-7b16-42ae-af67-9959ad68e1a9
md"""
## Define fluorophores and cross sections
"""

# ╔═╡ ec58050a-51bd-49fc-83f3-77980db52d7c
begin
	rusl = jwf.jwfsim.Fluorophore(mecrusldf.λ[1], mecrusldf.λ[end], mecrusldf.ϵ, 0.11)
	xsrusl = jwf.jwfsim.xsecfl(rusl)
	irsl = jwf.jwfsim.Fluorophore(mecirsldf.λ[1], mecirsldf.λ[end], mecirsldf.ϵ, 0.47)
	xsirsl = jwf.jwfsim.xsecfl(irsl)
	md"""
	#### Plot absorption cross section 
	"""
end


# ╔═╡ c5d9874d-fcc4-4e01-95e0-83ba95089884
begin
	pirxsec = plot(mecirsldf.λ, xsirsl.(mecirsldf.λ)*1e+16, label= "Ir abs xsec * 1e+16", lw=2)
	xlims!(200.,550.)
	pruxsec = plot(mecrusldf.λ, xsrusl.(mecrusldf.λ)*1e+16, label= "Ru abs xsec * 1e+16", lw=2)
	xlims!(200.,550.)
	plot(pirxsec, pruxsec)
	
end

# ╔═╡ f95424f4-738f-48c4-b295-6d342732118b
md"""
#### Plot excitation spectrum 
"""

# ╔═╡ 42e03b5f-2e2c-4773-8653-28a87039be44
md""" Select Ru concentration : $(@bind scru Select(ruc))"""

# ╔═╡ 2ef97c9f-1c42-4b83-80be-2b6afe89ee8b
md""" Select Ir concentration : $(@bind scir Select(irc))"""

# ╔═╡ fc6837b1-7abe-4200-a322-e0ce3f3abba1
begin
	pru = plot(ruslexdf.λ, ruslexdf[!, scru], label= "Ru exc spectrum", lw=2)
	pir = plot(irslexdf.λ, irslexdf[!, scir], label= "Ir exc spectrum", lw=2)
	plot(pru, pir)
	#xlims!(200.,550.)

end

# ╔═╡ 42a0fe7f-2945-4804-ae7e-bb5763497853
md"""
# Epifluorescent microscope (narrow field)
"""

# ╔═╡ f69269ad-592c-4ea8-b228-1b97177299b4
md"""
## Define laser
"""

# ╔═╡ feaddf95-827c-417d-b6c0-59e86f839ae9
md" set laser wavelength in nm $(@bind ll NumberField(10.0^2:10.0^3; default=325.0))"

# ╔═╡ 157ebdfe-082f-475f-b4ba-1c022817d259
md" set laser power in mW $(@bind lp NumberField(1.0:1.0:10.0^3; default=1.0))"

# ╔═╡ 22442f67-6769-4f83-88b6-060f5a204748
lsr = jwf.jwfsim.CLaser(ll*nm, lp*mW)

# ╔═╡ ff3eb715-77b3-4306-aa7e-015aeeb6568d
md" set laser diameter in nm $(@bind ld NumberField(0.1:0.1:5.0; default=0.9))"

# ╔═╡ 50abd398-0450-4b98-9973-fef8277f1f16
d0 = ld*mm/2.0


# ╔═╡ 15a1900f-e886-444d-a041-1f72ed9885b4
md"""
The laser is initially a gaussian beam of waist $(d0/2)
"""

# ╔═╡ 6bb4cd1f-7e2c-4b1d-acd8-d344b59534d9
glsr = jwf.jwfsim.GaussianLaser(lsr, d0)

# ╔═╡ 7cedb3f7-97d7-465a-b9f1-993945b482b7
md"""
## Define objective
"""

# ╔═╡ 3ff0051b-11f3-4d59-915e-d9b47e88de17
md"""
- Current objective is a LMM-40x-VUV-Vac
"""

# ╔═╡ 91ca0f59-2ac9-475b-a477-c2406cccc0c8
obj = jwf.jwfsim.Objective(0.5, 40.0, 5.0mm, 5.1mm, 7.8mm, 0.85)

# ╔═╡ 3bc0e387-d305-4f19-9f62-8948d7d88256
md"""
## Compute the waist of the focused beam and define a new gaussian beam
"""

# ╔═╡ d8462918-6440-43eb-b8d4-279b619816da
w0dl = jwf.jwfsim.w0f(glsr, obj) 

# ╔═╡ 3de9e9aa-cbf7-4d81-ac3a-38c549beb3ba
glsrdl = jwf.jwfsim.GaussianLaser(lsr, w0dl)

# ╔═╡ aa9e5f13-28e3-424b-a6b6-0cfb0c1c7625
jwf.jwfsim.spot_size(glsrdl)

# ╔═╡ e7e7ac60-bded-48a6-b446-5aec1833cf4e
gI = jwf.jwfsim.I(glsrdl)

# ╔═╡ 4afd906b-dc50-4c78-ac4c-5b03d9533d4e
begin
	zr= 0.1:0.1:5.0
	ZZ = [z*cm for z in zr]
	pI = plot(ZZ, gI.(0.0*cm, ZZ), label= "Power at R=0 as a function of Z", lw=2)
end

# ╔═╡ 31de8cca-d53b-4d23-892b-294e19e0daf3
md"""
## Number of photons per molecule per unit time
"""

# ╔═╡ 34552b8d-2e63-411a-b231-5af142b7215d
md"""
- Emitted for Ru and Ir 
"""

# ╔═╡ ea0ccd60-90e0-4d83-91d7-59baba4d0063
ngrusl = glsrdl.ρ0 * xsrusl(ll*nm)

# ╔═╡ 8529dca1-32e2-4969-91c9-c97a95facfd3
ngirsl = glsrdl.ρ0 * xsirsl(ll*nm)

# ╔═╡ 19c06e59-f289-4c5b-9f80-9da19b6cb5fd
ng = [ngrusl, ngirsl]

# ╔═╡ 42818859-e2ed-4fc0-8c3a-913ffdefc6d7
4.6/59.0 * sqrt(59.0/511.0)

# ╔═╡ d31ffe9e-e9af-46f3-a689-5f9811767646
md"""
## Setup narrow field
"""

# ╔═╡ fc7160cf-f81e-4778-85e3-a27383f7b5c2
load("nf_feb_2023.png")

# ╔═╡ 45c79d1c-650e-4db7-8e9c-4d7cf4851fa9
md"""
- Trasnmitted through the objective 
"""

# ╔═╡ 8499afe3-5313-44cb-920f-99e4c2b2f13c
tobj = jwf.jwfsim.transmission(obj)

# ╔═╡ b0e163eb-9ffc-4dbd-b886-1b1ca8971ded
ngtobj = ng *tobj

# ╔═╡ cb9b6dfa-d54c-4635-9bff-1a89b73ad93a
md"""
- CCD efficiency
"""

# ╔═╡ 7e9efea1-061b-4d57-a6fb-b66a955d2926
eccd = jwf.jwfsim.ϵccd()

# ╔═╡ 8cd2b3cb-660b-48a5-a0e7-402fb835b499
begin
	wl = 350.0:10.0:1000.0
	wλ=[w*nm for w in wl]
	plot(wl, eccd.(wλ), label= "CCD efficiency", lw=2)
	#xlims!(200.,550.)

end

# ╔═╡ 3efab997-d8c3-4d7c-a67b-92ddd37d0a19
md"""
# TPA
"""

# ╔═╡ 896068bb-19e3-45af-a873-08972586cd67
md"""
### Quartz TPA cross section

- The TPA coefficient of commercial fused silica at 264 nm is computed [here](https://pubs.aip.org/aip/apl/article/80/7/1114/510979/Two-photon-absorption-properties-of-commercial)

- Using ths value (~10$^{-11}$ cm/W), the wavelength of the measurement (270 nm) and the density and atomic number of quartz (fused silica) we find a cross section which is at least two orders of magnitud smaller that what we expect for our molecules

"""

# ╔═╡ c2587cb0-f5c9-44fe-9473-5a36f8ed2bfb
ru2p = jwf.jwfsim.Fluorophore(mecrusldf.λ[1], mecrusldf.λ[end], mecrusldf.ϵ, 0.11,
	1.0e-50*cm^4*s)

# ╔═╡ 6dc8f788-5bc8-4942-a02b-16bc1413c554
md"""
### TPA 

- The number na of photons absorbed per fluorophore per pulse is:

$n_a = \frac{P^2 \sigma}{\tau f^2} (\frac{A^2}{2\hbar c \lambda})^2$

Where:

- P is the beam power   (e.g, mW)
- σ the TPA absorption cross section (cm$^4$ s)
- τ is the pulse duration (~100-200 fs)
- A is the N.A. of the setup
- λ is the laser wavelength (e.g, ~800 nm)


"""

# ╔═╡ 32f40dce-1aa7-4e82-9601-f557cd79b6c2
tt2 = uconvert(NoUnits, 1.0*cm^4* s* mW^2* MHz^-2* J^-2* nm^-2* m^-2* fs^-1)

# ╔═╡ 0f271b9e-b504-46e6-98c3-291d232d134c
md"""
#### Generic case:

- P = 50 mW
- σ = 1 GM (1 GM = 10$^{-50}$ cm$^4$ s)
- τ = 100 fs
- f = 80 MHz
- A = 1.4 (oil immersion)
- λ = 630 nm)

One achieve na ~1 (2 photons must be absorbed) per pulse. Assuming Q = 1 the rate of emitted photons is $\sim 3.7 \cdot 10^7$
"""

# ╔═╡ 414e6a62-b9f8-4ed7-b90b-9196edb62bfa
ltpa = jwf.jwfsim.PLaser(630.0*nm, 50.0*mW, 80.0*MHz, 100.0*fs)

# ╔═╡ b94300b7-ceb0-4405-b9b7-98af9527c77c
uconvert(NoUnits,1.0e-50*cm^4*s*((50.0*mW)^2/(100*fs*(80*MHz)^2))* 1.4^4/(2* Unitful.c0 * Unitful.ħ * 630*nm)^2)

# ╔═╡ c7ccec76-43ad-4332-8e64-4e9182eff08c
md"""
#### Origami case:

- P = 20 mW
- σ = 1 GM (1 GM = 10$^{-50}$ cm$^4$ s)
- τ = 160 fs
- f = 50 MHz
- A = 0.9 (air)
- λ = 780 nm)

One achieves na ~026 (2 photons must be absorbed) per pulse. Assuming Q = 1 the rate of emitted photons is $\sim 6.7 \cdot 10^5$
"""

# ╔═╡ b7a89e2a-b2a6-431a-b97e-90a7d0ac85a5
orgm = jwf.jwfsim.PLaser(780.0*nm, 20.0*mW, 50.0*MHz, 160.0*fs)

# ╔═╡ 72be96df-ec55-45ae-b728-941a36fdced0
uconvert(NoUnits,1.0e-50*cm^4*s*((20.0*mW)^2/(160*fs*(50*MHz)^2))* 0.9^4/(2* Unitful.c0 * Unitful.ħ * 780*nm)^2)

# ╔═╡ b5df9cf1-0226-442b-88ae-bf52395f8810
md"""
#### Origami XP:

- P = 4000 mW
- σ = 1 GM (1 GM = 10$^{-50}$ cm$^4$ s)
- τ = 370 fs
- f = 80 MHz
- A = 0.9 (air)
- λ = 1030 nm)

One achieves na ~026 (2 photons must be absorbed) per pulse. Assuming Q = 1 the rate of emitted photons is $\sim 6.7 \cdot 10^5$
"""

# ╔═╡ fa6e8970-35c3-4ecd-a2bb-10e0ef2bb3a2
orxp = jwf.jwfsim.PLaser(1030.0*nm, 500.0*mW, 80.0*MHz, 370.0*fs)

# ╔═╡ b3b8c4c6-c05b-47ea-b2ef-36f0e820a429
md"""
#### Carmel:

- P = 1000 mW
- σ = 1 GM (1 GM = 10$^{-50}$ cm$^4$ s)
- τ = 90 fs
- f = 80 MHz
- A = 0.9 (air)
- λ = 780 nm

One achieves na ~026 (2 photons must be absorbed) per pulse. Assuming Q = 1 the rate of emitted photons is $\sim 6.7 \cdot 10^5$
"""

# ╔═╡ c63a0509-501a-4363-a0e0-c3ccbfcbaeca
carml = jwf.jwfsim.PLaser(780.0*nm, 700.0*mW, 80.0*MHz, 90.0*fs)

# ╔═╡ bd400f0e-427e-4de8-971f-0439f2d8703a
jwf.jwfsim.n_photons(carml)^2

# ╔═╡ bb25bf7f-28e7-419a-8249-b6484aae25e9
N_A

# ╔═╡ 1279174a-3e68-4dad-9b41-7c798807bd44
md"""
## Functions
``\sigma``
"""

# ╔═╡ 145caf25-418e-4198-9868-213ca0ee405a
md"""
Computes the two photon absorption (TPA cross section(σ) from the TPA coefficient (β)

$\sigma = \beta \times \frac{E}{N}$

Where $\beta$ is the TPA coefficient (in cm/W) N is the molecular density of the material (molecules/cm3) and E is the energy of the laser (in J). Then \sigma is the TPA cross section in $cm^4 \cdot s \cdot photons^{-1} \cdot molecules^{-1}$

Units:

$[cm \cdot J^{-1} \cdot s] \cdot[J \cdot photons^{-1}] \cdot[molecules^{-1} cm^3]$

"""

# ╔═╡ 756445a8-8c57-475f-9789-d911866368d4
"""
Computes the two photon absorption (TPA cross section(σ) from the TPA coefficient (β)

- λ is the wavelength of the laser light
- ρ is the density of the material
- A is the atomic weight 

"""
function tpa_σ_from_β(β::typeof(1.0cm/W), λ::typeof(1.0nm), ρ::typeof(1.0g/cm^3), A::typeof(1.0g/mol))
	
	N = N_A*ρ/A
	E = uconvert(J, λ, Spectral())
	β * E/N
	uconvert(cm^4*s, β * E/N)
end

# ╔═╡ c903957b-8a9d-4a8b-b2f3-56d4e1a21718
tpa_σ_from_β(1e-11*cm/W, 270.0*nm, 2.2*(g/cm^3),60.0*g/mol)


# ╔═╡ 41811b42-21e9-47ff-9b94-da4fbdb56f92
function na_tpa(f::jwf.jwfsim.Fluorophore, lb::jwf.jwfsim.Laser, NA::Float64)
	hbc = Unitful.c0 * Unitful.ħ
	t1 = lb.P^2 * f.σ /(lb.τ * lb.f^2)
	t2 = NA^2/(2.0 * hbc * lb.λ)
	uconvert(NoUnits, t1*t2^2)
end

# ╔═╡ 60bc5945-3608-4fb3-9252-d55553185154
na_tpa(ru2p,ltpa, 1.4)

# ╔═╡ 13963ccb-2cda-417e-b431-686be59e5e05
ltpa.f * na_tpa(ru2p,ltpa, 1.4)/2.0

# ╔═╡ 68708dc7-f0b1-44b3-bec0-6a5e860e0633
na_tpa(ru2p,orgm, 0.9)

# ╔═╡ 71422aee-46c9-4373-880f-e0747a477c58
orgm.f * na_tpa(ru2p,orgm, 0.9) /2.0

# ╔═╡ 09b005f8-a256-47ce-8785-92d24a45ed8c
na_tpa(ru2p,orxp, 0.9)

# ╔═╡ 54f0e83f-3855-4ff3-aa47-c994387ba78c
[na_tpa(ru2p,carml, NA) for NA in 0.5:0.05:1.00]

# ╔═╡ ee7b08bf-d9fb-4d0c-a7e0-a72de3508620
function fluorescence_2p_dl(C::typeof(1.0*mol/μm^3),
	                        m::jwf.jwfsim.Fluorophore,
                       	    lb::jwf.jwfsim.Laser)
	gp = 0.664
	n=1
	uconvert(Hz, 0.5 * m.Q * m.σ * C * N_A * (gp/(lb.f * lb.τ)) * 8 * n * jwf.jwfsim.n_photons(lb)^2/(π*lb.λ))
end

# ╔═╡ 96369fde-708d-4478-b6d5-82d4b5bf1552
fluorescence_2p_dl(N_A^-1/μm^3, ru2p, carml)

# ╔═╡ 6978a5db-fbcd-4862-b59f-e94830e0e28c
md"""
# Hector
"""

# ╔═╡ bcf23d84-d93d-43f5-a0cb-34c710e9160c
dz(R, G, M, m, γ) = - 2.0 * γ * sqrt(G *M) /(R * m)

# ╔═╡ 666e01f6-f1d5-4783-b2a0-6d24b666e023
dz(6.4e+6*m, 6.67e-11*N*m^2/kg^2, 6.0e+24*kg, 140.0e+3kg, 1.0e-8*kg/m)

# ╔═╡ 1841f2b4-5e1c-454c-be2a-70add51652d8
uconvert(s^-1, 1.0*N^(1.0/2.0)*kg^(-1.0/2.0)*m^-1)

# ╔═╡ 277be488-9840-4f62-874f-d4f9fbb893d8


# ╔═╡ Cell order:
# ╠═c72e53d2-a21b-11ed-29d3-c72f8dda6658
# ╠═7926a440-16a6-4e7a-b9a4-ebc1102699e1
# ╠═f9be2504-3c14-4ab2-bd67-90f3859cd96c
# ╠═842e709a-8c93-4764-844b-e9193d6a2ed5
# ╠═13513708-7ed8-4a6b-80c8-dd0456aea681
# ╠═6f0c7d83-04e6-456c-8a65-b6ecccf4ae0b
# ╠═b0d80e2d-caaa-449d-a509-3487536dc9ae
# ╠═78ff48a5-a903-44a3-af8d-ffc4c832cea5
# ╠═6d047015-90c5-4774-9fa6-3e5f49160c19
# ╠═abf3d6d5-ca26-434c-b248-4c2395c79451
# ╠═8bc26294-daad-409e-8eca-b03f1c79288e
# ╠═babadc24-94dd-4c00-80ae-42a0bbc5fcc7
# ╠═c5bacabf-70a3-4c29-9562-a1a030c769f2
# ╠═41577c16-0c4b-4849-a265-3acc017bbf4c
# ╠═72b5d01a-b6f1-49ed-9674-7a51ae9ba1df
# ╠═1d0779eb-73da-4163-b8bf-b4895ef4484e
# ╠═4ee2f115-6188-413d-a3fa-04a22773cd8e
# ╠═e02432fe-7b16-42ae-af67-9959ad68e1a9
# ╠═ec58050a-51bd-49fc-83f3-77980db52d7c
# ╠═c5d9874d-fcc4-4e01-95e0-83ba95089884
# ╠═f95424f4-738f-48c4-b295-6d342732118b
# ╠═42e03b5f-2e2c-4773-8653-28a87039be44
# ╠═2ef97c9f-1c42-4b83-80be-2b6afe89ee8b
# ╠═fc6837b1-7abe-4200-a322-e0ce3f3abba1
# ╠═42a0fe7f-2945-4804-ae7e-bb5763497853
# ╠═f69269ad-592c-4ea8-b228-1b97177299b4
# ╠═feaddf95-827c-417d-b6c0-59e86f839ae9
# ╠═157ebdfe-082f-475f-b4ba-1c022817d259
# ╠═22442f67-6769-4f83-88b6-060f5a204748
# ╠═ff3eb715-77b3-4306-aa7e-015aeeb6568d
# ╠═15a1900f-e886-444d-a041-1f72ed9885b4
# ╠═50abd398-0450-4b98-9973-fef8277f1f16
# ╠═6bb4cd1f-7e2c-4b1d-acd8-d344b59534d9
# ╠═7cedb3f7-97d7-465a-b9f1-993945b482b7
# ╠═3ff0051b-11f3-4d59-915e-d9b47e88de17
# ╠═91ca0f59-2ac9-475b-a477-c2406cccc0c8
# ╠═3bc0e387-d305-4f19-9f62-8948d7d88256
# ╠═d8462918-6440-43eb-b8d4-279b619816da
# ╠═3de9e9aa-cbf7-4d81-ac3a-38c549beb3ba
# ╠═aa9e5f13-28e3-424b-a6b6-0cfb0c1c7625
# ╠═e7e7ac60-bded-48a6-b446-5aec1833cf4e
# ╠═4afd906b-dc50-4c78-ac4c-5b03d9533d4e
# ╠═31de8cca-d53b-4d23-892b-294e19e0daf3
# ╠═34552b8d-2e63-411a-b231-5af142b7215d
# ╠═ea0ccd60-90e0-4d83-91d7-59baba4d0063
# ╠═8529dca1-32e2-4969-91c9-c97a95facfd3
# ╠═19c06e59-f289-4c5b-9f80-9da19b6cb5fd
# ╠═42818859-e2ed-4fc0-8c3a-913ffdefc6d7
# ╠═d31ffe9e-e9af-46f3-a689-5f9811767646
# ╠═fc7160cf-f81e-4778-85e3-a27383f7b5c2
# ╠═45c79d1c-650e-4db7-8e9c-4d7cf4851fa9
# ╠═8499afe3-5313-44cb-920f-99e4c2b2f13c
# ╠═b0e163eb-9ffc-4dbd-b886-1b1ca8971ded
# ╠═cb9b6dfa-d54c-4635-9bff-1a89b73ad93a
# ╠═7e9efea1-061b-4d57-a6fb-b66a955d2926
# ╠═8cd2b3cb-660b-48a5-a0e7-402fb835b499
# ╠═3efab997-d8c3-4d7c-a67b-92ddd37d0a19
# ╠═896068bb-19e3-45af-a873-08972586cd67
# ╠═c903957b-8a9d-4a8b-b2f3-56d4e1a21718
# ╠═c2587cb0-f5c9-44fe-9473-5a36f8ed2bfb
# ╠═6dc8f788-5bc8-4942-a02b-16bc1413c554
# ╠═32f40dce-1aa7-4e82-9601-f557cd79b6c2
# ╠═0f271b9e-b504-46e6-98c3-291d232d134c
# ╠═414e6a62-b9f8-4ed7-b90b-9196edb62bfa
# ╠═60bc5945-3608-4fb3-9252-d55553185154
# ╠═b94300b7-ceb0-4405-b9b7-98af9527c77c
# ╠═13963ccb-2cda-417e-b431-686be59e5e05
# ╠═c7ccec76-43ad-4332-8e64-4e9182eff08c
# ╠═b7a89e2a-b2a6-431a-b97e-90a7d0ac85a5
# ╠═68708dc7-f0b1-44b3-bec0-6a5e860e0633
# ╠═72be96df-ec55-45ae-b728-941a36fdced0
# ╠═71422aee-46c9-4373-880f-e0747a477c58
# ╠═b5df9cf1-0226-442b-88ae-bf52395f8810
# ╠═fa6e8970-35c3-4ecd-a2bb-10e0ef2bb3a2
# ╠═09b005f8-a256-47ce-8785-92d24a45ed8c
# ╠═b3b8c4c6-c05b-47ea-b2ef-36f0e820a429
# ╠═c63a0509-501a-4363-a0e0-c3ccbfcbaeca
# ╠═54f0e83f-3855-4ff3-aa47-c994387ba78c
# ╠═bd400f0e-427e-4de8-971f-0439f2d8703a
# ╠═96369fde-708d-4478-b6d5-82d4b5bf1552
# ╠═bb25bf7f-28e7-419a-8249-b6484aae25e9
# ╠═1279174a-3e68-4dad-9b41-7c798807bd44
# ╠═145caf25-418e-4198-9868-213ca0ee405a
# ╠═756445a8-8c57-475f-9789-d911866368d4
# ╠═41811b42-21e9-47ff-9b94-da4fbdb56f92
# ╠═ee7b08bf-d9fb-4d0c-a7e0-a72de3508620
# ╠═6978a5db-fbcd-4862-b59f-e94830e0e28c
# ╠═bcf23d84-d93d-43f5-a0cb-34c710e9160c
# ╠═666e01f6-f1d5-4783-b2a0-6d24b666e023
# ╠═1841f2b4-5e1c-454c-be2a-70add51652d8
# ╠═277be488-9840-4f62-874f-d4f9fbb893d8
