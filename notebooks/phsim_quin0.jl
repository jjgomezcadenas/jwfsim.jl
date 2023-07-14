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
	μW, mW, W, kW,
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
PlutoUI.TableOfContents(title="Quin", indent=true)

# ╔═╡ 78ff48a5-a903-44a3-af8d-ffc4c832cea5
md"""
# Simulation steps
"""

# ╔═╡ 6d047015-90c5-4774-9fa6-3e5f49160c19
md"""
## Absorption and excitation spectra for QUINN molecules
"""

# ╔═╡ aadf254a-47ea-417f-9f5e-734486445e32
md"""
### Define fluorophores and cross sections
"""

# ╔═╡ abf3d6d5-ca26-434c-b248-4c2395c79451
begin
	pdata = joinpath(ENV["JWfSim"], "data")
	mequin ="QUIN_1E-5M_molar_extinction_coefficient.csv"
	
end

# ╔═╡ babadc24-94dd-4c00-80ae-42a0bbc5fcc7
begin
	mequindfr = sort(jwf.jwfsim.load_df_from_csv(pdata, mequin, jwf.jwfsim.spG))
	names(mequindfr)
	mequindf = select(mequindfr, [:l, :c, :e])
	mequindf[!, "ϵ"] = mequindf[!, "e"] .* 1.0 * cm^-1 .* M^-1
	mequindf[!, "λ"] = mequindf[!, "l"] .* 1.0 * nm
	mequindf[!, "C"] = mequindf[!, "c"] .* 1E-5* M
	mqdf = select(mequindf, Not([:l, :c, :e]))
	rusl = jwf.jwfsim.Fluorophore(mqdf.λ[1], mqdf.λ[end], mqdf.ϵ, 0.03)
	xquin = jwf.jwfsim.xsecfl(rusl)
	md"""
	- Cross section at 375 nm $(xquin(375.0*nm))
	"""
end

# ╔═╡ e02432fe-7b16-42ae-af67-9959ad68e1a9
md"""
#### Absorption spectrum 
"""

# ╔═╡ ec58050a-51bd-49fc-83f3-77980db52d7c
begin
	pxsec = plot(mqdf.λ, xquin.(mqdf.λ)*1e+16, label= "QUIN abs xsec * 1e+16", lw=2)
	xlims!(200.,750.)
	plot(pxsec)
end


# ╔═╡ 42a0fe7f-2945-4804-ae7e-bb5763497853
md"""
# Setup
"""

# ╔═╡ f69269ad-592c-4ea8-b228-1b97177299b4
md"""
## Define laser
"""

# ╔═╡ feaddf95-827c-417d-b6c0-59e86f839ae9
md" set laser wavelength in nm $(@bind ll NumberField(10.0^2:10.0^3; default=370.0))"

# ╔═╡ 157ebdfe-082f-475f-b4ba-1c022817d259
md" set laser power in muW $(@bind lp NumberField(1.0:1.0:10.0^2; default=7.0))"

# ╔═╡ ad2254eb-efe5-44b5-a1d1-c023a4456801
md" set laser repetition rate in MHz $(@bind ff NumberField(0.5:0.5:10.0^2; default=1.0))"

# ╔═╡ 22442f67-6769-4f83-88b6-060f5a204748
lsr = jwf.jwfsim.PLaser(1.0*ll*nm, 1.0*lp*μW, 1.0*ff*MHz, 1e+3*fs)

# ╔═╡ 31c8c858-2178-4fce-a772-d116dca3f496
md"""
## Define objective
- Current objective is a LMM-40x-VUV-Vac
"""

# ╔═╡ a3d8ba37-07dd-4424-9c6f-e4b453205b90
obj = jwf.jwfsim.Objective(0.5, 40.0, 5.0mm, 5.1mm, 7.8mm, 0.85)

# ╔═╡ 0a2b607f-4e30-41be-a73d-c0daf5c8b855
md"""
## Photon rate per molecule (assuming diffractive limit formulas)
"""

# ╔═╡ ebb200f1-f5ab-4b66-a7ce-c2173fd1fc8f
begin
	dl  = jwf.jwfsim.diffractive_limit(lsr.λ, obj.NA; ff=1.22)
	zr  = jwf.jwfsim.diffractive_limit(lsr.λ, obj.NA)
	fov = jwf.jwfsim.Fov(dl, 2*zr)
	I   = uconvert(Hz/cm^2, jwf.jwfsim.photon_density(lsr, fov))
	ngdl = I * xquin(ll*nm)
	md"""
	##### Diffractive limit assumed
	- dl = $dl
	- zr = $zr
	- I = $I
	- Photon rate = $ngdl
	"""
end

# ╔═╡ 5a69f27c-86b2-48bb-99b1-357fbee639bf
md"""
### Photon rate per molecule using a gaussian beam

- Define the waist of a gaussian beam which fils the objective
"""

# ╔═╡ e672f7bc-5020-40ae-9706-a12640607b6f
w0 = jwf.jwfsim.w0f(lsr.λ, obj.d, obj.f)

# ╔═╡ 6bb4cd1f-7e2c-4b1d-acd8-d344b59534d9
glsr = jwf.jwfsim.GaussianLaser(lsr, w0)

# ╔═╡ 8375a215-ce9c-4a38-96a5-861e0f452282
Igl = jwf.jwfsim.I(glsr)

# ╔═╡ 1d8692e3-7cf4-4739-a750-3f713512979b
uconvert(kW/cm^2, Igl(0.0nm,0.0nm))

# ╔═╡ aa9e5f13-28e3-424b-a6b6-0cfb0c1c7625
ng = glsr.ρ0 * xquin(ll*nm)

# ╔═╡ 05e4e447-4012-46e2-a6cf-e14c49e31d40
md"""
	##### Assming a gaussian beam that fills the objective
	- dl = $(jwf.jwfsim.spot_size(glsr))
	- Photon rate = $ng
	"""

# ╔═╡ 31de8cca-d53b-4d23-892b-294e19e0daf3
md"""
#### Number of photons per molecule per unit time
"""

# ╔═╡ ea0ccd60-90e0-4d83-91d7-59baba4d0063
ng/ ngdl

# ╔═╡ b6339888-e34c-47e4-8af1-2e9393461433
md"""
### Number of photons per molecule per pulse
"""

# ╔═╡ 663f5495-7e50-4a58-b57f-db915a17fd15
ngp = uconvert(NoUnits, ng/lsr.f)

# ╔═╡ fc7160cf-f81e-4778-85e3-a27383f7b5c2
#load("nf_feb_2023.png")

# ╔═╡ 45c79d1c-650e-4db7-8e9c-4d7cf4851fa9
md"""
- photons Transmitted through the objective 
"""

# ╔═╡ 8499afe3-5313-44cb-920f-99e4c2b2f13c
tobj = jwf.jwfsim.transmission(obj)

# ╔═╡ b0e163eb-9ffc-4dbd-b886-1b1ca8971ded
ngtobj = ng *tobj

# ╔═╡ cb9b6dfa-d54c-4635-9bff-1a89b73ad93a
md"""
- APD efficiency
"""

# ╔═╡ 7e9efea1-061b-4d57-a6fb-b66a955d2926
eapd = jwf.jwfsim.ϵapd()

# ╔═╡ 8cd2b3cb-660b-48a5-a0e7-402fb835b499
begin
	wl = 350.0:10.0:1000.0
	wλ=[w*nm for w in wl]
	plot(wl, eapd.(wλ), label= "APD efficiency", lw=2)
	#xlims!(200.,550.)

end

# ╔═╡ ac8d75b1-4804-423b-8268-124ead6366eb
md"""
- Approximate average efficiency for QUIN ~60 %
"""

# ╔═╡ 9fda79b5-f02d-479a-a7e8-9e2e6ca3e8fc
ngtoapd = ng *tobj * 0.6

# ╔═╡ 38947c5e-aa3f-427d-bee1-a4532693381e
ngpx = uconvert(NoUnits, ngtoapd/lsr.f)

# ╔═╡ 7aef4cca-a6ce-4f16-ba7d-e7914d054247
md"""
- Number of photons per molecule per unit time =$ng
- Transmitted through the objective =$ngtobj
- Transmitted through the apd =$ngtoapd
- Number of photons per molecule per pulse = $ngp
- Detected, per pulse =$ngpx
"""

# ╔═╡ 61ffc4ea-5107-453a-8048-9bdc11a94659
md" set DC in Hz $(@bind dchz NumberField(1:0.5:10.0^2; default=20))"

# ╔═╡ d3810185-890e-4c52-a947-e53d6d9e0ebc
ngdcx = dchz/1e+6

# ╔═╡ 2e403946-4a25-45ee-935b-2f6af7f87658
md"""
- DC (Hz) =$dchz
- DC (per pulse) =$ngdcx
"""

# ╔═╡ 1ae83c20-d8ec-4717-9275-8828bc1e193a
md"""
### Distribution per pulse of signal and DC
"""

# ╔═╡ 2d68c170-0d91-4d84-975a-0165b86a2883
dsgn = Poisson(ngpx)

# ╔═╡ 4f7a9fc1-9534-42fa-890d-f298fba635e7
ddc = Poisson(ngdcx)

# ╔═╡ eedef1fc-4b5a-49b8-8f61-ae99cc27215e
md"""
### Time distribution of signal
"""

# ╔═╡ 57da5efb-d21b-4819-a623-32ea162e6456
dexpmus = Exponential(1.0)

# ╔═╡ ea7e2262-ca35-4537-b914-5e1971cbfc20
md"""
### Generate the events and the times of signal and DC
"""

# ╔═╡ 3d3fa4a0-be5c-450f-9b0a-8500878f38b1
md" set acquistion time (in seconds) $(@bind act NumberField(1:0.5:10.0^2; default=10))"

# ╔═╡ 0588ae25-2f71-48d1-9fbf-e682f6161f24
ntx = uconvert(NoUnits, act*s*lsr.f)

# ╔═╡ e38ad10c-f0c2-4fd8-b98e-27ef6b7b591e
evtsgn = rand(dsgn, Int(ntx))

# ╔═╡ 16f91e51-2f6b-4431-9994-b471a1f3994e
evtdc = rand(ddc, Int(ntx))

# ╔═╡ e7c14906-1994-4bfc-b9cd-2f9ccc8cd2ad
tsgn = rand(dexpmus, Int(ntx))

# ╔═╡ 540cebd9-7630-4f35-9bf6-091516eda6bb
tdc = rand(Int(ntx))

# ╔═╡ c90c9b5a-ff51-46b0-9ea1-3ec1c745d68a
md"""
#### Find the indexes of events which are not zero (thus one)
"""

# ╔═╡ 5a70603e-ce7a-4518-ad58-86c720825267
ixevtsgn = findall(!iszero, evtsgn)


# ╔═╡ f2f3e1fb-0737-45a7-94c2-787b23bc2267
ixevtdc = findall(!iszero, evtdc)


# ╔═╡ d77017b8-ce47-41d8-90c2-0cf36761382c
md"""
#### Find the times of events which are not zero (thus one)
"""

# ╔═╡ e4b5ac2c-0c18-417a-bb4e-71e8ea84ec77
dct = tdc[ixevtdc]

# ╔═╡ 3917a74a-1070-49ed-80d1-fa4f21695a0b
sgnt = tsgn[ixevtsgn]

# ╔═╡ c0637347-a5a6-4084-9d38-770cd0590eb3
histogram(sgnt)

# ╔═╡ de23b0b1-0cd4-4eca-8914-79b190e3b111
histogram(dct)

# ╔═╡ 0e6e6bfb-59c1-484a-a2a1-4243cfcce32b
begin
hsgn  = jwf.jwfsim.h1d(sgnt, 100, 0.0, 2.5)
ptsgn =jwf.jwfsim.plot_h1dx(hsgn; xs = "t (mus)", ys = "counts", 
                           xlim=true, xl = (0.0, 1.0),
						   ylim=true, yl=(1,200), markersize=2, fwl=true)
hdc  = jwf.jwfsim.h1d(dct, 100, 0.0, 2.5)
ptdc =plot(ptsgn, jwf.jwfsim.plot_h1dx(hdc; xs = "t (mus)", ys = "counts", 
                           xlim=true, xl = (0.0, 1.0),
						   ylim=true, yl=(1,200), markersize=2, fwl=true))

end

# ╔═╡ e83e2775-4a26-474d-b7d4-44d0829b2fe4
begin

end

# ╔═╡ a0abf20e-cbcd-43a1-9b7f-23993e838539
hsgn.centers

# ╔═╡ 873cc039-3349-4873-9b6c-17ddfb1d37be
pois1=Poisson(10.0)

# ╔═╡ 2c300fea-b7d8-4289-be2a-be1f1d845d16
dc = rand(pois1, 10000)

# ╔═╡ 6bd6d9df-cb79-4069-9e92-2b719d8adfdf
histogram(dc)

# ╔═╡ 1279174a-3e68-4dad-9b41-7c798807bd44
md"""
## Functions
``\sigma``
"""

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
# ╠═aadf254a-47ea-417f-9f5e-734486445e32
# ╠═abf3d6d5-ca26-434c-b248-4c2395c79451
# ╠═babadc24-94dd-4c00-80ae-42a0bbc5fcc7
# ╠═e02432fe-7b16-42ae-af67-9959ad68e1a9
# ╠═ec58050a-51bd-49fc-83f3-77980db52d7c
# ╠═42a0fe7f-2945-4804-ae7e-bb5763497853
# ╠═f69269ad-592c-4ea8-b228-1b97177299b4
# ╠═feaddf95-827c-417d-b6c0-59e86f839ae9
# ╠═157ebdfe-082f-475f-b4ba-1c022817d259
# ╠═ad2254eb-efe5-44b5-a1d1-c023a4456801
# ╠═22442f67-6769-4f83-88b6-060f5a204748
# ╠═31c8c858-2178-4fce-a772-d116dca3f496
# ╠═a3d8ba37-07dd-4424-9c6f-e4b453205b90
# ╠═0a2b607f-4e30-41be-a73d-c0daf5c8b855
# ╠═ebb200f1-f5ab-4b66-a7ce-c2173fd1fc8f
# ╠═5a69f27c-86b2-48bb-99b1-357fbee639bf
# ╠═e672f7bc-5020-40ae-9706-a12640607b6f
# ╠═6bb4cd1f-7e2c-4b1d-acd8-d344b59534d9
# ╠═8375a215-ce9c-4a38-96a5-861e0f452282
# ╠═1d8692e3-7cf4-4739-a750-3f713512979b
# ╠═aa9e5f13-28e3-424b-a6b6-0cfb0c1c7625
# ╠═05e4e447-4012-46e2-a6cf-e14c49e31d40
# ╠═31de8cca-d53b-4d23-892b-294e19e0daf3
# ╠═ea0ccd60-90e0-4d83-91d7-59baba4d0063
# ╠═b6339888-e34c-47e4-8af1-2e9393461433
# ╠═663f5495-7e50-4a58-b57f-db915a17fd15
# ╠═fc7160cf-f81e-4778-85e3-a27383f7b5c2
# ╠═45c79d1c-650e-4db7-8e9c-4d7cf4851fa9
# ╠═8499afe3-5313-44cb-920f-99e4c2b2f13c
# ╠═b0e163eb-9ffc-4dbd-b886-1b1ca8971ded
# ╠═cb9b6dfa-d54c-4635-9bff-1a89b73ad93a
# ╠═7e9efea1-061b-4d57-a6fb-b66a955d2926
# ╠═8cd2b3cb-660b-48a5-a0e7-402fb835b499
# ╠═ac8d75b1-4804-423b-8268-124ead6366eb
# ╠═9fda79b5-f02d-479a-a7e8-9e2e6ca3e8fc
# ╠═38947c5e-aa3f-427d-bee1-a4532693381e
# ╠═7aef4cca-a6ce-4f16-ba7d-e7914d054247
# ╠═61ffc4ea-5107-453a-8048-9bdc11a94659
# ╠═d3810185-890e-4c52-a947-e53d6d9e0ebc
# ╠═2e403946-4a25-45ee-935b-2f6af7f87658
# ╠═1ae83c20-d8ec-4717-9275-8828bc1e193a
# ╠═2d68c170-0d91-4d84-975a-0165b86a2883
# ╠═4f7a9fc1-9534-42fa-890d-f298fba635e7
# ╠═eedef1fc-4b5a-49b8-8f61-ae99cc27215e
# ╠═57da5efb-d21b-4819-a623-32ea162e6456
# ╠═ea7e2262-ca35-4537-b914-5e1971cbfc20
# ╠═3d3fa4a0-be5c-450f-9b0a-8500878f38b1
# ╠═0588ae25-2f71-48d1-9fbf-e682f6161f24
# ╠═e38ad10c-f0c2-4fd8-b98e-27ef6b7b591e
# ╠═16f91e51-2f6b-4431-9994-b471a1f3994e
# ╠═e7c14906-1994-4bfc-b9cd-2f9ccc8cd2ad
# ╠═540cebd9-7630-4f35-9bf6-091516eda6bb
# ╠═c90c9b5a-ff51-46b0-9ea1-3ec1c745d68a
# ╠═5a70603e-ce7a-4518-ad58-86c720825267
# ╠═f2f3e1fb-0737-45a7-94c2-787b23bc2267
# ╠═d77017b8-ce47-41d8-90c2-0cf36761382c
# ╠═e4b5ac2c-0c18-417a-bb4e-71e8ea84ec77
# ╠═3917a74a-1070-49ed-80d1-fa4f21695a0b
# ╠═c0637347-a5a6-4084-9d38-770cd0590eb3
# ╠═de23b0b1-0cd4-4eca-8914-79b190e3b111
# ╠═0e6e6bfb-59c1-484a-a2a1-4243cfcce32b
# ╠═e83e2775-4a26-474d-b7d4-44d0829b2fe4
# ╠═a0abf20e-cbcd-43a1-9b7f-23993e838539
# ╠═873cc039-3349-4873-9b6c-17ddfb1d37be
# ╠═2c300fea-b7d8-4289-be2a-be1f1d845d16
# ╠═6bd6d9df-cb79-4069-9e92-2b719d8adfdf
# ╠═1279174a-3e68-4dad-9b41-7c798807bd44
# ╠═277be488-9840-4f62-874f-d4f9fbb893d8
