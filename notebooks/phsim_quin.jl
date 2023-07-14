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
	using LsqFit
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

# ╔═╡ c6baa96a-2793-4f82-a932-be627f989964
molecules = ["ANN155", "Ruth", "IrQuin"]

# ╔═╡ e9dd3d52-e54a-4315-8964-0b4c9f770fea
states = ["free", "Ba"]

# ╔═╡ 812d1a54-72b5-4a7f-b9ed-afc91b342e0e
mpath = Dict("IrQuin"=>"IrQUIN_UVvis_1E-5M_molar_extinction_coefficient.csv",
			"Ruth"=>"ZFFL28_RUTH_AAN247_cross section.csv",
            "ANN155"=>"AAN155_1E-5M_molar_extinction_coefficient.csv")

# ╔═╡ 9e7a849e-cd01-478a-98ec-85af0c5303af
epath = Dict("ANN155"=>"AAN155_1E-5_ACN_exc370.csv")

# ╔═╡ 0b28bd82-2b12-41c4-802b-e456fb691f19
mpar = Dict("IrQuinfree"=>(λ1=1600.0ns, f1=0.92, Q=0.024),
	"IrQuinBa"=>(λ1=1600.0ns, f1=0.92, Q=0.05),
	"ANN155free"=>(λ1=3.4ns, f1=1.0, Q=0.043),
	"ANN155Ba"=>(λ1=364.0ns, f1=0.73, λ2=1011.0ns, f2=0.22, Q=0.11),
	"Ruthfree"=>(λ1=3.4ns, f1=1.0, Q=0.14),
	"RuthBa"=>(λ1=1094.0ns, f1=0.96, Q=0.13)
)

# ╔═╡ 2b3561b7-36b0-485b-8532-de1d89702423
md""" Select molecules : $(@bind xmol Select(molecules))"""

# ╔═╡ f9b30836-d841-452a-b23c-a4c256ce3235
md""" Select state : $(@bind xstate Select(states))"""

# ╔═╡ abf3d6d5-ca26-434c-b248-4c2395c79451
begin
	pdata = joinpath(ENV["JWfSim"], "data")
	mequin =mpath[xmol]
	exct = epath[xmol]
end

# ╔═╡ 38ded754-c916-44a4-b673-02fcb9055cf8
begin
	xmolf = string(xmol, "free")
	xmolb = string(xmol, "Ba")
	md"""
	- molecule: free = $(xmolf);  chelated = $(xmolb)
	"""
end

# ╔═╡ babadc24-94dd-4c00-80ae-42a0bbc5fcc7
begin
	mequindfr = sort(jwf.jwfsim.load_df_from_csv(pdata, mequin, jwf.jwfsim.spG))
	names(mequindfr)
	mequindf = select(mequindfr, [:lnm, :eM5, :eM5Ba])
	mequindf[!, "ϵ"] = mequindf[!, "eM5"] .* 1.0 * cm^-1 .* M^-1
	mequindf[!, "ϵba"] = mequindf[!, "eM5Ba"] .* 1.0 * cm^-1 .* M^-1
	mequindf[!, "λ"] = mequindf[!, "lnm"] .* 1.0 * nm
	mqdf = select(mequindf, Not([:lnm, :eM5, :eM5Ba]))
	istep = mqdf.λ[2] - mqdf.λ[1] 
	
	quinf = jwf.jwfsim.Fluorophore(mqdf.λ[1], mqdf.λ[end], mqdf.ϵ, mpar[xmolf].Q)
	
	quinb = jwf.jwfsim.Fluorophore(mqdf.λ[1], mqdf.λ[end], mqdf.ϵba, mpar[xmolb].Q)
	
	xquinf = jwf.jwfsim.xsecfl(quinf, istep/nm)
	xquinb = jwf.jwfsim.xsecfl(quinb, istep/nm)
	
	md"""
	- Cross section at 375 nm: 
	- Free = $(xquinf(375.0*nm))
	- Ba2+ = $(xquinb(375.0*nm))
	"""
end

# ╔═╡ e02432fe-7b16-42ae-af67-9959ad68e1a9
md"""
#### Absorption spectrum 
"""

# ╔═╡ ec58050a-51bd-49fc-83f3-77980db52d7c
begin
	lblf = string(xmolf,": σ x 10^{16}")
	pxsecf = plot(mqdf.λ, xquinf.(mqdf.λ)*1e+16, label= lblf, lw=2)
	xlims!(210.,750.)
	lblb = string(xmolb,": σ x 10^{16}")
	pxsecb = plot(mqdf.λ, xquinb.(mqdf.λ)*1e+16, label= lblb, lw=2)
	xlims!(210.,750.)
	plot(pxsecf, pxsecb)
end


# ╔═╡ 7dad7b06-d969-43bf-8a15-30f30a05a72b
md"""
#### Emission spectrum 
"""

# ╔═╡ a1e65010-0e0f-4e6e-917a-b643a92b9ec4
begin
	excdf = sort(jwf.jwfsim.load_df_from_csv(pdata, exct, jwf.jwfsim.spG))
end

# ╔═╡ 7048b1a0-bd71-4f29-be1d-b18cd7191ec1
begin
	lblar =string(xmolf,"Ar")
	lblarba =string(xmolb,"ArBa")
	pe375 = plot(excdf.lnm, excdf.e375, label= xmolf, lw=2)
	pe375ba = plot!(pe375,excdf.lnm, excdf.e375Ba, label= xmolb, lw=2)
	pe375ar = plot(excdf.lnm, excdf.e375Ar, label= lblar, lw=2)
	pe375arba = plot!(pe375ar,excdf.lnm, excdf.e375ArBa, label= lblarba, lw=2)
	plot(pe375ba, pe375arba)
end

# ╔═╡ 42a0fe7f-2945-4804-ae7e-bb5763497853
md"""
## Setup
"""

# ╔═╡ f69269ad-592c-4ea8-b228-1b97177299b4
md"""
### Define lasers
"""

# ╔═╡ 15ebb995-0b3b-45c0-bed9-b6ecf6d101f8
lasers = ["pq375l", "epl375"]

# ╔═╡ 960472cc-bdfa-40c4-b4d6-31bf8e840598
md""" Select laser : $(@bind xlsr Select(lasers))"""

# ╔═╡ ad2254eb-efe5-44b5-a1d1-c023a4456801
md" set laser repetition rate in MHz $(@bind ff NumberField(0.5:0.5:10.0^2; default=1.0))"

# ╔═╡ 31c8c858-2178-4fce-a772-d116dca3f496
md"""
### Define objectives
- Current objective is a LMM-40x-VUV-Vac
"""

# ╔═╡ 8d891f09-9745-4b83-8435-028f64f7ecb0
objs = ["MUE12900","LMM-40x-VUV-Vac"]

# ╔═╡ ae7c9b4b-38f0-4bd6-ba62-4ebcbf4b04fc
md""" Select objective : $(@bind xobj Select(objs))"""

# ╔═╡ be9dfc6d-c627-4a38-8b46-4d9b1310162c
if xobj == "LMM-40x-VUV-Vac"
	obj = jwf.jwfsim.Objective(0.5, 40.0, 5.0mm, 5.1mm, 7.8mm, 0.85)
	tobj = jwf.jwfsim.tLMM40VUV()
elseif xobj == "MUE12900"
	obj = jwf.jwfsim.Objective(0.9, 100.0)
	tobj = jwf.jwfsim.tmue12900()
end

# ╔═╡ 1eb59032-324e-4f36-8c6e-647f84e5c570
md"""
## Detected photons
"""

# ╔═╡ 45c79d1c-650e-4db7-8e9c-4d7cf4851fa9
md"""
- The number of photons detected will be
$n_d = n \times \epsilon \times T_i \times T_o \times T_s \times Q$

Where

- ``\epsilon`` is the geometrical efficiency (depends on NA)
- ``T_i`` is the input transmitance of the objective (at input ``\lambda``)
- ``T_o`` is the output transmitance of the objective (at output ``\lambda``)
- ``T_s`` is the optical transmitance (filters, coupling to fibers, etc.)
- Q is the (average) effiency of APD (or PMT)


"""

# ╔═╡ be76b61f-843e-4dbc-b13a-4aebb4cfba5f
md" Set optical transmitance of the system  $(@bind ts NumberField(0.1:0.1:1; default=0.5))"

# ╔═╡ 8499afe3-5313-44cb-920f-99e4c2b2f13c
begin
	ϵ = jwf.jwfsim.geometrical_transmission(obj)
	md"""
	Efficiency due to NA = $(round(ϵ, digits=2))
	"""
end

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

# ╔═╡ 61ffc4ea-5107-453a-8048-9bdc11a94659
md" set DC in Hz $(@bind dchz NumberField(1:0.5:10.0^2; default=100))"

# ╔═╡ 1ae83c20-d8ec-4717-9275-8828bc1e193a
md"""
## Distribution per pulse of signal and DC

- Both the signal and the background follow a Poisson distribution, with λ equal to the mean number of signal (or background) events
"""

# ╔═╡ eedef1fc-4b5a-49b8-8f61-ae99cc27215e
md"""
### Time distribution of signal and DC

- The time distribution of the signal is an exponential with the lifetime characteristic of the molecule. The distribution of the background is flat. 
"""

# ╔═╡ b102fe6d-0654-4d7b-afc7-56f3c0f1e883
md" Set Lifetime in μs $(@bind xl NumberField(0.3:0.1:5.0; default=1.0))"

# ╔═╡ 57da5efb-d21b-4819-a623-32ea162e6456
dexpmus = Exponential(xl)

# ╔═╡ ea7e2262-ca35-4537-b914-5e1971cbfc20
md"""
### Generate the events and the times of signal and DC
"""

# ╔═╡ 3d3fa4a0-be5c-450f-9b0a-8500878f38b1
md" set acquistion time (in seconds) $(@bind act NumberField(1:0.5:10.0^2; default=1))"

# ╔═╡ d77017b8-ce47-41d8-90c2-0cf36761382c
md"""
#### Find the times of events which are not zero (thus one)
"""

# ╔═╡ 1279174a-3e68-4dad-9b41-7c798807bd44
md"""
## Functions
"""

# ╔═╡ 5c3ac9c2-8b16-4072-ba9a-1d336cea52bf
function float_to_fstr(val, fmt)
	ex = quote
	@sprintf $fmt $val
	end
	eval(ex)
end

# ╔═╡ b1d523c2-83cc-4893-b1a7-c62d3fd40a9d
function tfit(xdata::Vector{Float64}, ydata::Vector{Float64}, 
	          i0::Integer, l0::Integer; pa0::Vector{Float64}=[100.0, 0.5])
	ffun(t, N1, λ) = N1 * exp(-t/λ)  
	mffun(t, p) = p[1] * exp.(-t/p[2]) 

	fit = curve_fit(mffun, xdata[i0:l0], ydata[i0:l0], pa0)
	coef(fit), stderror(fit), ffun.(xdata, coef(fit)...)
end

# ╔═╡ e5719716-65fc-4166-8d54-19bbb7dc025f
function sgndc(dsgn, ddc, dexpmus, ntx)
	evtsgn = rand(dsgn, Int(ntx))
	evtdc = rand(ddc, Int(ntx))
	tsgn = rand(dexpmus, Int(ntx))
	tdc = rand(Int(ntx))
	evtsgn, evtdc, tsgn, tdc 
end

# ╔═╡ 277be488-9840-4f62-874f-d4f9fbb893d8
function epl375(f=1.0MHz)
    λ = 370.0nm
    fr = 20.0MHz
    pr = 150.0μW 
    P = pr/(fr/f)
    tau = 75.0ps
    jwf.jwfsim.PLaser(λ, P, f, tau)
end

# ╔═╡ cbe94ca2-0b1f-4246-8843-6b8e4e5e55cc
function pqLDH375B(f=1.0MHz)
    λ = 376.0nm
    fr = 80.0MHz
    pr = 7.4mW 
    P = pr/(fr/f)
    tau = 56.0ps
    jwf.jwfsim.PLaser(λ, P, f, tau)
end

# ╔═╡ 24c43b9a-cff6-462f-9d13-95e3bdb71792
if xlsr == "epl375"
	lsr = epl375(ff*MHz)
elseif xlsr == "pq375l"
	lsr = pqLDH375B(ff*MHz)
end

# ╔═╡ e672f7bc-5020-40ae-9706-a12640607b6f
begin
	w0 = jwf.jwfsim.w0f(lsr.λ, obj.d, obj.f)
	glsr = jwf.jwfsim.GaussianLaser(lsr, w0)
	Igl = jwf.jwfsim.I(glsr)
	ρ = uconvert(kW/cm^2, Igl(0.0nm,0.0nm))
end

# ╔═╡ 0a2b607f-4e30-41be-a73d-c0daf5c8b855
md"""
### Photon rate per molecule 
- Gaussian beam filling the objective.
- Objective diameter = $(obj.d)
- Objective focal = $(obj.f)
- Waist of focused beam = $(w0)
- Power density at focus =$ρ
- Photon density at focus =$(glsr.ρ0)
"""

# ╔═╡ aa9e5f13-28e3-424b-a6b6-0cfb0c1c7625
begin
	λl = 375.0nm
	ng = glsr.ρ0 * xquinf(λl)
	ngp = uconvert(NoUnits, ng/lsr.f)
	ngba = glsr.ρ0 * xquinb(λl)
	ngbap = uconvert(NoUnits, ngba/lsr.f)
end

# ╔═╡ 05e4e447-4012-46e2-a6cf-e14c49e31d40
md"""
	##### Assuming a gaussian beam that fills the objective
	- dl = $(jwf.jwfsim.spot_size(glsr))
	- Free molecule
	- Photon rate = $ng
	- Photon rate per pulse = $ngp
	- Ba2+ molecule
	- Photon rate = $ngba
	- Photon rate per pulse = $ngbap
	"""

# ╔═╡ a57f3c34-3319-4275-bdda-982ff5f67445
begin
	ngtoapd = ϵ * ng * tobj(375.0nm) *  tobj(600.0nm) * ts * eapd.(600.0nm) 
	ngpx = uconvert(NoUnits, ngtoapd/lsr.f)
	ngbatoapd = ϵ * ngba * tobj(375.0nm) *  tobj(600.0nm) * ts * eapd.(600.0nm) 
	ngbapx = uconvert(NoUnits, ngbatoapd/lsr.f)
end

# ╔═╡ 7aef4cca-a6ce-4f16-ba7d-e7914d054247
md"""
- **Free**
- Number of photons per molecule per unit time =$ng
- Transmitted through the objective =$(ϵ * ng * tobj(375.0nm) *  tobj(600.0nm))
- Arriving to the APD =$(ϵ * ng * tobj(375.0nm) *  tobj(600.0nm) * ts)
- Transmitted through the apd =$ngtoapd
- Detected, per pulse =$ngpx
- **Ba2**
- Number of photons per molecule per unit time =$ngba
- Transmitted through the objective =$(ϵ * ngba * tobj(375.0nm) *  tobj(600.0nm))
- Arriving to the APD =$(ϵ * ngba * tobj(375.0nm) *  tobj(600.0nm) * ts)
- Transmitted through the apd =$ngbatoapd
- Detected, per pulse =$ngbapx
"""

# ╔═╡ 2d68c170-0d91-4d84-975a-0165b86a2883
dsgn = Poisson(ngpx)

# ╔═╡ d3810185-890e-4c52-a947-e53d6d9e0ebc
ngdcx = uconvert(NoUnits, dchz * Hz/lsr.f)

# ╔═╡ 2e403946-4a25-45ee-935b-2f6af7f87658
md"""
- DC (Hz) =$dchz
- DC (per pulse) =$(round(ngdcx, digits=5))
"""

# ╔═╡ 4f7a9fc1-9534-42fa-890d-f298fba635e7
ddc = Poisson(ngdcx)

# ╔═╡ 0588ae25-2f71-48d1-9fbf-e682f6161f24
ntx = uconvert(NoUnits, act*s*lsr.f)

# ╔═╡ e38ad10c-f0c2-4fd8-b98e-27ef6b7b591e
evtsgn, evtdc, tsgn, tdc = sgndc(dsgn, ddc, dexpmus, ntx)

# ╔═╡ 8e185d65-607c-435c-9f1b-da738060ddd1
function tsgndc(evtsgn, evtdc)
	ixevtsgn = findall(!iszero, evtsgn)
	ixevtdc = findall(!iszero, evtdc)
	sgnt = tsgn[ixevtsgn]
	dct = tdc[ixevtdc]
	sgndct = copy(sgnt)
	append!(sgndct,dct)
	sgnt, dct, sgndct
end

# ╔═╡ 94e22133-5d39-4b53-abe6-87ff249b825b
sgnt, dct, sgndct = tsgndc(evtsgn, evtdc)

# ╔═╡ 0e6e6bfb-59c1-484a-a2a1-4243cfcce32b
begin
hsgn  = jwf.jwfsim.h1d(sgnt, 100, 0.0, 2.5)
ptsgn =jwf.jwfsim.plot_h1dx(hsgn; xs = "t (mus)", ys = "counts", 
                           xlim=true, xl = (0.0, 1.0),
						   ylim=false, yl=(1,400), markersize=2, fwl=true)
hdc  = jwf.jwfsim.h1d(dct, 100, 0.0, 2.5)
ptdc =plot(ptsgn, jwf.jwfsim.plot_h1dx(hdc; xs = "t (mus)", ys = "counts", 
                           xlim=true, xl = (0.0, 1.0),
						   ylim=true, yl=(1,200), markersize=2, fwl=true))

end

# ╔═╡ b8550b70-67c0-49d2-9900-5d262bbf397a
begin
	ptsgndc =jwf.jwfsim.plot_h1dx(hsgn, hdc; xs = "t (mus)", ys = "counts", 
                           xlim=true, xl = (0.0, 1.0), 
                           label1="signal", label2="dc", legend=true,
						   ylim=false, yl=(1,400), markersize=2, fwl=true)
	xdata = hsgn.centers
	ydata = hsgn.weights
	cft, stdft, yft = tfit(xdata, ydata, 1, 100; pa0=[100.0, 0.5])
	fts = float_to_fstr(cft[2], "%7.1f")
	lbl = string("fitted lifetime (mus) =", fts)
	plot(ptsgndc, xdata, yft, label=lbl)
end


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
# ╠═c6baa96a-2793-4f82-a932-be627f989964
# ╠═e9dd3d52-e54a-4315-8964-0b4c9f770fea
# ╠═812d1a54-72b5-4a7f-b9ed-afc91b342e0e
# ╠═9e7a849e-cd01-478a-98ec-85af0c5303af
# ╠═0b28bd82-2b12-41c4-802b-e456fb691f19
# ╠═2b3561b7-36b0-485b-8532-de1d89702423
# ╠═f9b30836-d841-452a-b23c-a4c256ce3235
# ╠═abf3d6d5-ca26-434c-b248-4c2395c79451
# ╠═38ded754-c916-44a4-b673-02fcb9055cf8
# ╠═babadc24-94dd-4c00-80ae-42a0bbc5fcc7
# ╠═e02432fe-7b16-42ae-af67-9959ad68e1a9
# ╠═ec58050a-51bd-49fc-83f3-77980db52d7c
# ╠═7dad7b06-d969-43bf-8a15-30f30a05a72b
# ╠═a1e65010-0e0f-4e6e-917a-b643a92b9ec4
# ╠═7048b1a0-bd71-4f29-be1d-b18cd7191ec1
# ╠═42a0fe7f-2945-4804-ae7e-bb5763497853
# ╠═f69269ad-592c-4ea8-b228-1b97177299b4
# ╠═15ebb995-0b3b-45c0-bed9-b6ecf6d101f8
# ╠═960472cc-bdfa-40c4-b4d6-31bf8e840598
# ╠═ad2254eb-efe5-44b5-a1d1-c023a4456801
# ╠═24c43b9a-cff6-462f-9d13-95e3bdb71792
# ╠═31c8c858-2178-4fce-a772-d116dca3f496
# ╠═8d891f09-9745-4b83-8435-028f64f7ecb0
# ╠═ae7c9b4b-38f0-4bd6-ba62-4ebcbf4b04fc
# ╠═be9dfc6d-c627-4a38-8b46-4d9b1310162c
# ╠═0a2b607f-4e30-41be-a73d-c0daf5c8b855
# ╠═e672f7bc-5020-40ae-9706-a12640607b6f
# ╠═aa9e5f13-28e3-424b-a6b6-0cfb0c1c7625
# ╠═05e4e447-4012-46e2-a6cf-e14c49e31d40
# ╠═1eb59032-324e-4f36-8c6e-647f84e5c570
# ╠═45c79d1c-650e-4db7-8e9c-4d7cf4851fa9
# ╠═be76b61f-843e-4dbc-b13a-4aebb4cfba5f
# ╠═8499afe3-5313-44cb-920f-99e4c2b2f13c
# ╠═cb9b6dfa-d54c-4635-9bff-1a89b73ad93a
# ╠═7e9efea1-061b-4d57-a6fb-b66a955d2926
# ╠═8cd2b3cb-660b-48a5-a0e7-402fb835b499
# ╠═a57f3c34-3319-4275-bdda-982ff5f67445
# ╠═7aef4cca-a6ce-4f16-ba7d-e7914d054247
# ╠═61ffc4ea-5107-453a-8048-9bdc11a94659
# ╠═d3810185-890e-4c52-a947-e53d6d9e0ebc
# ╠═2e403946-4a25-45ee-935b-2f6af7f87658
# ╠═1ae83c20-d8ec-4717-9275-8828bc1e193a
# ╠═2d68c170-0d91-4d84-975a-0165b86a2883
# ╠═4f7a9fc1-9534-42fa-890d-f298fba635e7
# ╠═eedef1fc-4b5a-49b8-8f61-ae99cc27215e
# ╠═b102fe6d-0654-4d7b-afc7-56f3c0f1e883
# ╠═57da5efb-d21b-4819-a623-32ea162e6456
# ╠═ea7e2262-ca35-4537-b914-5e1971cbfc20
# ╠═3d3fa4a0-be5c-450f-9b0a-8500878f38b1
# ╠═0588ae25-2f71-48d1-9fbf-e682f6161f24
# ╠═e38ad10c-f0c2-4fd8-b98e-27ef6b7b591e
# ╠═d77017b8-ce47-41d8-90c2-0cf36761382c
# ╠═94e22133-5d39-4b53-abe6-87ff249b825b
# ╠═0e6e6bfb-59c1-484a-a2a1-4243cfcce32b
# ╠═b8550b70-67c0-49d2-9900-5d262bbf397a
# ╠═1279174a-3e68-4dad-9b41-7c798807bd44
# ╠═5c3ac9c2-8b16-4072-ba9a-1d336cea52bf
# ╠═b1d523c2-83cc-4893-b1a7-c62d3fd40a9d
# ╠═e5719716-65fc-4166-8d54-19bbb7dc025f
# ╠═8e185d65-607c-435c-9f1b-da738060ddd1
# ╠═277be488-9840-4f62-874f-d4f9fbb893d8
# ╠═cbe94ca2-0b1f-4246-8843-6b8e4e5e55cc
