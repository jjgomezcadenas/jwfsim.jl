### A Pluto.jl notebook ###
# v0.19.27

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
	using EasyFit
	using LaTeXStrings
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
PlutoUI.TableOfContents(title="ANN-155 simulation", indent=true)

# ╔═╡ 78ff48a5-a903-44a3-af8d-ffc4c832cea5
md"""
# Simulation steps
"""

# ╔═╡ 6d047015-90c5-4774-9fa6-3e5f49160c19
md"""
## Absorption and excitation spectra for ANN-155 and cousins 
"""

# ╔═╡ aadf254a-47ea-417f-9f5e-734486445e32
md"""
### Define fluorophores and cross sections
"""

# ╔═╡ c6baa96a-2793-4f82-a932-be627f989964
begin
	molecules = ["ANN155", "Ruth", "IrQuin"]
	mpath = Dict("IrQuin"=>"IrQUIN_UVvis_1E-5M_molar_extinction_coefficient.csv",
			"Ruth"=>"ZFFL28_RUTH_AAN247_cross section.csv",
            "ANN155"=>"AAN155_1E-5M_molar_extinction_coefficient.csv")
	epath = Dict("ANN155"=>"AAN155_1E-5_ACN_exc370.csv")
	mpar = Dict("IrQuinfree"=>(λ1=1600.0ns, f1=0.92, Q=0.024),
	"IrQuinBa"=>(λ1=1600.0ns, f1=0.92, Q=0.05),
	"ANN155free"=>(λ1=280.0ns, f1=0.61, λ2=538.0ns, f2=0.39, Q=0.043),
	"ANN155Ba"=>(λ1=319.0ns, f1=0.34, λ2=1002.0ns, f2=0.66, Q=0.11),
	"Ruthfree"=>(λ1=3.4ns, f1=1.0, Q=0.14),
	"RuthBa"=>(λ1=1094.0ns, f1=0.96, Q=0.13)
)
	nothing
end

# ╔═╡ 2b3561b7-36b0-485b-8532-de1d89702423
md""" Select molecules : $(@bind xmol Select(molecules))"""

# ╔═╡ abf3d6d5-ca26-434c-b248-4c2395c79451
begin
	pdata = joinpath(ENV["JWfSim"], "data")
	mequin =mpath[xmol]
	exct = epath[xmol]
	nothing
end

# ╔═╡ 38ded754-c916-44a4-b673-02fcb9055cf8
begin
	xmolf = string(xmol, "free")
	xmolb = string(xmol, "Ba")
	md"""
	#### Selected molecules
	- molecule: free = $(xmolf);  chelated = $(xmolb)
	"""
end

# ╔═╡ a276091d-0c83-4a75-be20-51ed938d7415
md"""
#### Lifetimes

- Free: λ1 = $(round(mpar[xmolf].λ1/ns)) ns,  λ2 = $(round(mpar[xmolf].λ2/ns)) ns, N1 = $(round(mpar[xmolf].f1, digits=3)), N2 = $(round(mpar[xmolf].f2, digits=3)) 

- Ba2+: λ1 = $(round(mpar[xmolb].λ1/ns)) ns,  λ2 = $(round(mpar[xmolb].λ2/ns)) ns, N1 = $(round(mpar[xmolb].f1, digits=3)), N2 = $(round(mpar[xmolb].f2, digits=3)) 



"""

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
	#### Cross section at 375 nm: 
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
	nothing
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
lasers = ["pq375l", "epl375"];

# ╔═╡ 960472cc-bdfa-40c4-b4d6-31bf8e840598
md""" Select laser : $(@bind xlsr Select(lasers))"""

# ╔═╡ ad2254eb-efe5-44b5-a1d1-c023a4456801
md" set laser repetition rate in MHz $(@bind ff NumberField(0.1:0.1:10.0; default=0.2))"

# ╔═╡ 31c8c858-2178-4fce-a772-d116dca3f496
md"""
### Define objectives
"""

# ╔═╡ 8d891f09-9745-4b83-8435-028f64f7ecb0
objs = ["LMM-40x-VUV-Vac", "MUE12900"];

# ╔═╡ ae7c9b4b-38f0-4bd6-ba62-4ebcbf4b04fc
md""" Select objective : $(@bind xobj Select(objs))"""

# ╔═╡ be9dfc6d-c627-4a38-8b46-4d9b1310162c
begin
	if xobj == "LMM-40x-VUV-Vac"
		obj = jwf.jwfsim.Objective(0.5, 40.0, 5.0mm, 5.1mm, 7.8mm, 0.85)
		tobj = jwf.jwfsim.tLMM40VUV()
	elseif xobj == "MUE12900"
		obj = jwf.jwfsim.Objective(0.9, 100.0)
		tobj = jwf.jwfsim.tmue12900()
	end
	ϵ = jwf.jwfsim.geometrical_transmission(obj)
	md"""
	- Optical efficiency due to NA of objective = $(round(ϵ, digits=2))
	- Tranmission of UV light = $(round(tobj(375.0nm), digits=2))
	- Tranmission of red light = $(round(tobj(600.0nm), digits=2))
	"""
end

# ╔═╡ 1eb59032-324e-4f36-8c6e-647f84e5c570
md"""
### Detected photon rate
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

# ╔═╡ cb9b6dfa-d54c-4635-9bff-1a89b73ad93a
md"""
- APD efficiency
"""

# ╔═╡ 7e9efea1-061b-4d57-a6fb-b66a955d2926
eapd = jwf.jwfsim.ϵapd();

# ╔═╡ 8cd2b3cb-660b-48a5-a0e7-402fb835b499
begin
	wl = 350.0:10.0:1000.0
	wλ=[w*nm for w in wl]
	plot(wl, eapd.(wλ), label= "APD efficiency", lw=2)
	#xlims!(200.,550.)

end

# ╔═╡ 61ffc4ea-5107-453a-8048-9bdc11a94659
md" set DC in Hz $(@bind dchz NumberField(1:0.5:10.0^2; default=100))"

# ╔═╡ 1a36aa14-a7ea-4db0-8ce6-4f52e0c1bd04
md" - set number of bkgnd molecules $(@bind nbkg NumberField(1:10^4; default=1))"

# ╔═╡ 4ce0fdf3-dc93-456d-8229-c05b0ef54268
md"""
## Simulation ingredients
"""

# ╔═╡ 3d3fa4a0-be5c-450f-9b0a-8500878f38b1
md" - set acquistion time (in seconds) $(@bind act NumberField(1:0.5:10^2; default=10.0))"

# ╔═╡ 1ae83c20-d8ec-4717-9275-8828bc1e193a
md"""
### Distribution per pulse of signal and DC

- Both the signal and the background follow a Poisson distribution, with λ equal to the mean number of signal (or background) events. The DC follows a flat (uniform) distribution 
"""

# ╔═╡ dd347540-1ee9-40e9-8f69-8ecda506cbbf


# ╔═╡ eedef1fc-4b5a-49b8-8f61-ae99cc27215e
md"""
### Time distribution of signal (free and chelated molecules) and DC

- The time distribution of the signal (free and chelated) is an exponential with the lifetime characteristic of the molecule. The distribution of the background is flat. 
"""

# ╔═╡ 70a95466-f7b4-4e5d-a0bf-b3c7a3150dc6
md"""
- Define functions **expf** and **expb** describing the time distribution of signals (f stands for free, b for barium), as well as models **mexpf** and **mexpb** which allow to throw random numbers according to the distributions specified for each one
"""

# ╔═╡ ea7e2262-ca35-4537-b914-5e1971cbfc20
md"""
### Generate the times of signal and DC
"""

# ╔═╡ cef81cd5-447d-459e-be70-4a7c531ed504
#md"""
#### Nominal time distributions for free and Ba2+
#"""

# ╔═╡ a23a9fbe-772a-4cc0-81de-f3b5dc58d766
md"""
#### Simulated response for free, Ba2+ and DC
"""

# ╔═╡ d6d35e53-a712-4e89-9e70-8c166f6b7550
md"""
#### Fit the reference distributions with a double exponential
"""

# ╔═╡ 9be5bbaa-32c4-41fb-b917-ed747dbeb213
md"""
## Extracting signal
"""

# ╔═╡ bf5af920-6a40-42af-8697-8fb089976a9f
md"""
### Integrate reference distributions for free and chelated molecules
"""

# ╔═╡ 1f19dfd5-cf6a-48a8-a65d-0fbc752bae99
md"Select background molecules $(@bind num_bkg Slider(1:10^4))"

# ╔═╡ 2e4620c0-9abb-4359-987b-839ebeccef96
md"""
##### Results
"""

# ╔═╡ de748626-2bad-4637-b948-bd2a2ad72e97


# ╔═╡ 0fd69074-e66b-40dc-8810-5821bd1c6ad8
md"""
- Generate data
"""

# ╔═╡ 1279174a-3e68-4dad-9b41-7c798807bd44
md"""
## Functions
"""

# ╔═╡ a8bf7b5f-286e-4bb7-9dc5-91d01878ed05
function fsgn(mpar, xmolf, xmolb)
	expf(t) = mpar[xmolf].f1 * exp(-t/mpar[xmolf].λ1) + mpar[xmolf].f2 * exp(-t/mpar[xmolf].λ2)
	
	expb(t) = mpar[xmolb].f1 * exp(-t/mpar[xmolb].λ1) + mpar[xmolb].f2 * exp(-t/mpar[xmolb].λ2)
	
	expf, expb	
end

# ╔═╡ b2bac877-7098-43d0-8241-c01b6e219bd1
function mfsgn(mpar, xmolf, xmolb)
	function combine_model(λ1, λ2, p)
		expf1 = Exponential(λ1)
		expf2 = Exponential(λ2)
		expf = MixtureModel([expf1,expf2],[p,1-p])
	end

	mexpf = combine_model(mpar[xmolf].λ1/ns, mpar[xmolf].λ2/ns, mpar[xmolf].f1)
	mexpb = combine_model(mpar[xmolb].λ1/ns, mpar[xmolb].λ2/ns, mpar[xmolb].f1)
	mexpf, mexpb	
end

# ╔═╡ 299ee6bd-b011-4937-a46c-c8febfcf61c7
begin
	expf, expb = fsgn(mpar, xmolf, xmolb)
	mexpf, mexpb = mfsgn(mpar, xmolf, xmolb)
	nothing
end

# ╔═╡ 7f6b34e2-ade3-4a0a-928b-66bee72456ab
function tevents(evtsgn, evtsgnba, evtdc, mexpf, mexpb, lsr)
	tsgn = rand(mexpf, evtsgn)
	tsgnba = rand(mexpb, evtsgnba)
	tdc = rand(evtdc) * uconvert(ns, 1/lsr.f)/ns
	tsgn, tsgnba, tdc 
end

# ╔═╡ c51c7590-35aa-4ba2-8539-61b42dfab5b2
function nevents(dsgn, dsgnba, ddc, ntx::Int)
	function tsgndc(evtsgn, evtsgnba, evtdc)
		ixevtsgn = findall(!iszero, evtsgn)
		ixevtsgnba = findall(!iszero, evtsgnba)
		ixevtdc = findall(!iszero, evtdc)
		length(ixevtsgn), length(ixevtsgnba), length(ixevtdc)
	end
	
	evtsgn = rand(dsgn, ntx)
	evtsgnba = rand(dsgnba, ntx)
	evtdc = rand(ddc, ntx)
	tsgndc(evtsgn, evtsgnba, evtdc)
end

# ╔═╡ 8e185d65-607c-435c-9f1b-da738060ddd1


# ╔═╡ b1d523c2-83cc-4893-b1a7-c62d3fd40a9d
function tfit(xdata::Vector{Float64}, ydata::Vector{Float64}; pa0::Vector{Float64}=[100.0, 1])
	
	ffun(t, N1, m) = N1 * (expf(t*ns) + m * expb(t*ns))	  
	mffun(t, p) = p[1] * (expf.(t*ns) + p[2] * expb.(t*ns) ) 

	fit = curve_fit(mffun, xdata, ydata, pa0)
	coef(fit), stderror(fit), ffun.(xdata, coef(fit)...)
end

# ╔═╡ b5345760-ae72-4125-a502-39318a87f5a3
function lbls(sgnfit, sgnfitba)
	function getlbls(λ1, λ2, sn1, sn2)
		sexp1 = round(λ1, digits=0)
		sexp2 = round(λ2, digits=0)
		sp1 = round(sn1 / (sn1+sn2), sigdigits=2)
		sp2 = round(1 -sp1, sigdigits=2)
		sexp1, sexp2, sp1, sp2
	end
	lfexp1, lfexp2, lfp1, lfp2 = getlbls(sgnfit.b[1], sgnfit.b[2], sgnfit.a[1], 
		                                 sgnfit.a[2])

	lbexp1, lbexp2, lbp1, lbp2 = getlbls(sgnfitba.b[1], sgnfitba.b[2], sgnfitba.a[1], 
		                                 sgnfitba.a[2])
		
	lbl1a = L"\lambda = %$lfexp1 \, \mu s \, (%$lfp1)" 
	lbl1b = L"\lambda2 = %$lfexp2 \, \mu s \, (%$lfp2)"
	lbl1 = string(lbl1a, "\n", lbl1b)
	lbl2a = L"\lambda1 = %$lbexp1 \, \mu s \,(%$lbp1)"
	lbl2b = L"\lambda2 = %$lbexp2 \, \mu s \, (%$lbp2)"
	string(lbl1a, "\n", lbl1b), string(lbl2a, "\n", lbl2b)
end

# ╔═╡ 5c3ac9c2-8b16-4072-ba9a-1d336cea52bf
function float_to_fstr(val, fmt)
	ex = quote
	@sprintf $fmt $val
	end
	eval(ex)
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
begin
	if xlsr == "epl375"
		lsr = epl375(ff*MHz)
	elseif xlsr == "pq375l"
		lsr = pqLDH375B(ff*MHz)
	end
	md"""
	- Laser wavelength = $(lsr.λ)
	- Laser frequency  = $(lsr.f) : time range = $(uconvert(μs, 1/lsr.f))
	- Laser power      = $(round(uconvert(μW, lsr.P)/μW)) μW
	"""
end

# ╔═╡ e672f7bc-5020-40ae-9706-a12640607b6f
begin
	w0 = jwf.jwfsim.w0f(lsr.λ, obj.d, obj.f)
	glsr = jwf.jwfsim.GaussianLaser(lsr, w0)
	Igl = jwf.jwfsim.I(glsr)
	ρ = uconvert(kW/cm^2, Igl(0.0nm,0.0nm))
	nothing
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
	nothing
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
	eff = ϵ  * tobj(375.0nm) *  tobj(600.0nm) * ts * eapd.(600.0nm)
	ngtoapd = eff * ng  
	ngpx = uconvert(NoUnits, ngtoapd/lsr.f)
	ngbatoapd = eff * ngba  
	ngbapx = uconvert(NoUnits, ngbatoapd/lsr.f)
	nothing
end

# ╔═╡ 7aef4cca-a6ce-4f16-ba7d-e7914d054247
md"""
- Overall optical efficiency of the system =$(round(eff, digits=3))
- **Free molecules**
- Number of photons per molecule per unit time =$ng
- Transmitted through the objective =$(ϵ * ng * tobj(375.0nm) *  tobj(600.0nm))
- Arriving to the APD =$(ϵ * ng * tobj(375.0nm) *  tobj(600.0nm) * ts)
- Transmitted through the apd =$ngtoapd
- Detected, per pulse =$ngpx
- **Ba2 molecules**
- Number of photons per molecule per unit time =$ngba
- Transmitted through the objective =$(ϵ * ngba * tobj(375.0nm) *  tobj(600.0nm))
- Arriving to the APD =$(ϵ * ngba * tobj(375.0nm) *  tobj(600.0nm) * ts)
- Transmitted through the apd =$ngbatoapd
- Detected, per pulse =$ngbapx
"""

# ╔═╡ d3810185-890e-4c52-a947-e53d6d9e0ebc
ngdcx = uconvert(NoUnits, dchz * Hz/lsr.f);

# ╔═╡ 2e403946-4a25-45ee-935b-2f6af7f87658
md"""
- DC (Hz) =$dchz
- DC (per pulse) =$(round(ngdcx, digits=5))
"""

# ╔═╡ 0588ae25-2f71-48d1-9fbf-e682f6161f24
begin
	ntx = uconvert(NoUnits, act*s*lsr.f)
	md"""
	- Number of pulses during acquisition time =$ntx
	"""
end

# ╔═╡ 2d68c170-0d91-4d84-975a-0165b86a2883
begin
	dsgn = Poisson(ngpx)
	dsgnba = Poisson(ngbapx)
	ddc = Poisson(ngdcx)
	dfree = Poisson(ngpx * num_bkg) # distribution of number of free molecules
	# evtsgn, evtsgnba, evtdc are references (correspond to 1 molecule free + 1 Ba2)
	evtsgn, evtsgnba, evtdc = nevents(dsgn, dsgnba, ddc, Int(ntx))
	# evtfr is reference (correspond to nbkg free molecule)
	evtfr, _, _ = nevents(dfree, dsgnba, ddc, Int(ntx))
	# dtfr, dtba, dtdc2 are measurements (correspond to nbkg free molecule)
	dtfr, dttba, dtdc = nevents(dfree, dsgnba, ddc, Int(ntx))
	md"""
	#### Number of detected photons per molecule during run:
	##### Reference: number of Ba2+ molecules = 1, number of Free = $nbkg
	- free = $evtsgn
	- Ba2+ = $evtsgnba
	- DC   = $evtdc

	##### Data: number of Ba2+ molecules = 1, number of Free = $num_bkg
	- free = $dtfr, substraction error: = $(round(sqrt(dtfr)))
	- Ba2+ = $dttba
	- DC   = $dtdc
	"""
end

# ╔═╡ 0822b115-575e-4915-87e9-1264abf6813a
begin
	# reference 1:1
	tsgn, tsgnba, tdc  =tevents(evtsgn, evtsgnba, evtdc, mexpf, mexpb,  lsr)
	# reference 1:nbkg
	tfree, _, _  =tevents(evtfr, evtsgnba, evtdc, mexpf, mexpb,  lsr)
	# data
	tfr, tba, tdc2  =tevents(dtfr, dttba, dtdc, mexpf, mexpb,  lsr)
	tdata = append!(tfr, tba, tdc2)
end

# ╔═╡ 12d7b082-0a19-495c-b5db-feb4ab700267
begin
	tlns =(uconvert(ns, 1/lsr.f))/ns
	xtt = collect(0:tlns) *ns
	#pexpf = plot(xtt, expf.(xtt))
	#pexpba = plot(xtt, expb.(xtt))
	#plot(pexpf, pexpba)
	nothing
end

# ╔═╡ 0e6e6bfb-59c1-484a-a2a1-4243cfcce32b
begin
	hsgn  = jwf.jwfsim.h1d(tsgn, 100, 0.0, tlns)
	ptsgn =jwf.jwfsim.plot_h1dx(hsgn; xs = "t (mus)", ys = "counts", 
	                           xlim=true, xl = (0.0, tlns), label="Free",legend=true,
							   ylim=false, yl=(1,400), markersize=2, fwl=false)
	hfree  = jwf.jwfsim.h1d(tfree, 100, 0.0, tlns)
	ptfree =jwf.jwfsim.plot_h1dx(hfree; xs = "t (mus)", ys = "counts", 
	                           xlim=true, xl = (0.0, tlns), 
	                           label="Free x nbkg",legend=true,
							   ylim=false, yl=(1,400), markersize=2, fwl=false)
	hsgnba  = jwf.jwfsim.h1d(tsgnba, 100, 0.0, tlns)
	ptsgnba =jwf.jwfsim.plot_h1dx(hsgnba; xs = "t (mus)", ys = "counts", 
	                           xlim=true, xl = (0.0, tlns),label="Ba2+",legend=true,
							   ylim=false, yl=(1,400), markersize=2, fwl=false)
	hdc  = jwf.jwfsim.h1d(tdc, 100, 0.0, tlns)
	ptdc = jwf.jwfsim.plot_h1dx(hdc; xs = "t (mus)", ys = "counts", 
                           xlim=true, xl = (0.0, tlns),
						   ylim=false, yl=(1,200), markersize=2, fwl=true)
	hdata  = jwf.jwfsim.h1d(tdata, 100, 0.0, tlns)
	pxdata =jwf.jwfsim.plot_h1dx(hdata; xs = "t (mus)", ys = "counts", 
	                           xlim=true, xl = (0.0, tlns), label="data",legend=true,
							   ylim=false, yl=(1,400), markersize=2, fwl=false)
plot(ptsgn, ptsgnba, ptfree, pxdata)

end

# ╔═╡ 4012d326-3a8e-426e-899f-5dddf7674b35
begin
	vdcmean =mean(hdc.weights)
	txsgn = hsgn.centers
	#vsgn = hsgn.weights .- vdcmean
	#vsgnba = hsgnba.weights .- vdcmean;
	vsgn = hsgn.weights 
	vsgnba = hsgnba.weights 
	sgnfitba = fitexp(txsgn, vsgnba,n=2)
	sgnfit = fitexp(txsgn, vsgn,n=2)
	lbl1, lbl2 = lbls(sgnfit, sgnfitba)
	nothing
end

# ╔═╡ 903ac25f-493f-4b70-87b1-480b8d8396fb
sgnfit

# ╔═╡ d2ed2a93-e059-4900-b0dc-3ffae20e4b34
sgnfitba

# ╔═╡ 648b7ab8-2b08-4f54-a312-f19596b13847
begin
	
	zfsgn = plot(ptsgn, txsgn, sgnfit.ypred, legend=true, label=lbl1)
	zfsgnba = plot(ptsgnba, txsgn, sgnfitba.ypred, legend=true, label=lbl2)
	plot(zfsgn, zfsgnba )
end

# ╔═╡ b56c6440-7fcc-4425-b0bd-5dde6795ae2c
hsgn.centers

# ╔═╡ b5de073e-3008-4ed8-b1d9-42b0ba14a6e5
df = tfit(hdata.centers , hdata.weights; pa0=[100.0, 1])

# ╔═╡ dfb29cee-3a98-4fd5-aead-2da2e77e777a
function hsums(hst; bw = 10, step=9)
	ebins = [hst.weights[1+n*bw:end] for n in 0:step]
	esum = [sum(e) for e in ebins]
	xenr =[hsgn.centers[1+n*bw] for n in 0:step]
	esum, xenr
end

# ╔═╡ d16741d3-f259-4cfe-a69f-a20f12d0ea88
fesum, xe = hsums(hsgn, bw = 5, step=19)

# ╔═╡ 074e55f8-d0a5-48bf-a86b-3b80aee57b6a
besum, _ = hsums(hsgnba, bw = 5, step=19)

# ╔═╡ f88cf9d2-f019-4e4e-b9ba-4522ca509e1a
bdc, _ = hsums(hdc, bw = 5, step=19)

# ╔═╡ 39c2e55d-e006-4882-9fd3-da4c8fd228c1
bfr, _ = hsums(hfree, bw = 5, step=19)

# ╔═╡ b719e5b3-9e1f-4137-83ee-19113b7d76e2
bdata, _ = hsums(hdata, bw = 5, step=19)

# ╔═╡ 93e2f0dd-81b8-40a9-9599-7a78e73e5b4a
cdata = bdata .- bfr .- bdc

# ╔═╡ 665affd3-addf-4bfc-9b96-d3312d229617
md"""
- number of background molecules selected = $num_bkg
- total photons measured $(sum(bdata))
- total photons corrected $(sum(cdata))
- total photons 1 molecule Ba2+ $(sum(besum))
"""

# ╔═╡ b742d915-38bf-4675-9e2f-a177db6a1a21
scatter(xe, bdata, yerr=sqrt.(abs.(bdata)), 
    xlabel ="t (ns)", ylabel= "counts", label="Data", markersize=2)

# ╔═╡ 8264cdba-b59a-435d-8b96-1a6279e2da92
begin
	#scatter(xe, fesum, yerr=sqrt.(abs.(fesum)), 
    #xlabel ="t (ns)", ylabel= "counts", label="free", markersize=2)
	plot(xe, fesum, xlabel ="t (ns)", ylabel= "counts", label="free")
	plot!(xe, besum, label="Ba2+")
	plot!(xe, bdc, label="DC")
	
	scatter!(xe, cdata, yerr=sqrt.(abs.(bdata)), 
    xlabel ="t (ns)", ylabel= "counts", label="Data Corrected", markersize=2)
	#plot!(xe, bdata,  lw=2, label="Data")
end

# ╔═╡ Cell order:
# ╟─c72e53d2-a21b-11ed-29d3-c72f8dda6658
# ╟─7926a440-16a6-4e7a-b9a4-ebc1102699e1
# ╟─f9be2504-3c14-4ab2-bd67-90f3859cd96c
# ╟─842e709a-8c93-4764-844b-e9193d6a2ed5
# ╟─13513708-7ed8-4a6b-80c8-dd0456aea681
# ╠═6f0c7d83-04e6-456c-8a65-b6ecccf4ae0b
# ╟─b0d80e2d-caaa-449d-a509-3487536dc9ae
# ╠═78ff48a5-a903-44a3-af8d-ffc4c832cea5
# ╠═6d047015-90c5-4774-9fa6-3e5f49160c19
# ╠═aadf254a-47ea-417f-9f5e-734486445e32
# ╠═c6baa96a-2793-4f82-a932-be627f989964
# ╠═2b3561b7-36b0-485b-8532-de1d89702423
# ╟─abf3d6d5-ca26-434c-b248-4c2395c79451
# ╟─38ded754-c916-44a4-b673-02fcb9055cf8
# ╠═a276091d-0c83-4a75-be20-51ed938d7415
# ╟─babadc24-94dd-4c00-80ae-42a0bbc5fcc7
# ╟─e02432fe-7b16-42ae-af67-9959ad68e1a9
# ╟─ec58050a-51bd-49fc-83f3-77980db52d7c
# ╟─7dad7b06-d969-43bf-8a15-30f30a05a72b
# ╟─a1e65010-0e0f-4e6e-917a-b643a92b9ec4
# ╟─7048b1a0-bd71-4f29-be1d-b18cd7191ec1
# ╠═42a0fe7f-2945-4804-ae7e-bb5763497853
# ╟─f69269ad-592c-4ea8-b228-1b97177299b4
# ╟─15ebb995-0b3b-45c0-bed9-b6ecf6d101f8
# ╟─960472cc-bdfa-40c4-b4d6-31bf8e840598
# ╟─ad2254eb-efe5-44b5-a1d1-c023a4456801
# ╠═24c43b9a-cff6-462f-9d13-95e3bdb71792
# ╠═31c8c858-2178-4fce-a772-d116dca3f496
# ╟─8d891f09-9745-4b83-8435-028f64f7ecb0
# ╟─ae7c9b4b-38f0-4bd6-ba62-4ebcbf4b04fc
# ╟─be9dfc6d-c627-4a38-8b46-4d9b1310162c
# ╟─0a2b607f-4e30-41be-a73d-c0daf5c8b855
# ╟─e672f7bc-5020-40ae-9706-a12640607b6f
# ╟─aa9e5f13-28e3-424b-a6b6-0cfb0c1c7625
# ╟─05e4e447-4012-46e2-a6cf-e14c49e31d40
# ╠═1eb59032-324e-4f36-8c6e-647f84e5c570
# ╟─45c79d1c-650e-4db7-8e9c-4d7cf4851fa9
# ╟─be76b61f-843e-4dbc-b13a-4aebb4cfba5f
# ╟─cb9b6dfa-d54c-4635-9bff-1a89b73ad93a
# ╟─7e9efea1-061b-4d57-a6fb-b66a955d2926
# ╟─8cd2b3cb-660b-48a5-a0e7-402fb835b499
# ╟─a57f3c34-3319-4275-bdda-982ff5f67445
# ╠═7aef4cca-a6ce-4f16-ba7d-e7914d054247
# ╟─61ffc4ea-5107-453a-8048-9bdc11a94659
# ╟─d3810185-890e-4c52-a947-e53d6d9e0ebc
# ╠═2e403946-4a25-45ee-935b-2f6af7f87658
# ╠═1a36aa14-a7ea-4db0-8ce6-4f52e0c1bd04
# ╠═4ce0fdf3-dc93-456d-8229-c05b0ef54268
# ╠═3d3fa4a0-be5c-450f-9b0a-8500878f38b1
# ╠═0588ae25-2f71-48d1-9fbf-e682f6161f24
# ╠═1ae83c20-d8ec-4717-9275-8828bc1e193a
# ╠═dd347540-1ee9-40e9-8f69-8ecda506cbbf
# ╠═2d68c170-0d91-4d84-975a-0165b86a2883
# ╟─eedef1fc-4b5a-49b8-8f61-ae99cc27215e
# ╠═70a95466-f7b4-4e5d-a0bf-b3c7a3150dc6
# ╠═299ee6bd-b011-4937-a46c-c8febfcf61c7
# ╠═ea7e2262-ca35-4537-b914-5e1971cbfc20
# ╠═0822b115-575e-4915-87e9-1264abf6813a
# ╟─cef81cd5-447d-459e-be70-4a7c531ed504
# ╟─12d7b082-0a19-495c-b5db-feb4ab700267
# ╠═a23a9fbe-772a-4cc0-81de-f3b5dc58d766
# ╠═0e6e6bfb-59c1-484a-a2a1-4243cfcce32b
# ╠═d6d35e53-a712-4e89-9e70-8c166f6b7550
# ╠═4012d326-3a8e-426e-899f-5dddf7674b35
# ╠═648b7ab8-2b08-4f54-a312-f19596b13847
# ╠═903ac25f-493f-4b70-87b1-480b8d8396fb
# ╠═d2ed2a93-e059-4900-b0dc-3ffae20e4b34
# ╠═9be5bbaa-32c4-41fb-b917-ed747dbeb213
# ╠═bf5af920-6a40-42af-8697-8fb089976a9f
# ╠═b56c6440-7fcc-4425-b0bd-5dde6795ae2c
# ╠═d16741d3-f259-4cfe-a69f-a20f12d0ea88
# ╠═074e55f8-d0a5-48bf-a86b-3b80aee57b6a
# ╠═f88cf9d2-f019-4e4e-b9ba-4522ca509e1a
# ╠═39c2e55d-e006-4882-9fd3-da4c8fd228c1
# ╠═b719e5b3-9e1f-4137-83ee-19113b7d76e2
# ╠═93e2f0dd-81b8-40a9-9599-7a78e73e5b4a
# ╠═1f19dfd5-cf6a-48a8-a65d-0fbc752bae99
# ╠═665affd3-addf-4bfc-9b96-d3312d229617
# ╠═b742d915-38bf-4675-9e2f-a177db6a1a21
# ╠═2e4620c0-9abb-4359-987b-839ebeccef96
# ╠═8264cdba-b59a-435d-8b96-1a6279e2da92
# ╠═de748626-2bad-4637-b948-bd2a2ad72e97
# ╠═0fd69074-e66b-40dc-8810-5821bd1c6ad8
# ╠═b5de073e-3008-4ed8-b1d9-42b0ba14a6e5
# ╠═1279174a-3e68-4dad-9b41-7c798807bd44
# ╠═dfb29cee-3a98-4fd5-aead-2da2e77e777a
# ╠═a8bf7b5f-286e-4bb7-9dc5-91d01878ed05
# ╠═b2bac877-7098-43d0-8241-c01b6e219bd1
# ╠═7f6b34e2-ade3-4a0a-928b-66bee72456ab
# ╠═c51c7590-35aa-4ba2-8539-61b42dfab5b2
# ╠═8e185d65-607c-435c-9f1b-da738060ddd1
# ╠═b1d523c2-83cc-4893-b1a7-c62d3fd40a9d
# ╠═b5345760-ae72-4125-a502-39318a87f5a3
# ╠═5c3ac9c2-8b16-4072-ba9a-1d336cea52bf
# ╠═277be488-9840-4f62-874f-d4f9fbb893d8
# ╠═cbe94ca2-0b1f-4246-8843-6b8e4e5e55cc
