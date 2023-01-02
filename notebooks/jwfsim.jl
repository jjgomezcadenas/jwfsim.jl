### A Pluto.jl notebook ###
# v0.19.14

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

# ╔═╡ 625c1998-828c-11ed-2c7b-f731be119aaf
using Pkg; Pkg.activate(ENV["JWfSim"])

# ╔═╡ 4b3981d2-d1f9-4f56-a44d-bcc55d0b1d87
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Images
	using Colors
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
end

# ╔═╡ d75ba333-2445-4828-a66c-3c33d96de16f
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 3007fa59-d2f0-4248-820d-744f04dc984d
import PhysicalConstants.CODATA2018: N_A,h,SpeedOfLightInVacuum

# ╔═╡ be18cfbb-3c4a-4e29-8903-0261ed57051e
c_0=SpeedOfLightInVacuum

# ╔═╡ 47414778-e472-4d4f-b05b-1213cc328a56
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

# ╔═╡ 558dbadc-cb83-4cfb-a05f-1f57b62bcf1d
jwf = ingredients("../src/jwfsim.jl")

# ╔═╡ 6497410d-d1eb-4c06-9f98-1da578fb683f
PlutoUI.TableOfContents(title="Wide Field Simulation", indent=true)

# ╔═╡ 40bcafee-88a9-4c7b-a611-d0599e4567e9
md"""
# Contributors

- J.J. Gómez Cadenas (jjgomezcadenas@gmail.com)
- M. Elorza
- P. Herrero
"""

# ╔═╡ 27b5dc4d-b688-4414-8fcb-97c4ef1d680d
md"""
# Introduction 
Program Structure:

The simulation has several distinct ingredients:

1. Simulation of the spatial distribution of the molecules.
2. Simulation of the molecule response (e.g, cross section, quantum yield)
3. Simulation of the laser setup (how many photons reach the sample per unit area, given a laser power and the wide field setup)

"""

# ╔═╡ d8ea1458-304a-462d-abdc-24ba32d3d8b0
md"""
# Set up
"""

# ╔═╡ 92db42ec-1f18-4277-8612-7613c7c32c40
md""" 
## Camera specs

 In order to simulate the images we need to know the number of pixels, their sizes, binning, noise, QE...
"""


# ╔═╡ d7f3fe4f-455a-435b-98f7-60c00e96c1eb
begin
	ntpxx=2048
	ntpxy=2048
	stpxx=2
	stpxy=2
	binning=16
	ssensx=stpxx*ntpxx
	ssensy=stpxy*ntpxy
	npxx=Int(ntpxx/sqrt(binning))
	npxy=Int(ntpxy/sqrt(binning))
	spxx=stpxx*sqrt(binning)
	spxy=stpxy*sqrt(binning)
	md"""
- Pixel array = $ntpxx x $ntpxy
- Pixel size = $stpxx x $stpxy (microns)
- Sensor size= $ssensx x $ssensy (microns)
- Binning= $binning
- Pixel array (with binning) = $npxx x $npxy
- Pixel size (with binning)= $spxx x $spxy (microns)

"""
	
end

# ╔═╡ f506db9e-3fe3-4963-856a-be91d91e8266
md""" 
## Laser specs
"""

# ╔═╡ 9eea388f-c453-480d-b022-2340fe880b89
begin
	pwr=1mW
	emwl=375nm
	sbx=2
	sby=2
	md"""
	- Power = $pwr mW
	- Emission wavelength = 375 nm
	- Beam width (1 sigma assuming gaussian beam) =  $sbx x $sby mm 
	"""
	
end

# ╔═╡ 63e4aa43-b606-4077-9205-e49fd7ac3b21
begin
center = [0, 0]
sigmas = [sbx,sby]
p = MvNormal(center,sigmas)

X = range(-3*sbx, 3*sbx, length=100)
Y = range(-3*sby, 3*sby, length=100)
Z = [pdf(p, [x,y]) for y in Y, x in X] # Note x-y "for" ordering
p0=plot(X,Y,Z,st=:surface,xlabel="x (mm)",ylabel="y (mm)")

p1=contourf(X, Y, Z, xlabel="x (mm)",ylabel="y (mm)")
plot(p0,p1,size=(900,300))
end

# ╔═╡ 03357bcc-b54d-4104-b2c5-cbc5ce31b721
md""" 
## Other specs

Magnification, power losses throug the set up...
"""

# ╔═╡ cc82232e-fc86-413f-8b73-8044b4de6260
begin
	mag=100
	md"""
	- Magnification= x $mag
	"""
end

# ╔═╡ 790af81b-45ce-4fae-ab0b-7fb52b87a2af
md"""
# Simulation of the spatial distribution of the molecules:

Ingredientes:

1. Molecules per micron square (nmu2). The number of molecules in each area of $\mu$m x $\mu$m.
2. Target surface, e.g, the surface where the molecules are deposited (tipically 1cm2)
"""

# ╔═╡ 018997ec-0c72-4677-b6d1-e8dafa1c2811
schema = ["Molecules per micron square", "Molarity"]

# ╔═╡ 0782a2c9-7d5a-49bc-8ddb-71c57238de8b
md""" Select model input : $(@bind sch Select(schema))"""

# ╔═╡ 55352380-ff6d-42f6-863a-e5de8cced158
md""" Select length (in mm) of (square) glass: $(@bind smm NumberField(0.0:0.1:10.0, default=1))"""

# ╔═╡ abe2369a-8a9a-4e0e-9707-4184e6d6f339
if sch=="Molecules per micron square"
	md""" Select the molecules per unit area (in molecules per $\mu$m$^2$): $(@bind a NumberField(0.0:100.0, default=1.0))"""
end

# ╔═╡ 00cfc948-31e4-4d17-bf44-0a25c911cf68
if sch=="Molarity"
	md""" Select the molarity of the solution (in -log(M)): $(@bind b NumberField(5.0:15.0, default=12))"""
end

# ╔═╡ 39ec5a5f-8e06-4143-8016-39fc9020aaa1
md"""
Once we have the mean value of the number of molecules per micron square (which we assume to be the granularity of our simulation, and call nmu2), we can simulate the molecule distribution on the sample.

"""

# ╔═╡ 766d17b5-499a-4a28-8628-3031b023a851
md"""
# Simulate image
"""

# ╔═╡ a24aab57-2145-40de-bfd0-b768fef6d01d
md"""
## FOV
The field of view of the camera scales with the magnification of the system. The size  of the field of view will be the sensor size divided by the magnfication, and the resolution will be given by the pixel size divided by the magnification.
"""

# ╔═╡ f7f4ba8b-e1a6-4632-a785-767b2bfb646f
begin
	sfovx=ssensx/mag
	sfovy=ssensy/mag
	spxfovx=spxx/mag
	spxfovy=spxy/mag
	md"""
	- Size of the FOV = $sfovx x $sfovy (microns)
	- Resolution = $spxfovx x $spxfovy (microns)
	"""
end

# ╔═╡ 245c0ff7-22b2-47aa-969d-5cf4eca9831e
md""" ### Molecule distribution in th FOV"""

# ╔═╡ a45829e4-ac16-4e5d-9ce2-6007d62f37e8
md""" ### Illumination in th FOV"""

# ╔═╡ 907729d2-6b55-4afb-a3fd-5a450dd57a63
begin
center1 = [npxx/2, npxy/2]
sigmas1 = [sbx*1000/spxfovx/mag,sby*1000/spxfovy/mag]
pp = MvNormal(center1,sigmas1)

X1 = range(0, npxx, length=npxx)
Y1 = range(0, npxy, length=npxy)
Z1 = [pdf(pp, [x,y]) for y in Y1, x in X1] # Note x-y "for" ordering
p01=plot(X1,Y1,Z1,st=:surface,xlabel="x (mm)",ylabel="y (mm)")

p11=contourf(X1, Y1, Z1, xlabel="x (mm)",ylabel="y (mm)")
plot(p01,p11,size=(900,300))
end

# ╔═╡ adebe61f-6355-46f3-b6d8-86fd50484366
begin
sigmafovx=sigmas1[1]
sigmafovy=sigmas1[2]

md""" Sigma of the spot in the FOV is $sigmafovx microns in x axis and $sigmafovy microns in y axis"""
end

# ╔═╡ 2e1ff7fa-1f63-4d52-87e2-b0d61fc1f85d
heatmap(Z1)

# ╔═╡ c1af380d-8f25-4e3b-b5b3-5ac78f0a9913
md""" ### Molecule response
The number of emitted photons per second from the molecules will be $N_{emmitted}(x,y)=N_{incident}(x,y)N_{molecules}\sigma$ where $N_{incident}$ is the number of photons hitting the sample at (x,y) position per unit of time, $N_{molecules}$ the number of molecules por unit area, and $\sigma$ the interaccion cross section. If the number of incident photons is too high there could be a saturation on the emitted photons.
"""

# ╔═╡ 090f9903-8eb8-4846-909e-2cd2ffeb536e
md""" - $N_{incident}(x,y)=\frac{P(x,y)}{e(\lambda)}=\frac{P_{total}G(x,y)}{e(\lambda)}$ where $P_{total}$ is the laser power, $e(\lambda)$ the energy per photon, and $G(x,y)$ the power distribution (in this case gaussian).
"""

# ╔═╡ f778f467-001d-4484-a20f-fd85e10558b1
md""" - $N_{molecules}$ is an input of the notebook. 
"""

# ╔═╡ 58163cae-5a87-4de9-8747-f47dc0398389
md""" - $\sigma$ is unknown (as far as I know), but it should be a function of the photon energy.  
"""

# ╔═╡ 8c9160ca-0048-4de6-8a39-ad00fe44f680
begin
md""" In this case:"""
	
end

# ╔═╡ f07e7613-1402-4abc-8b06-dc98df2b8c98
begin
ei=uconvert(eV,h*c_0/emwl)
md""" $e(\lambda)$=$ei """
end

# ╔═╡ bd7a2018-d882-4891-8f15-d90e92cb1077
begin
md""" $P_{total}$=$pwr"""
end

# ╔═╡ c2289697-ab87-4f11-81e6-bee3439ca3cb
md"""
# Functions
"""

# ╔═╡ 8d920a11-2050-4562-b55c-0b916cabf42c
function toimage(data,xran,yran,spxfovx,spxfovy)
	x=data[1,:]
	y=data[2,:]
	fit(Histogram,(x, y), (xran[1]:spxfovx:xran[2], yran[1]:spxfovy:yran[2])).weights
end

# ╔═╡ 13b7f5db-8593-4511-b540-6f40e8eb3498
function poss_fov(poss,xran,yran,spxfovx,spxfovy)
	possfov=[0,0]
	for pos in eachrow(poss)
		if (xran[1]<pos[1])&&(pos[1]<xran[2])
			if (yran[1]<pos[2])&&(pos[1]<yran[2])
				possfov=hcat(possfov,pos)
			end
		end
	end
	possfov=possfov[:,2:end]
end

# ╔═╡ 736f187e-e5c0-4458-8438-5860ad0f2f98
function toncm3(r::Float64)
	uconvert(cm^-3,10^-r*M*N_A)/cm^-3
end

# ╔═╡ 5ab2275a-a576-4c2a-ae01-23748092d51a
uconvert(cm^-3,10^-5*M*N_A)

# ╔═╡ cc6a3d4b-b032-4092-b29a-0f7b6854852b
function toncm2(r::Float64, s::Float64)
	(uconvert(cm^-3,10^-r*M*N_A) /(s*cm^-1))/cm^-2 
end

# ╔═╡ 18fc7474-c7b8-4092-9dfa-92f07ce3cb2a
begin
if sch=="Molecules per micron square"
	nmu2=a
end
if sch=="Molarity"
	nmu2=toncm2(b,smm)/1e8
end
md""" Molecules per micron square (nmu2)=$nmu2 """
end


# ╔═╡ 84e73dfd-a9a3-4bd4-93ce-2fc4a4e78f1b
begin
smu=smm*1000
ntot=Int(trunc(nmu2*smu^2))
md""" Total number of generated molecules = $ntot """
end

# ╔═╡ b32cd1dd-a93a-41e6-a9cd-92bb9295f3ad
begin
poss=smu*rand(Float64, (ntot,2))
md"""
Positions (x,y) in microns of the generated molecules are stored in 'poss' array
"""
end

# ╔═╡ d9e3e919-f1b6-4b0b-9e07-04d5855f3587
md""" Select x0 position on the sample (in microns). It ranges from 0 to $smu : $(@bind x0 NumberField(0:smu, default=smu/2))"""

# ╔═╡ 50c17b83-57f5-48dc-aaa1-22e16f627022
md""" Select y0 position on the sample (in microns). It ranges from 0 to $smu : $(@bind y0 NumberField(0:smu, default=smu/2))"""

# ╔═╡ 245205d3-ef94-46cf-ba0a-48810c201610
begin
xran=[x0-sfovx/2,x0+sfovx/2]
yran=[y0-sfovy/2,y0+sfovy/2]
end

# ╔═╡ da45d069-b08d-4221-9e7d-32371090e94d
begin
possfov=poss_fov(poss, xran, yran, spxfovx,spxfovy)
poss_array=toimage(possfov,xran,yran,spxfovx,spxfovy)
heatmap(poss_array)

end

# ╔═╡ 3e1f4f31-a563-4063-8ef9-653a508b9d86
begin
	nmu2u=nmu2/μm/μm
	md""" $N_{molecules}$=$nmu2u """
end

# ╔═╡ 03bf3bfb-2ec3-4ffa-9e83-37e21634a59d
function tonpers(r::Float64, s::Float64, unit)
	toncm2(r::Float64, s::Float64) / (uconvert(cm^-2, unit)/cm^-2)
end

# ╔═╡ a988f5a4-ce60-4820-abea-3e3c79854d89
function populate_fov(fovl::Integer, nmu2::Float64)
	FOV = zeros(fovl, fovl)
	for i in 1:fovl
		for j in 1:fovl
			FOV[i,j] = pois_rand(nmu2)
		end
	end
	FOV
end

# ╔═╡ c45c60cd-b581-4be1-bfbf-fca725221e4b


# ╔═╡ f893442d-230c-4ecd-b618-2ad89393963b
md"""
# Tests
"""

# ╔═╡ c0f7e1ad-37d2-4b42-8683-5d2088a8c4ea
"""
Given cc in M and ll in cm, compute number of molecules per square micron
"""
function nmcm3(cc::Float64=1e-12, ll::Float64=1)
	ltcm3 = 10^3  # L to cm3
	cmtmu = 10^4
 	ncm3 = cc / (ltcm3) * N_A * mol  # n/cm3
	ncm2 = ncm3 /ll # n/cm2
	nmu2 = ncm2 / cmtmu^2
end

# ╔═╡ 5850de84-9c5c-4036-b941-97e604d9b4a0
nmcm3(1e-12, 1.0) ≈ tonpers(12.0, 1.0, 1.0μm^-2)

# ╔═╡ 9541c8a6-318a-4266-b388-4b69ce7eb5b1
pois_rand

# ╔═╡ 0ab9fa82-6a93-48c9-9419-eb2e15f1df99
5nm/2W

# ╔═╡ Cell order:
# ╠═625c1998-828c-11ed-2c7b-f731be119aaf
# ╠═4b3981d2-d1f9-4f56-a44d-bcc55d0b1d87
# ╠═d75ba333-2445-4828-a66c-3c33d96de16f
# ╠═3007fa59-d2f0-4248-820d-744f04dc984d
# ╠═be18cfbb-3c4a-4e29-8903-0261ed57051e
# ╠═47414778-e472-4d4f-b05b-1213cc328a56
# ╠═558dbadc-cb83-4cfb-a05f-1f57b62bcf1d
# ╠═6497410d-d1eb-4c06-9f98-1da578fb683f
# ╟─40bcafee-88a9-4c7b-a611-d0599e4567e9
# ╟─27b5dc4d-b688-4414-8fcb-97c4ef1d680d
# ╟─d8ea1458-304a-462d-abdc-24ba32d3d8b0
# ╟─92db42ec-1f18-4277-8612-7613c7c32c40
# ╠═d7f3fe4f-455a-435b-98f7-60c00e96c1eb
# ╟─f506db9e-3fe3-4963-856a-be91d91e8266
# ╠═9eea388f-c453-480d-b022-2340fe880b89
# ╠═63e4aa43-b606-4077-9205-e49fd7ac3b21
# ╟─03357bcc-b54d-4104-b2c5-cbc5ce31b721
# ╟─cc82232e-fc86-413f-8b73-8044b4de6260
# ╟─790af81b-45ce-4fae-ab0b-7fb52b87a2af
# ╟─018997ec-0c72-4677-b6d1-e8dafa1c2811
# ╟─0782a2c9-7d5a-49bc-8ddb-71c57238de8b
# ╟─55352380-ff6d-42f6-863a-e5de8cced158
# ╟─abe2369a-8a9a-4e0e-9707-4184e6d6f339
# ╟─00cfc948-31e4-4d17-bf44-0a25c911cf68
# ╟─18fc7474-c7b8-4092-9dfa-92f07ce3cb2a
# ╟─39ec5a5f-8e06-4143-8016-39fc9020aaa1
# ╟─84e73dfd-a9a3-4bd4-93ce-2fc4a4e78f1b
# ╠═b32cd1dd-a93a-41e6-a9cd-92bb9295f3ad
# ╟─766d17b5-499a-4a28-8628-3031b023a851
# ╟─a24aab57-2145-40de-bfd0-b768fef6d01d
# ╟─f7f4ba8b-e1a6-4632-a785-767b2bfb646f
# ╠═d9e3e919-f1b6-4b0b-9e07-04d5855f3587
# ╠═50c17b83-57f5-48dc-aaa1-22e16f627022
# ╟─245205d3-ef94-46cf-ba0a-48810c201610
# ╟─245c0ff7-22b2-47aa-969d-5cf4eca9831e
# ╟─da45d069-b08d-4221-9e7d-32371090e94d
# ╟─a45829e4-ac16-4e5d-9ce2-6007d62f37e8
# ╟─907729d2-6b55-4afb-a3fd-5a450dd57a63
# ╟─adebe61f-6355-46f3-b6d8-86fd50484366
# ╟─2e1ff7fa-1f63-4d52-87e2-b0d61fc1f85d
# ╟─c1af380d-8f25-4e3b-b5b3-5ac78f0a9913
# ╟─090f9903-8eb8-4846-909e-2cd2ffeb536e
# ╟─f778f467-001d-4484-a20f-fd85e10558b1
# ╟─58163cae-5a87-4de9-8747-f47dc0398389
# ╟─8c9160ca-0048-4de6-8a39-ad00fe44f680
# ╟─f07e7613-1402-4abc-8b06-dc98df2b8c98
# ╟─bd7a2018-d882-4891-8f15-d90e92cb1077
# ╟─3e1f4f31-a563-4063-8ef9-653a508b9d86
# ╠═c2289697-ab87-4f11-81e6-bee3439ca3cb
# ╠═8d920a11-2050-4562-b55c-0b916cabf42c
# ╠═13b7f5db-8593-4511-b540-6f40e8eb3498
# ╠═736f187e-e5c0-4458-8438-5860ad0f2f98
# ╠═5ab2275a-a576-4c2a-ae01-23748092d51a
# ╠═cc6a3d4b-b032-4092-b29a-0f7b6854852b
# ╠═03bf3bfb-2ec3-4ffa-9e83-37e21634a59d
# ╠═a988f5a4-ce60-4820-abea-3e3c79854d89
# ╠═c45c60cd-b581-4be1-bfbf-fca725221e4b
# ╠═f893442d-230c-4ecd-b618-2ad89393963b
# ╠═c0f7e1ad-37d2-4b42-8683-5d2088a8c4ea
# ╠═5850de84-9c5c-4036-b941-97e604d9b4a0
# ╠═9541c8a6-318a-4266-b388-4b69ce7eb5b1
# ╠═0ab9fa82-6a93-48c9-9419-eb2e15f1df99
