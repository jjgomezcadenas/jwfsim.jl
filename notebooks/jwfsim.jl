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
import PhysicalConstants.CODATA2018: N_A

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

- In order to simulate the images we need to know the number of pixels, their sizes, binning, noise, QE...
"""


# ╔═╡ d7f3fe4f-455a-435b-98f7-60c00e96c1eb
begin
	ntpxx=2048
	ntpxy=2048
	stpxx=2
	stpxy=2
	binning=16
	ssensx=stpxx*ntpxx/1000
	ssenst=stpxy*ntpxy/1000
	npxx=ntpxx/sqrt(binning)
	npxy=ntpxy/sqrt(binning)
	spxx=stpxx*sqrt(binning)
	spxy=stpxy*sqrt(binning)
	md"""
- Pixel array = $ntpxx x $ntpxy
- Pixel size = $stpxx x $stpxy (microns)
- Sensor size= $ssensx x $ssensx (mm)
- Binning= $binning
- Pixel array (with binning) = $npxx x $npxy
- Pixel size (with binning)= $spxx x $spxy (microns)

"""
	
end

# ╔═╡ 03357bcc-b54d-4104-b2c5-cbc5ce31b721
md""" 
## Other specs

- Magnification, power losses throug the set up...
"""

# ╔═╡ cc82232e-fc86-413f-8b73-8044b4de6260
begin
	mag=100
	md"""
	-Magnification= x $mag
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

- Once we have simulated the positions of the molecules in our sample we can 

"""

# ╔═╡ c2289697-ab87-4f11-81e6-bee3439ca3cb
md"""
# Functions
"""

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

# ╔═╡ Cell order:
# ╠═625c1998-828c-11ed-2c7b-f731be119aaf
# ╠═4b3981d2-d1f9-4f56-a44d-bcc55d0b1d87
# ╠═d75ba333-2445-4828-a66c-3c33d96de16f
# ╠═3007fa59-d2f0-4248-820d-744f04dc984d
# ╠═47414778-e472-4d4f-b05b-1213cc328a56
# ╠═558dbadc-cb83-4cfb-a05f-1f57b62bcf1d
# ╠═6497410d-d1eb-4c06-9f98-1da578fb683f
# ╠═40bcafee-88a9-4c7b-a611-d0599e4567e9
# ╟─27b5dc4d-b688-4414-8fcb-97c4ef1d680d
# ╟─d8ea1458-304a-462d-abdc-24ba32d3d8b0
# ╟─92db42ec-1f18-4277-8612-7613c7c32c40
# ╟─d7f3fe4f-455a-435b-98f7-60c00e96c1eb
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
# ╟─b32cd1dd-a93a-41e6-a9cd-92bb9295f3ad
# ╠═766d17b5-499a-4a28-8628-3031b023a851
# ╠═c2289697-ab87-4f11-81e6-bee3439ca3cb
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
