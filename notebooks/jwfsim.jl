### A Pluto.jl notebook ###
# v0.19.19

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

# ╔═╡ 790af81b-45ce-4fae-ab0b-7fb52b87a2af
md"""
# Simulation of the spatial distribution of the molecules:

Ingredientes:

1. Molecular density (number of molecules/unit volume). Tipically the input is molecular concentration in M (mol/liter), e.g, micromolar = 10^-6 mol/l 
2. Target surface, e.g, the surface where the molecules are deposited (tipically 1cm2)
"""

# ╔═╡ b65232d2-f8ed-4afc-969a-417d78405477
md"""
## Molecular concentration and target length

- It is conveniente to use -log(ρ), where ρ is the molecular concentration (in M), thus, for example, to specifiy 1 μM (10^-6 M), -log(ρ) = 6

- Also convenient to assume that the target is a square slice, and specify just its longitudinal dimensions (e.g, in cm)
"""

# ╔═╡ abe2369a-8a9a-4e0e-9707-4184e6d6f339
md""" Select -log concentration (-log(ρ)) in M: $(@bind logr NumberField(1.0:15.0, default=5.0))"""


# ╔═╡ 55352380-ff6d-42f6-863a-e5de8cced158
md""" Select length (in cm) of (square) glass: $(@bind scm NumberField(1:10.0, default=1.00))"""

# ╔═╡ 766d17b5-499a-4a28-8628-3031b023a851
md"""
## Simulate image

- once we have the mean value of the number of molecules per micron square (which we assume to be the granularity of our simulation, and call nmu2), we need to define the field of view (FOV) in pixels (we can take 1 pixel = 1 micron)

"""

# ╔═╡ c1b960ad-b7f5-46ff-985f-0f304bdb0a48
md""" Select FOV size (in pixels = μm) : $(@bind fovl NumberField(10:1000, default=30))"""

# ╔═╡ 4153e45e-3973-4b37-adfc-49feda60e3a2
md"""
Next step is to generate in each pixel of the FOV (with dimensions fovl^2), a population of molecules. This is done throwing Poisson numbers (with mu=nmu2) in each pixel. 
"""

# ╔═╡ c2289697-ab87-4f11-81e6-bee3439ca3cb
md"""
# Functions
"""

# ╔═╡ 736f187e-e5c0-4458-8438-5860ad0f2f98
function toncm3(r::Float64)
	uconvert(cm^-3,10^-r*M*N_A)/cm^-3
end

# ╔═╡ de442740-e0c0-4b7c-93d4-07974858afb1
md"""
- The selected concentration correspondes to $(round((toncm3(logr)), sigdigits=3)) N/cm3 (N = molecules)
"""

# ╔═╡ cc6a3d4b-b032-4092-b29a-0f7b6854852b
function toncm2(r::Float64, s::Float64)
	(uconvert(cm^-3,10^-r*M*N_A) /(s*cm^-1))/cm^-2 
end

# ╔═╡ 03bf3bfb-2ec3-4ffa-9e83-37e21634a59d
function tonpers(r::Float64, s::Float64, unit)
	toncm2(r::Float64, s::Float64) / (uconvert(cm^-2, unit)/cm^-2)
end

# ╔═╡ 1dfc12c7-ad06-48f7-a4f4-7ed945db50b0
md"""
- The selected number of molecules is to $(round(toncm2(logr, scm), sigdigits=3)) N/cm2 (N = molecules)
- This correspondes to $(round(tonpers(logr, scm, 1.0μm^-2))) N/μ^2
"""

# ╔═╡ 538477a1-748f-443c-a36b-cba2f1a7645b
nmu2 = tonpers(logr, scm, 1.0μm^-2)

# ╔═╡ 588cd3a3-9b34-4fbb-9959-0f9219ef43b7
pois_rand(nmu2)

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

# ╔═╡ 4c6f43a2-5354-4dc5-99bf-0f136c648c33
nfov = populate_fov(fovl, nmu2)

# ╔═╡ 629c7167-4314-4668-9433-c4ab54dc85b8
length(filter(x -> x == 1.0, nfov)) 

# ╔═╡ fba1abc3-09c8-4f33-88a9-25e8a9748b69
nfov[nfov .== 1.0]

# ╔═╡ a9873b66-6336-406f-a9b1-0fe8a1fc0dc3
heatmap(nfov)

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

# ╔═╡ Cell order:
# ╠═625c1998-828c-11ed-2c7b-f731be119aaf
# ╠═4b3981d2-d1f9-4f56-a44d-bcc55d0b1d87
# ╠═d75ba333-2445-4828-a66c-3c33d96de16f
# ╠═3007fa59-d2f0-4248-820d-744f04dc984d
# ╠═47414778-e472-4d4f-b05b-1213cc328a56
# ╠═558dbadc-cb83-4cfb-a05f-1f57b62bcf1d
# ╠═6497410d-d1eb-4c06-9f98-1da578fb683f
# ╠═40bcafee-88a9-4c7b-a611-d0599e4567e9
# ╠═27b5dc4d-b688-4414-8fcb-97c4ef1d680d
# ╠═790af81b-45ce-4fae-ab0b-7fb52b87a2af
# ╠═b65232d2-f8ed-4afc-969a-417d78405477
# ╠═abe2369a-8a9a-4e0e-9707-4184e6d6f339
# ╠═de442740-e0c0-4b7c-93d4-07974858afb1
# ╠═55352380-ff6d-42f6-863a-e5de8cced158
# ╠═1dfc12c7-ad06-48f7-a4f4-7ed945db50b0
# ╠═766d17b5-499a-4a28-8628-3031b023a851
# ╠═c1b960ad-b7f5-46ff-985f-0f304bdb0a48
# ╠═4153e45e-3973-4b37-adfc-49feda60e3a2
# ╠═538477a1-748f-443c-a36b-cba2f1a7645b
# ╠═588cd3a3-9b34-4fbb-9959-0f9219ef43b7
# ╠═4c6f43a2-5354-4dc5-99bf-0f136c648c33
# ╠═629c7167-4314-4668-9433-c4ab54dc85b8
# ╠═fba1abc3-09c8-4f33-88a9-25e8a9748b69
# ╠═a9873b66-6336-406f-a9b1-0fe8a1fc0dc3
# ╠═c2289697-ab87-4f11-81e6-bee3439ca3cb
# ╠═736f187e-e5c0-4458-8438-5860ad0f2f98
# ╠═cc6a3d4b-b032-4092-b29a-0f7b6854852b
# ╠═03bf3bfb-2ec3-4ffa-9e83-37e21634a59d
# ╠═a988f5a4-ce60-4820-abea-3e3c79854d89
# ╠═f893442d-230c-4ecd-b618-2ad89393963b
# ╠═c0f7e1ad-37d2-4b42-8683-5d2088a8c4ea
# ╠═5850de84-9c5c-4036-b941-97e604d9b4a0
