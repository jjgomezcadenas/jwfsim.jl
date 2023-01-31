### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 6a3aa124-9c93-11ed-0dbf-d142a80cde9f
using Pkg; Pkg.activate(ENV["JWfSim"])


# ╔═╡ 4e9f0e0d-ee89-40d2-b6df-6e7f0d0e7d6d
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

# ╔═╡ d57a99a6-008d-4d00-8a00-ee1315b55445
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M


# ╔═╡ a35c53c3-a162-492d-a9ff-c0aa96429b55
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 4f68d411-cb03-42ea-b9be-e37fdf273060
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

# ╔═╡ 8b2dd016-f4ba-4872-895f-cbd54d53252d
jwf = ingredients("../src/jwfsim.jl")

# ╔═╡ 605a3cb2-feed-42ae-9816-62f6f6bc7d62
PlutoUI.TableOfContents(title="Ru & Ir", indent=true)

# ╔═╡ 539b8d7b-e2f0-4a58-a868-03249ebb7ccd
md"""
# Ru cross sections
"""

# ╔═╡ af94bd47-0817-48bf-b8e8-b53843ba7109
md"""
## Load Ru cross section data frames
"""

# ╔═╡ f8e49a64-379f-4b23-9e87-4fe3a32e6c0e
pdata = joinpath(ENV["JWfSim"], "data")

# ╔═╡ 41e2bf4a-7543-4481-9f78-bff2629b7b55
mecrusl ="VS014_RuSL_1E-5_molar_extinction coefficient.csv"

# ╔═╡ 22edeacc-9aff-45df-9c61-750f8386b84f
#run(`ls /Users/jjgomezcadenas/Projects/jwfsim/data/`)

# ╔═╡ bb8b43d3-41c2-4c6f-999d-e1146d2b351b
begin
	mecrusldf = sort(jwf.jwfsim.load_df_from_csv(pdata, mecrusl, jwf.jwfsim.spG))
	names(mecrusldf)
	mecrusldf[!, "ϵ"] = mecrusldf[!, "e_M_1_cm_1"] .* cm^-1 .* M^-1
	mecrusldf[!, "λ"] = mecrusldf[!, "λ_nm"] .* nm
	mecrusldf
end

# ╔═╡ e8aecd1c-a752-41f0-a536-37abecc9b906
md"""
## Interpolate function to data frame
"""

# ╔═╡ f6788e65-3b13-4982-a521-510970dac519
md"""
## Define fluorophore
"""

# ╔═╡ 463982fb-42cf-488d-ae12-b234abab2d12
md"""
## Obtain cross section
"""

# ╔═╡ c7154634-0a70-4eb6-a1d6-b36e122dabcb
rusl2 = jwf.jwfsim.Fluorophore(mecrusldf.λ[1], mecrusldf.λ[end], mecrusldf.ϵ, 0.11)

# ╔═╡ aaadb753-1124-44aa-a063-089d5c5a44a8
xsrusl = jwf.jwfsim.xsecfl(rusl2)

# ╔═╡ 9e476b6d-9a07-4d25-9444-43974e276c07
N_A

# ╔═╡ 4d6ffc87-f112-48d4-9eca-6bb19475f0a7
md"""
# Ir cross sections
"""

# ╔═╡ aca38e80-6a03-4ab7-95e7-586814848004
md"""
## Load Ir cross section data frames
"""

# ╔═╡ 8927a004-8324-4003-92e3-60d80862dc8c
mecirsl ="VS020_IrSL_1E-5_molar_extinction_coefficient.csv"

# ╔═╡ 9a1385d6-4775-4edc-b067-84e46ac0938e
begin
	mecirsldf = sort(jwf.jwfsim.load_df_from_csv(pdata, mecirsl, jwf.jwfsim.spG))
	mecirsldf[!, "ϵ"] = mecirsldf[!, "e_M_1_cm_1"] .* cm^-1 .* M^-1
	mecirsldf[!, "λ"] = mecirsldf[!, "λ_nm"] .* nm
	mecirsldf
end

# ╔═╡ 8fc41ed4-3c3c-4baf-9e72-525444b1979e
irsl = jwf.jwfsim.Fluorophore(mecirsldf.λ[1], mecirsldf.λ[end], mecirsldf.ϵ, 0.37)

# ╔═╡ 20de0098-563f-4e15-a17a-325a2d459c02
xsirsl = jwf.jwfsim.xsecfl(irsl)

# ╔═╡ ab5ce525-a5b1-4a2c-83b1-128dea1d0087
md"""
## Define fluorophore & obtain cross section
"""

# ╔═╡ a4f39abc-d187-48c4-91a2-cddaa3e5d0cd
begin
	plot(mecirsldf.λ, xsirsl.(mecirsldf.λ)*1e+16, label= "Xsec * 1e+16", lw=2)
	xlims!(200.,550.)
	
end

# ╔═╡ 020231a9-eb96-4a20-9d10-0a88d1568c81
xsirsl(325*nm)

# ╔═╡ 87536d90-a8d5-400c-9d68-e63ab1aaf05f
xsrusl(325*nm)

# ╔═╡ 6408aec3-2d07-406b-8235-8a0642b5d283
md"""
# Functions
"""

# ╔═╡ 10a5a5e3-ead5-414c-a10c-1b706825a48f
typeof(mecrusldf.ϵ[1])

# ╔═╡ 05dd6ada-4c43-4869-9d18-4f36fc5d3c87
typeof(mecrusldf.λ[end])

# ╔═╡ 6b09dab6-31ef-40c8-a36c-ed21888bdcc8
begin
	function icam(λ0::Real, λl::Real, vϵ::Vector{Float64})
		function gfpdf_(fi, xmin::Real, xmax::Real)
			function fn(λ::Unitful.Length)
				x = uconvert(NoUnits, λ/nm)
			   	if x < xmin || x > xmax
					return 0.0
				else
					return fi(x) * (cm^-1*M^-1)
				end
			end
			return fn
		end
		wl=λ0:λl
	    li = LinearInterpolation(wl, vϵ)
		return gfpdf_(li, λ0, λl)
	end
end

# ╔═╡ 36b52566-0fe0-499d-9ffa-0cae92e1b563
feps = icam(mecrusldf.λ[1]/nm, mecrusldf.λ[end]/nm, mecrusldf.ϵ ./(cm^-1*M^-1))

# ╔═╡ 170b34c1-f299-43a5-8f83-d8cb5586ad66
begin
	plot(mecrusldf.λ, mecrusldf.ϵ, label= "ϵ RuSL", lw=2)
	plot!(mecrusldf.λ, feps.(mecrusldf.λ), label= "interpolated", lw=2)
	xlims!(200.,550.)
end

# ╔═╡ 28662c26-b85b-407e-a9ad-7a6a75d49159
feps.(mecrusldf.λ)

# ╔═╡ b465a60e-db85-4423-9f0d-388eab0fe959
feps(325*nm)

# ╔═╡ af4da0ad-ae2e-44db-90a5-a10b28408db1
σ = log(10) * uconvert(cm^2/mol, feps(mecrusldf.λ[1])) / N_A

# ╔═╡ 5c12dfd4-b270-4acd-9255-ff69b252f7a5
struct Fluorophore
    λ0::Unitful.Length 
	λl::Unitful.Length
    ϵ::Vector{typeof(1.0/(cm*M))}
    Q::Float64
end

# ╔═╡ a9edc88c-ccd7-4f36-97ab-23c5fe82bc88
rusl = Fluorophore(mecrusldf.λ[1], mecrusldf.λ[end], mecrusldf.ϵ, 0.11)

# ╔═╡ cd727e71-8287-428d-bb5a-b13d164c6ff5
function xsecfl(fl::Fluorophore)
	function fn(λ::Unitful.Length)
		return log(10) * uconvert(cm^2/mol, fl.Q * feps(λ)) / N_A
	end

	feps = icam(fl.λ0/nm, fl.λl/nm, fl.ϵ ./(cm^-1*M^-1))
	return fn
	
end

# ╔═╡ 5dda25ae-fde8-43dc-8153-edac389f4d7d
xsfl = xsecfl(rusl)

# ╔═╡ 0ca757d4-c040-4722-bc1d-acb37e1de0bd
begin
	plot(mecrusldf.λ, xsfl.(mecrusldf.λ)*1e+16, label= "Xsec * 1e+16", lw=2)
	plot!(mecrusldf.λ, xsrusl.(mecrusldf.λ)*1e+16, label= "Xsec * 1e+16", lw=2)
	xlims!(200.,550.)
	
end

# ╔═╡ 8cf1ea7c-075e-4a57-ac14-243a04c8f0d2
xsfl(mecrusldf.λ[1])

# ╔═╡ Cell order:
# ╠═6a3aa124-9c93-11ed-0dbf-d142a80cde9f
# ╠═4e9f0e0d-ee89-40d2-b6df-6e7f0d0e7d6d
# ╠═d57a99a6-008d-4d00-8a00-ee1315b55445
# ╠═a35c53c3-a162-492d-a9ff-c0aa96429b55
# ╠═4f68d411-cb03-42ea-b9be-e37fdf273060
# ╠═8b2dd016-f4ba-4872-895f-cbd54d53252d
# ╠═605a3cb2-feed-42ae-9816-62f6f6bc7d62
# ╠═539b8d7b-e2f0-4a58-a868-03249ebb7ccd
# ╠═af94bd47-0817-48bf-b8e8-b53843ba7109
# ╠═f8e49a64-379f-4b23-9e87-4fe3a32e6c0e
# ╠═41e2bf4a-7543-4481-9f78-bff2629b7b55
# ╠═22edeacc-9aff-45df-9c61-750f8386b84f
# ╠═bb8b43d3-41c2-4c6f-999d-e1146d2b351b
# ╠═e8aecd1c-a752-41f0-a536-37abecc9b906
# ╠═36b52566-0fe0-499d-9ffa-0cae92e1b563
# ╠═170b34c1-f299-43a5-8f83-d8cb5586ad66
# ╠═28662c26-b85b-407e-a9ad-7a6a75d49159
# ╠═b465a60e-db85-4423-9f0d-388eab0fe959
# ╠═f6788e65-3b13-4982-a521-510970dac519
# ╠═a9edc88c-ccd7-4f36-97ab-23c5fe82bc88
# ╠═463982fb-42cf-488d-ae12-b234abab2d12
# ╠═5dda25ae-fde8-43dc-8153-edac389f4d7d
# ╠═c7154634-0a70-4eb6-a1d6-b36e122dabcb
# ╠═aaadb753-1124-44aa-a063-089d5c5a44a8
# ╠═0ca757d4-c040-4722-bc1d-acb37e1de0bd
# ╠═9e476b6d-9a07-4d25-9444-43974e276c07
# ╠═af4da0ad-ae2e-44db-90a5-a10b28408db1
# ╠═8cf1ea7c-075e-4a57-ac14-243a04c8f0d2
# ╠═4d6ffc87-f112-48d4-9eca-6bb19475f0a7
# ╠═aca38e80-6a03-4ab7-95e7-586814848004
# ╠═8927a004-8324-4003-92e3-60d80862dc8c
# ╠═9a1385d6-4775-4edc-b067-84e46ac0938e
# ╠═8fc41ed4-3c3c-4baf-9e72-525444b1979e
# ╠═20de0098-563f-4e15-a17a-325a2d459c02
# ╠═ab5ce525-a5b1-4a2c-83b1-128dea1d0087
# ╠═a4f39abc-d187-48c4-91a2-cddaa3e5d0cd
# ╠═020231a9-eb96-4a20-9d10-0a88d1568c81
# ╠═87536d90-a8d5-400c-9d68-e63ab1aaf05f
# ╠═6408aec3-2d07-406b-8235-8a0642b5d283
# ╠═10a5a5e3-ead5-414c-a10c-1b706825a48f
# ╠═05dd6ada-4c43-4869-9d18-4f36fc5d3c87
# ╠═6b09dab6-31ef-40c8-a36c-ed21888bdcc8
# ╠═cd727e71-8287-428d-bb5a-b13d164c6ff5
# ╠═5c12dfd4-b270-4acd-9255-ff69b252f7a5
