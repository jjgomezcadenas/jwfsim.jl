using Pkg; Pkg.activate(ENV["JWfSim"])

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

end

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M


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

jwf = ingredients("../src/jwfsim.jl")

PlutoUI.TableOfContents(title="Notebook title", indent=true)

md"""
## Load Ru cross section data frames
"""

pdata = joinpath(ENV["JWfSim"], "data")

mecrusl ="VS014_RuSL_1E-5_molar_extinction coefficient.csv"

begin
	mecrusldf = sort(jwf.jwfsim.load_df_from_csv(pdata, mecrusl, jwf.jwfsim.spG))
	names(mecrusldf)
	mecrusldf[!, "ϵ"] = mecrusldf[!, "e_M_1_cm_1"] .* cm^-1 .* M^-1
	mecrusldf[!, "λ"] = mecrusldf[!, "λ_nm"] .* nm
	mecrusldf
end