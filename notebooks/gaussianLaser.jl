### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ afdebed6-9589-11ed-1225-8975381de797
using Pkg; Pkg.activate(ENV["JWfSim"])

# ╔═╡ 32b15fac-bc76-4e63-a128-4c3ba48271cd
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

# ╔═╡ 92abd9eb-ed5d-47fd-b07a-6b37a4cdf006
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 943fbe38-bf25-4b46-887e-f8ded565c509
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

# ╔═╡ ed784b77-1703-4725-a4a5-e5b284f1c3bd
jwf = ingredients("../src/jwfsim.jl")

# ╔═╡ 127b308d-3fe4-4c2f-afd7-ba0d76869302
PlutoUI.TableOfContents(title="Gaussian Beam", indent=true)

# ╔═╡ 33461aa7-2a43-46a7-97f9-8d2b7992d88c
md"### Gaussian beam formulas
$u(r,z) = A_0 \frac{w_0}{w(z)} e^{-\rho^2/w^2(z)} e^{-j\chi}$

where $\chi$ defines the beam phase

and 

$w(z) = w_0 \sqrt{1 + (\frac{z}{z_0})^2}$

$R(z) = z (1 + \frac{z}{z_0})^2$

$w_0 =\sqrt{\frac{\lambda z_0}{\pi}}$

thus

$z_0 = \frac{\pi w_o^2}{\lambda}$

The beam intensity 

$I(r,z) =|u(r,z)|^2 =  I_0 [\frac{w_0}{w(z)}]^2 e^{-2\rho^2/w^2(z)}$

$z_r = \pi  \frac{w_0^2}{\lambda}$
The formulas verify the following limits:
$I(r, 0) = I_0 e^{-2(r/w_0)^2}$

$I(0, z) = I_0 [\frac{w_0}{w(z)}]^2 = \frac{I_0}{1 + (z/z_0)^2}$

Which has maximum value ($I_0$) at $z=0$ and then decreases as z increases. For $z=z_0$ the intensity along the axis is reduced by half.

For $z << z_r, w(z) = w_0$, and:
$I(0, z) = I_0$
$I(r, z) = I_0 e^{-2(r/w_0)^2}$

Thus, the beam intensity falls exponentially with the radius r and is a constant ($I_0$) along the longitudinal axis, for $r=0$.

The beam power is obtained integrating the irradiance at any transverse plane. For simplicity we can choose $z=0$:
$P = \int_0^\infty I(r,0) 2 \pi r dr = 2 \pi I_0 \int e^{-2(r/w_0)^2} r dr$ 

Define $y = \frac{2 r^2}{w_0^2}$. Then $dy = \frac{4 r}{w_0^2} dr$, $dr = \frac{w_0^2}{4 r} dy$ and:

$P = 2 \pi I_0 \frac{w_0^2}{4} \int_0^\infty e^{-y} dy$ 

Since $\int_0^\infty e^-y dy = 1$, we find:

$P = I_0 \frac{\pi w_0^2}{2}$

or:

$I_0 = \frac{2 P}{\pi w_0^2}$ 

and:

The beam intensity 

$I(r,z) =  \frac{2 P}{\pi w^2(z)} e^{-2\rho^2/w^2(z)}$
"

# ╔═╡ 2f2ba87d-c436-4904-9f57-5d493db0a208
md"""
## Example 1
- Consider a 1 mW He-Ne laser beam, with λ= 633 nm and spot size  $2W_0=0.1$mm.
- a) Determine the angular divergence of the beam, its depth of focus and is diameter at $z=3.5 \times 10^5$km (e.g., distance Earth to the moon)
- b) What if the radius of curvature of the wavefront at $z=0$, $z= z_0$ and $z= 2z_0$?
- c) What is the optical intensity (in mW/cm$^2$) at the beam center ($z=0, \rho=0$) and at $z= z_0$? Compare this with the intensity of a point-like source located at $z=0$  of a 100 W which produces a isotropical spherical wave. 
"""

# ╔═╡ 722f6bb8-3aa5-492e-9edc-b6dc1bd07d0b
md"""
- First define a laser of 633 nm and 1 mW power
"""

# ╔═╡ 1dd58aa6-e9cd-43e7-aa17-11c29ba10bb8
lsrhn = jwf.jwfsim.CLaser(633*nm, 1*mW)

# ╔═╡ eb2c0698-de01-4083-a233-0544f7d2a5b5
md"""
- Now define a gaussian laser of $W_0=500 \mu m$. 
"""

# ╔═╡ 29f89f52-a7a9-46ac-b233-e4509d346708
glsrhn = jwf.jwfsim.GaussianLaser(lsrhn, 500.0*μm)

# ╔═╡ f00bd514-46ca-4597-ba38-1b4baa85d8d7
md"""
- The angular divergence of the beam is $(jwf.jwfsim.angular_divergence(glsrhn)) rads
"""

# ╔═╡ 98aefa58-42a5-4879-a59e-39c3b5c8bfa0
md"""
- depth of focus is $(jwf.jwfsim.depth_of_focus(glsrhn))
"""

# ╔═╡ 9e6f7132-f796-400c-80a1-53a3bd903af1
md"""
Beam radius as a function of z
"""

# ╔═╡ 375681c3-f909-4d28-b4f3-60b210d15fe4
brz = jwf.jwfsim.W(glsrhn)

# ╔═╡ 9d4aa699-cc5c-4ce8-a61e-4815d5c8a082
md"""
- beam diameter at z = 0 => $(brz(0.0*km))
- beam diameter at the moon => $(uconvert(km, brz(3.5e+5*km)))
"""

# ╔═╡ d263d3ed-c631-441d-8021-3c75f8228fcb
md"""

- Radius of curvature
"""

# ╔═╡ 6e4fd8b8-8eed-4cb8-80f4-fcf7eef564db
rz = jwf.jwfsim.R(glsrhn)

# ╔═╡ fbdc97ec-beac-4edc-964e-5297dced5519
md"""

- Radius of curvature a z = 0 => $(rz(0.0*cm))
- Radius of curvature a z = z0 => $(rz(glsrhn.z0))
- Radius of curvature a z = 2z0 => $(rz(2*glsrhn.z0))
"""

# ╔═╡ 8b0bf8c4-bef4-4f7a-80d1-5d3db768a83c
md"""
- Optical intensity
"""

# ╔═╡ 9ef8b4ac-c34a-4aa1-905e-f376d1b7213d
bI = jwf.jwfsim.I(glsrhn)

# ╔═╡ d80b84e2-8220-4361-bff5-9717708c2977
begin

X1 = range(-2*glsrhn.w0, 2*glsrhn.w0, length=100)
Y1 = range(-2*glsrhn.w0, 2*glsrhn.w0, length=100)

G = [bI(sqrt(x^2 + y^2), 0.0*μm) for y in Y1, x in X1] # Note x-y "for" ordering
p01=plot(X1,Y1,G,st=:surface,xlabel="x",ylabel="y")

p11=contourf(X1, Y1, G, xlabel="x",ylabel="y")
plot(p01,p11,size=(800,400))
end

# ╔═╡ 78f39ebd-01b6-46ba-b23d-6724f7d58140
md"""

- Optical intensity at ρ=0 and z = 0 => $(bI(0.0*cm, 0.0*cm))
- Optical intensity at ρ=0 and z = z0 => $(bI(0.0*cm, glsrhn.z0))
"""

# ╔═╡ c2154ed9-cac7-4729-94bb-6511962f9213
begin

Z1 = range(-3*glsrhn.z0, 3*glsrhn.z0, length=100)

G1 = [bI(0.0*μm, z) for z in Z1] 
plot(Z1,G1, xlabel="z",ylabel="I", lw =2, size=(800,400))
end

# ╔═╡ b1b06dcb-2379-4699-93d4-b565079b9e00
md"""
## Transporting a gaussian beam
"""

# ╔═╡ a4297af5-da85-4feb-b426-665033451cce
md"""
### Focusing length well outside the depth of focus

$z-f >> z_0$

- Lens of magnification *M* and focal distance *f*
$w_0' = M w_0$

$\frac{1}{z'} + \frac{1}{z} = \frac{1}{f}$

$M = |\frac{f}{z-f}|$
"""

# ╔═╡ 58a97035-f56e-4696-a2e0-b98872b6b4d4
load("gaussianBeamFillingLens.png")

# ╔═╡ e7c08954-20bd-4560-bc2a-c6521630221a
md"""
### Waist fills the focusing lens

$w_0' = \frac{\lambda}{\pi w_0}$

$z' = f$
"""

# ╔═╡ c24a184e-ef7a-4e9d-8e4e-68c62abb84f2
load("gaussianBeamFillingLens2.png")

# ╔═╡ 7d7c53ab-5349-443e-834f-42d5905a2fc7
md"""
### Wide field setup

1. Use a focusing length to focus laser into a point at f'
2. Plase MO at a distance $f_{MO}$ to achieve a paralell beam of diameter $d_0$

"""

# ╔═╡ 0e5da5ee-d0d9-4a66-81d0-8aea06188669
load("WideField.png")

# ╔═╡ ca3d0f9e-331a-4076-9d6d-bc42cd8877c6
md"""
### Wide field setup

$w_0 = \frac{2\lambda}{\pi D}$
$d_0 = 2 w_0' = \frac{1.22\lambda f_{MO}}{w_0}$
$w_0' =\frac{1.22\pi D}{4} \frac{f_{MO}}{f'}$
"""

# ╔═╡ f36c1871-8d07-465f-b563-63e8f249e4b9
md"""
## TOPATU WFS

- Laser: $\lambda=405$ nm; Variable power between 10 mW and 100 mW. 
- Setup 
 
 $f'=75-200$mm 

 $d_{beam} = 0.94$mm, 

 $D = G \times d_{beam}$, where G ϵ[2,8] 
 
 MO: 100x, 50x, 20x. Assume $f_{MO} = \frac{f_{TL}}{M} = \frac{200 mm}{M}$

- Considering FWHM replace 1.22 by 1.05 then:

$w_0' =\frac{1.05\pi \times 0.94 G}{4} \frac{f_{TL}}{M f'} \sim \frac{50 \pi G}{M f'}$

where the dimensions of $w_0'$ are the same than those of $d_{beam}$, e.g., mm
"""

# ╔═╡ d5eee90c-42b2-47a8-8289-44f2e34c05bd
w0wf(G::Float64, M::Float64, fmm::Float64) = 50 * π *G * mm/(M* fmm)

# ╔═╡ 14aad59d-a04c-4b93-85c4-f103d58055e6
md"""
- Wide field for M = 100
"""

# ╔═╡ 9ccec6ab-0c97-4271-b87e-8327f800c1fa
begin

GG = range(0,10, length=100)
WW75 = [w0wf(g, 100.0, 75.0) for g in GG] 
WW100 = [w0wf(g, 100.0, 100.0) for g in GG]
WW150 = [w0wf(g, 100.0, 150.0) for g in GG]
plot(GG,uconvert.(μm, WW75), xlabel="G",ylabel="w0'", lw =2, label="f'=75", size=(800,400))
plot!(GG,uconvert.(μm, WW100), xlabel="G",ylabel="w0'", lw =2, label="f'=100", size=(800,400))
plot!(GG,uconvert.(μm, WW150), xlabel="G",ylabel="w0'", lw =2, label="f'=150", size=(800,400))
end

# ╔═╡ f7a6f0bb-34fd-4f73-88d5-6d263a4d54df
md"""
- Wide field for M = 50
"""

# ╔═╡ b1240822-5295-4b4e-9d96-5cce44d30578
begin


W2W75 = [w0wf(g, 50.0, 75.0) for g in GG] 
W2W100 = [w0wf(g, 50.0, 100.0) for g in GG]
W2W150 = [w0wf(g, 50.0, 150.0) for g in GG]
plot(GG,uconvert.(μm, W2W75), xlabel="G",ylabel="w0'", lw =2, label="f'=75", size=(800,400))
plot!(GG,uconvert.(μm, W2W100), xlabel="G",ylabel="w0'", lw =2, label="f'=100", size=(800,400))
plot!(GG,uconvert.(μm, W2W150), xlabel="G",ylabel="w0'", lw =2, label="f'=150", size=(800,400))
end

# ╔═╡ c8765711-9d04-42c8-83a4-4daf3e9a6f36
md"""
### Gaussian beam power for M = 100, f' = 75, G=5 
"""

# ╔═╡ 745bd4ee-716b-4be3-90dd-837c1e211563
w0_m100_f75_g5 = uconvert(μm, w0wf(5.0, 100.0, 75.0))

# ╔═╡ dbdf2cfa-657a-48b6-8b1c-d5dd70b4aaea
lwf = jwf.jwfsim.CLaser(405*nm, 10*mW)

# ╔═╡ fa46a4b2-4165-4fc6-93b4-29d03bdcd91e
glwf_w0_m100_f75_g5 = jwf.jwfsim.GaussianLaser(lwf, w0_m100_f75_g5)

# ╔═╡ dcbdbfb2-d0eb-434d-8b83-2e220c7b93df
md"""
### Gaussian beam power for M = 50, f' = 75, G=8 
"""

# ╔═╡ 35f06064-a1f2-47bc-a23e-020d0a64a965
w0_m50_f75_g8 = uconvert(μm, w0wf(8.0, 50.0, 75.0))

# ╔═╡ 68405a48-a910-4df3-b79a-52ccf2e3e7d5
glwf_w0_m50_f75_g8 = jwf.jwfsim.GaussianLaser(lwf, w0_m50_f75_g8)

# ╔═╡ a2817440-1118-4760-a9f5-9ac4a672e079
md"""
#### Power ratio from w0'~100 μm to w0' ∼ 330 μm: 

$(glwf_w0_m50_f75_g8.ρ0/glwf_w0_m100_f75_g5.ρ0) 
"""

# ╔═╡ Cell order:
# ╠═afdebed6-9589-11ed-1225-8975381de797
# ╠═32b15fac-bc76-4e63-a128-4c3ba48271cd
# ╠═92abd9eb-ed5d-47fd-b07a-6b37a4cdf006
# ╠═943fbe38-bf25-4b46-887e-f8ded565c509
# ╠═ed784b77-1703-4725-a4a5-e5b284f1c3bd
# ╠═127b308d-3fe4-4c2f-afd7-ba0d76869302
# ╠═33461aa7-2a43-46a7-97f9-8d2b7992d88c
# ╠═2f2ba87d-c436-4904-9f57-5d493db0a208
# ╠═722f6bb8-3aa5-492e-9edc-b6dc1bd07d0b
# ╠═1dd58aa6-e9cd-43e7-aa17-11c29ba10bb8
# ╠═eb2c0698-de01-4083-a233-0544f7d2a5b5
# ╠═29f89f52-a7a9-46ac-b233-e4509d346708
# ╠═f00bd514-46ca-4597-ba38-1b4baa85d8d7
# ╠═98aefa58-42a5-4879-a59e-39c3b5c8bfa0
# ╠═9e6f7132-f796-400c-80a1-53a3bd903af1
# ╠═375681c3-f909-4d28-b4f3-60b210d15fe4
# ╠═9d4aa699-cc5c-4ce8-a61e-4815d5c8a082
# ╠═d263d3ed-c631-441d-8021-3c75f8228fcb
# ╠═6e4fd8b8-8eed-4cb8-80f4-fcf7eef564db
# ╠═fbdc97ec-beac-4edc-964e-5297dced5519
# ╠═8b0bf8c4-bef4-4f7a-80d1-5d3db768a83c
# ╠═9ef8b4ac-c34a-4aa1-905e-f376d1b7213d
# ╠═d80b84e2-8220-4361-bff5-9717708c2977
# ╠═78f39ebd-01b6-46ba-b23d-6724f7d58140
# ╠═c2154ed9-cac7-4729-94bb-6511962f9213
# ╠═b1b06dcb-2379-4699-93d4-b565079b9e00
# ╠═a4297af5-da85-4feb-b426-665033451cce
# ╠═58a97035-f56e-4696-a2e0-b98872b6b4d4
# ╠═e7c08954-20bd-4560-bc2a-c6521630221a
# ╠═c24a184e-ef7a-4e9d-8e4e-68c62abb84f2
# ╠═7d7c53ab-5349-443e-834f-42d5905a2fc7
# ╠═0e5da5ee-d0d9-4a66-81d0-8aea06188669
# ╠═ca3d0f9e-331a-4076-9d6d-bc42cd8877c6
# ╠═f36c1871-8d07-465f-b563-63e8f249e4b9
# ╠═d5eee90c-42b2-47a8-8289-44f2e34c05bd
# ╠═14aad59d-a04c-4b93-85c4-f103d58055e6
# ╠═9ccec6ab-0c97-4271-b87e-8327f800c1fa
# ╠═f7a6f0bb-34fd-4f73-88d5-6d263a4d54df
# ╠═b1240822-5295-4b4e-9d96-5cce44d30578
# ╠═c8765711-9d04-42c8-83a4-4daf3e9a6f36
# ╠═745bd4ee-716b-4be3-90dd-837c1e211563
# ╠═dbdf2cfa-657a-48b6-8b1c-d5dd70b4aaea
# ╠═fa46a4b2-4165-4fc6-93b4-29d03bdcd91e
# ╠═dcbdbfb2-d0eb-434d-8b83-2e220c7b93df
# ╠═35f06064-a1f2-47bc-a23e-020d0a64a965
# ╠═68405a48-a910-4df3-b79a-52ccf2e3e7d5
# ╠═a2817440-1118-4760-a9f5-9ac4a672e079
