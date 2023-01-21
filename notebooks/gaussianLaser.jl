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
md"## Gaussian beam formulas
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
## Example 
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

$w_0 = \frac{2\lambda f}{\pi D}$
$d_0 = 2 w_0' = \frac{1.22\lambda f_{MO}}{w_0}$
$w_0' =\frac{1.22\pi D}{4} \frac{f_{MO}}{f'}$
$w_0' =\frac{1.05\pi D}{4} \frac{f_{MO}}{f'}$ (FWHM)
"""

# ╔═╡ f36c1871-8d07-465f-b563-63e8f249e4b9
md"""
## TOPATU WFS

- Laser: $\lambda=405$ nm; Variable power between 10 mW and 100 mW. 
- Setup 
$f'=75-200 ~{\rm mm}$
$d_{beam} = 0.94~{\rm mm}$ 
$D = G \times d_{beam}$ where G ϵ[2,8] 
 
 - For MUE31900 100x/0.6NA, which generates the widefield, $f_{MO}$ is 3mm

- MUE21200 20x/0.4NA, $f_{MO}$ is 10mm

- MUE20501 50x/0.55NA, $f_{MO}$ is 4mm

- The default value of f' is 150 mm. 
"""

# ╔═╡ 741243bd-f78a-46c9-95bf-f3b38a707248
md" set laser wavelength in nm $(@bind ll NumberField(10.0^2:10.0^3; default=405.0))"

# ╔═╡ 1cc13377-4b7e-479c-ae8e-3940b412a6c8
md" set laser power in mW $(@bind lp NumberField(10.0^1:10.0^3; default=10.0))"

# ╔═╡ ef4e3211-810a-4e33-aa2f-80af009af2c4
lsrwf = jwf.jwfsim.CLaser(ll*nm, lp*mW)

# ╔═╡ c2b0317b-2e55-456f-a523-71e12a28f205
md" set laser diameter in nm $(@bind ld NumberField(0.1:0.1:5.0; default=0.9))"

# ╔═╡ 0da8bc60-d543-4911-9584-536e7497e469
md"""
The laser is initially a gaussian beam of waist d0/2
"""

# ╔═╡ edf44dcb-48cd-4672-a241-1a932a47edcf
d0 = ld*mm/2.0

# ╔═╡ 164ca5d2-44f5-48f0-a075-8c8865b1d416
glsrwf1 = jwf.jwfsim.GaussianLaser(lsrwf, d0)

# ╔═╡ aa5579d3-d3c1-486d-a212-d58ece870db7
md"""
The laser is now expanded by a factor G and focused by a lens of focal length f
"""

# ╔═╡ fdca70a4-fd96-4527-b6b0-9a9b4e48061e
md" set laser expansion $(@bind le NumberField(1.0:0.1:10.0; default=8.0))"

# ╔═╡ 6b0890b7-0124-43b1-a0c6-cd582419d319
md" set focal length of focusing lens $(@bind fp NumberField(1.0:1.0:200.0; default=150.0))"

# ╔═╡ 2ad24978-223d-4c9c-9b78-739e617aac52
md""" 
The new waist at f is 
$w_0 = \frac{2\lambda f}{\pi G d_0}$
"""

# ╔═╡ 16d21d78-9182-4068-8308-735f7b7517e6
ff = fp * mm

# ╔═╡ 6a1fb576-0d8b-4fc0-b30d-067d1a2089f2
w0f = uconvert(μm, 2.0 * lsrwf.λ * ff/(π * le * d0))

# ╔═╡ 7890e3a4-a280-442b-aa87-ce69c42a5fc4
uconvert(μm, jwf.jwfsim.w0f(lsrwf.λ, le * d0, ff))

# ╔═╡ 84896adb-0a0d-48ff-b723-409b99b951f6
md"""
Finally a paralell beam is generated by using an infinity corrected objective
"""

# ╔═╡ ff1bf531-e1de-430a-aab0-41d050123fdd
schema = ["MUE31900", "MUE21200", "MUE20501"]

# ╔═╡ 6893c3d1-7b31-4090-981a-3eb0e0ba7fc4
md""" Select objective : $(@bind obj Select(schema))"""

# ╔═╡ 2cf0dfef-1158-4926-a3d5-ff07c823e173
FMO =Dict("MUE31900"=>3*mm, "MUE21200"=>10*mm, "MUE20501"=>4*mm)

# ╔═╡ 14aad59d-a04c-4b93-85c4-f103d58055e6
md"""
- Wide field waist for selected parameters (FWHM)
- Laser wavelength = $(lsrwf.λ)
- Laser diameter = $d0
- Laser expansion factor = $le
- Fousing length f' = $ff
- Objective = $obj
- Objective fMO = $(FMO[obj])
- w0f = $(uconvert(μm, jwf.jwfsim.w0f(lsrwf.λ, le * d0, ff)))
- w0fMO = $(uconvert(μm, jwf.jwfsim.w0wf(d0, ff, FMO[obj], le, 1.05)))
"""

# ╔═╡ 050e7fac-2867-4f69-a287-f4abae92a7ac
w0M0 = uconvert(μm, jwf.jwfsim.w0wf(d0, ff, FMO[obj], le, 1.05))
             

# ╔═╡ 9ccec6ab-0c97-4271-b87e-8327f800c1fa
begin

GG = range(0,10, length=10)
WW = [jwf.jwfsim.w0wf(d0, ff, FMO[obj], g, 1.05) for g in GG] 

plot(GG,uconvert.(μm, WW), xlabel="G",ylabel="w0MO", lw =2, label="f'=$ff, obj=$obj", size=(800,400))
end

# ╔═╡ e924e28a-0829-4209-a4d6-8a87d800b6eb
md"""
### Wide field gaussian beam
"""

# ╔═╡ 3b2b9d9e-9ae3-415d-9bb8-fdae2b87c084
gwf = jwf.jwfsim.GaussianLaser(lsrwf, w0M0)

# ╔═╡ 688f1328-90a8-4512-9349-8d18bc0bdd35
md"""
## Imaging Ru and Ir sensors
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
# ╠═741243bd-f78a-46c9-95bf-f3b38a707248
# ╠═1cc13377-4b7e-479c-ae8e-3940b412a6c8
# ╠═ef4e3211-810a-4e33-aa2f-80af009af2c4
# ╠═c2b0317b-2e55-456f-a523-71e12a28f205
# ╠═0da8bc60-d543-4911-9584-536e7497e469
# ╠═edf44dcb-48cd-4672-a241-1a932a47edcf
# ╠═164ca5d2-44f5-48f0-a075-8c8865b1d416
# ╠═aa5579d3-d3c1-486d-a212-d58ece870db7
# ╠═fdca70a4-fd96-4527-b6b0-9a9b4e48061e
# ╠═6b0890b7-0124-43b1-a0c6-cd582419d319
# ╠═2ad24978-223d-4c9c-9b78-739e617aac52
# ╠═16d21d78-9182-4068-8308-735f7b7517e6
# ╠═6a1fb576-0d8b-4fc0-b30d-067d1a2089f2
# ╠═7890e3a4-a280-442b-aa87-ce69c42a5fc4
# ╠═84896adb-0a0d-48ff-b723-409b99b951f6
# ╠═ff1bf531-e1de-430a-aab0-41d050123fdd
# ╠═6893c3d1-7b31-4090-981a-3eb0e0ba7fc4
# ╠═2cf0dfef-1158-4926-a3d5-ff07c823e173
# ╠═14aad59d-a04c-4b93-85c4-f103d58055e6
# ╠═050e7fac-2867-4f69-a287-f4abae92a7ac
# ╠═9ccec6ab-0c97-4271-b87e-8327f800c1fa
# ╠═e924e28a-0829-4209-a4d6-8a87d800b6eb
# ╠═3b2b9d9e-9ae3-415d-9bb8-fdae2b87c084
# ╠═688f1328-90a8-4512-9349-8d18bc0bdd35
