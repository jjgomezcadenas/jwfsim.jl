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
	stpxx=6.5μm
	stpxy=6.5μm
	binning=16
	ssensx=stpxx*ntpxx
	ssensy=stpxy*ntpxy
	npxx=Int(ntpxx/sqrt(binning))
	npxy=Int(ntpxy/sqrt(binning))
	spxx=stpxx*sqrt(binning)
	spxy=stpxy*sqrt(binning)
	QE=0.8
	readout_n=0.8
	dc=0.06s^-1
	readout_nt=readout_n*binning
	dct=dc*binning
	md"""
- Pixel array = $ntpxx x $ntpxy
- Pixel size = $stpxx x $stpxy 
- Sensor size= $ssensx x $ssensy 
- Binning= $binning
- Pixel array (with binning) = $npxx x $npxy
- Pixel size (with binning)= $spxx x $spxy 

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
	sbx=2.0mm
	sby=2.0mm
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
p = MvNormal(center,sigmas/mm)

X = range(-3*sigmas[1], 3*sigmas[1], length=100)
Y = range(-3*sigmas[1], 3*sigmas[1], length=100)
Z = [pdf(p, [x,y]) for y in Y/mm, x in X/mm] # Note x-y "for" ordering
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
md""" Select length (in mm) of (square) glass: $(@bind smm_0 NumberField(0.0:0.1:10.0, default=1))"""

# ╔═╡ 003fab4d-5181-4b93-a889-7594598e3ccf
smm=smm_0*mm

# ╔═╡ abe2369a-8a9a-4e0e-9707-4184e6d6f339
if sch=="Molecules per micron square"
	md""" Select the molecules per unit area (in molecules per $\mu$m$^2$): $(@bind a_0 NumberField(0.0:100.0, default=1.0))"""
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
	- Size of the FOV = $sfovx x $sfovy 
	- Resolution = $spxfovx x $spxfovy 
	"""
end

# ╔═╡ 245c0ff7-22b2-47aa-969d-5cf4eca9831e
md""" ## Molecule distribution in th FOV"""

# ╔═╡ a45829e4-ac16-4e5d-9ce2-6007d62f37e8
md""" ## Illumination in the FOV"""

# ╔═╡ 907729d2-6b55-4afb-a3fd-5a450dd57a63
begin
center1 = [npxx/2, npxy/2]
sigmas1 = [uconvert(NoUnits,sbx/spxfovx)/mag,uconvert(NoUnits,sby/spxfovy)/mag]
pp = MvNormal(center1,sigmas1)

X1 = range(0, npxx, length=npxx)
Y1 = range(0, npxy, length=npxy)
G = [pdf(pp, [x,y]) for y in Y1, x in X1] # Note x-y "for" ordering
p01=plot(X1,Y1,G,st=:surface,xlabel="x (mm)",ylabel="y (mm)")

p11=contourf(X1, Y1, G, xlabel="x (mm)",ylabel="y (mm)")
plot(p01,p11,size=(900,300))
end

# ╔═╡ adebe61f-6355-46f3-b6d8-86fd50484366
begin
sigmafovx=sigmas1[1]
sigmafovy=sigmas1[2]

md""" Sigma of the spot in the FOV is $sigmafovx microns in x axis and $sigmafovy microns in y axis"""
end

# ╔═╡ 2e1ff7fa-1f63-4d52-87e2-b0d61fc1f85d
heatmap(G)

# ╔═╡ c1af380d-8f25-4e3b-b5b3-5ac78f0a9913
md""" ## Molecule response
The number of emitted photons per second from the molecules will be $N_{emmitted}(x,y)=N_{incident}(x,y)N_{molecules}(x,y)\sigma$ where $N_{incident}$ is the number of photons hitting the sample at (x,y) position per unit of time, $N_{molecules}$ the number of molecules por unit area, and $\sigma$ the interaccion cross section. If the number of incident photons is too high there could be a saturation on the emitted photons.
"""

# ╔═╡ 090f9903-8eb8-4846-909e-2cd2ffeb536e
md""" - $N_{incident}(x,y)=\frac{P(x,y)}{e(\lambda)}=\frac{P_{total}G(x,y)}{e(\lambda)}$ where $P_{total}$ is the laser power, $e(\lambda)$ the energy per photon, and $G(x,y)$ the power distribution (in this case gaussian).
"""

# ╔═╡ f778f467-001d-4484-a20f-fd85e10558b1
md""" - $N_{molecules}(x,y)$ is the local molecule density. It's value will be computed in each pixel. Dividing the number of molecules in each pixel by its area.
"""

# ╔═╡ 58163cae-5a87-4de9-8747-f47dc0398389
md""" - $\sigma$ is unknown (as far as I know), but it should be a function of the photon energy.  
"""

# ╔═╡ b3c0f88a-074f-4ca6-8642-16801816cc36
md""" Select cross section -log$\sigma$ in (for $\sigma$ in cm$^2$): $(@bind r NumberField(10.0:20.0, default=16))"""

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

# ╔═╡ d32c1cc1-57b1-4d8e-a816-cf118e3353f1
begin
	sigma=(10^(-r))*cm^2
	md""" $\sigma$=$sigma"""
end

# ╔═╡ 1cc91044-ff75-4cae-9513-fa99b19a4697
md"""
### Incident photons
"""

# ╔═╡ cb384e7d-c506-4206-82f6-a34ac837fe66
begin
N_inc=uconvert.(s^-1,G*pwr/ei)
heatmap(N_inc,colorbar_exponentformat="power")
end

# ╔═╡ c063e645-f73d-41d0-9faa-253972d73f1c
md"""
### Emitted photons"""

# ╔═╡ 40645b82-f470-4cab-bb63-d775e289d8b0
md"""
## Diffraction limit
"""

# ╔═╡ 04e8c605-f597-46af-88ed-c3c154bee3e2
md"""
The maximum resolution that can be obtained with an opcical system is the diffraction limited resolution. In order to simulate that limit a gaussian filter has been aplied to the images. The sigma of that kernel must be similar to the diffraction limit. 
"""

# ╔═╡ cfa76021-fdd0-4719-8c40-0b4cac126239
dl=0.3*μm

# ╔═╡ 7d7fdb32-f273-42f4-ab16-a9ab3693fb97
md"""
## Camera response
"""

# ╔═╡ 5de5c815-f011-4dd0-9c18-6f769a7fd36a
md"""
The number of generated electrons will be $N_e=QE(\lambda)N_{em}T_{exp}+N_{dark}T_{exp}+N_{readout}$ where $QE$ is the quantum efficiency of the camera for the emission wavelength, $T_{exp}$ the exposition time, $N_{em}$ the emitted photons by the sample per second (shown above), $N_{dark}$ the dark current of the camera, and $N_{readout}$ the number of electrons generated in the readout process in each pixel.
"""

# ╔═╡ b95fed67-5021-4070-abef-d5150c6a82aa
md""" Set exposure time in seconds: $(@bind texp0 NumberField(0.0:20.0, default=5))"""

# ╔═╡ 8e6aac57-72cf-40d7-90de-12f5717e3420
texp=texp0*s

# ╔═╡ 3d7cd04d-89b8-4c58-9507-2b2477155504
md"""
# Simulate temporal evolution
"""

# ╔═╡ 7ae42d0e-190e-40bf-872e-29524b5a0dca
md""" To simulate the temporal evolution of the images we will simulate the temporal evolution of the molecules and the obtain the camera response of the corresponding molecule distribution. """

# ╔═╡ 445860e3-7dcf-4f0a-bd87-049e181a7c0d
md""" We will asume that each molecule has $\tau$ life-time for photobleaching. This means that after $t$, the molecule will have $P_{live}(t)=e^{-t/\tau}$ probability to be 'alive'. So the probability to die is $P_{die}(t)\equiv P_0=1-e^{-t/\tau}$"""

# ╔═╡ 531b42a1-a373-4e8a-9ed9-6228f07ecfbe
md""" In order to take into acount the effect of the illumination (the more photons hitting a molecule, the more likely is photobleaching to happen) we will scale the photobleaching probability with the illumination distribution $G$. This way, the probability for the photobleaching to happen will be $P_{die}(t,x,y)=\frac{G(x,y)}{G_{max}}(1-e^{-t/\tau})=\frac{G(x,y)}{G_{max}}P_0$"""

# ╔═╡ 9f89d2bf-74c3-4a1d-ab7f-64bfde1f4dc5
md""" Defining a frame rate of $1/\Delta t$ we can simulate the molecule evolution using Montecarlo methods."""

# ╔═╡ 656f048c-d0bb-4eb1-b2ea-b35ec1839da6
md""" Set life time $\tau$ in seconds: $(@bind tau0 NumberField(10:300.0, default=10))"""

# ╔═╡ 49ae9e64-7858-4bab-b89a-e6fb9028e168
tau=tau0*s

# ╔═╡ d5f267b2-6c6a-438c-97dd-75a05080cf87
md""" Set $\Delta t$ in miliseconds: $(@bind dt0 NumberField(1:1000.0, default=500))"""

# ╔═╡ 6c97d42b-f0e1-4e8d-a9cd-af6d79c8db72
dt=dt0*ms

# ╔═╡ a4528adf-f185-4bce-b695-12b0e1cecfab
prob0=1-exp(-uconvert(NoUnits,dt/tau))

# ╔═╡ 48ad5b24-bb6e-4628-8df0-0e00a4dc336f
md""" Set the number of frames $N_t$ : $(@bind Nt NumberField(1:200, default=100))"""

# ╔═╡ d5dfad8c-6b6c-4d08-93fc-7971c7bec339
md""" Select frame number : $(@bind Nframe NumberField(1:Nt, default=Nt))"""

# ╔═╡ 092ff59e-3aae-4dde-a778-973eb83acba1
md"""
# Data analysis
"""

# ╔═╡ 095c751f-b5de-42ed-b3cc-e9d518fe0bc8
md""" Select threshold : $(@bind thresh NumberField(1:1000, default=1000))"""

# ╔═╡ f6aa8eee-4537-4b0e-a6b0-6efd2b27c9f0
md""" Select cluster : $(@bind clustn NumberField(1:100,default=1))"""

# ╔═╡ c2289697-ab87-4f11-81e6-bee3439ca3cb
md"""
# Functions
"""

# ╔═╡ 32caa85b-f46a-4dc4-b3d8-d756bbd7e1f0
minimum(G)

# ╔═╡ 1e44b00e-70f5-4c46-a99e-939dbf5bf7f6
function evol(poss_array0,G,prob)
	g0=maximum(G)
	poss_array1=zeros(size(poss_array0))
	for i in range(1,size(poss_array0)[1])
		for j in range(1,size(poss_array0)[2])
			n=poss_array0[i,j]
			g=G[i,j]
			if n!=0
				count=0
				for m in range(1,n)
					r=rand()
					count+=1*(r>g/g0*prob0)
				end
				poss_array1[i,j]=count
			end
		end
	end
	poss_array1
end

# ╔═╡ 564e24e7-015f-4b5e-86a5-107909120b7f
function evol_data(poss_array0,G,prob,N)
	data=poss_array0
	prev=poss_array0
	for i in range(1,N-1)
		poss_array=evol(prev,G,prob)
		data=cat(data,poss_array,dims=3)
		prev=poss_array
	end
	data
end

# ╔═╡ ed8d60e5-fe91-4048-945f-692b791bbdae
function frame(poss_array,spxfovx,spxfovy,G,pwr,ei,sigma,dl,texp,dc,readout_n)
	N_inc=uconvert.(s^-1,G*pwr/ei)
	pxarea=uconvert(m^2,spxfovx*spxfovy)
	N_mol=poss_array/pxarea
	N_em=uconvert.(s^-1,N_inc.*N_mol*sigma)
	sigma_dl=[dl/spxfovx,dl/spxfovy]
	N_em_dl = imfilter(N_em,Kernel.gaussian((sigma_dl[1],sigma_dl[2])))
	N_e=QE*N_em_dl*texp.+(dct*texp.+readout_nt)*rand()
end

# ╔═╡ 7eada1cb-8890-4470-a9e2-f42620c9aa45
rand()

# ╔═╡ 8d920a11-2050-4562-b55c-0b916cabf42c
function toimage(data,xran0,yran0,spxfovx0,spxfovy0)
	u=unit(data[1,1])
	xran,yran,spxfovx,spxfovy=xran0/u,yran0/u,spxfovx0/u,spxfovy0/u
	x=data[1,:]/unit(data[1,1])
	y=data[2,:]/unit(data[1,1])
	fit(Histogram,(x, y), (xran[1]:spxfovx:xran[2], yran[1]:spxfovy:yran[2])).weights
end

# ╔═╡ 13b7f5db-8593-4511-b540-6f40e8eb3498
function poss_fov(poss0,xran0,yran0)
	poss=poss0/unit(poss0[1,1])
	xran=xran0/unit(xran0[1])
	yran=yran0/unit(yran0[1])
	possfov=[0,0]
	for pos in eachrow(poss)
		if (xran[1]<pos[1])&&(pos[1]<xran[2])
			if (yran[1]<pos[2])&&(pos[1]<yran[2])
				possfov=hcat(possfov,pos)
			end
		end
	end
	possfov=possfov[:,2:end]*unit(poss0[1,1])
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
	nmu2=a=a_0/μm/μm
end
if sch=="Molarity"
	nmu2=toncm2(b,smm/mm)/1e8/μm/μm
end
md""" Molecules per micron square (nmu2)=$nmu2 """
end


# ╔═╡ 84e73dfd-a9a3-4bd4-93ce-2fc4a4e78f1b
begin
smu=uconvert(μm,smm)
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
md""" Select x0 position on the sample (in microns). It ranges from 0 to $smu : $(@bind x0_0 NumberField(0:smu/μm, default=smu/2/μm))"""

# ╔═╡ 3388900d-6fa1-4c0b-86b0-9ca541ad39ed
x0=x0_0*μm

# ╔═╡ 50c17b83-57f5-48dc-aaa1-22e16f627022
md""" Select y0 position on the sample (in microns). It ranges from 0 to $smu : $(@bind y0_0 NumberField(0:smu/μm, default=smu/2/μm))"""

# ╔═╡ c4248750-7259-49ea-9aae-11ed7e3bf089
y0=y0_0*μm

# ╔═╡ 245205d3-ef94-46cf-ba0a-48810c201610
begin
xran=[x0-sfovx/2,x0+sfovx/2]
yran=[y0-sfovy/2,y0+sfovy/2]
poss
end

# ╔═╡ da45d069-b08d-4221-9e7d-32371090e94d
begin
possfov=poss_fov(poss, xran, yran)
poss_array=toimage(possfov,xran,yran,spxfovx,spxfovy)
heatmap(poss_array)
end

# ╔═╡ 1d5d2334-7973-48c3-a9c4-e0617379ec2f
begin
	pxarea=uconvert(m^2,spxfovx*spxfovy)
	N_mol=poss_array/pxarea
	N_em=uconvert.(s^-1,N_inc.*N_mol*sigma)
	heatmap(N_em)
end

# ╔═╡ 76820384-74f6-46d9-9b70-83d93f1856ed
begin
	sigma_dl=[dl/spxfovx,dl/spxfovy]
	N_em_dl = imfilter(N_em,Kernel.gaussian((sigma_dl[1],sigma_dl[2])))
	heatmap(N_em_dl)
end

# ╔═╡ 631fd633-0519-42b0-89a7-eae72bba519a
begin
	N_e=QE*N_em_dl*texp.+(dct*texp.+readout_nt)*rand()
	heatmap(N_e)
end

# ╔═╡ 1b13988e-8e3b-475e-8e08-8a92e4a638e4
begin
	datas=evol_data(poss_array,G,prob0,Nt)
	md""" The temporal evolution of the molecule distribution is stored in an $npxx x $npxy x $Nt array (named datas). Each [:,:,i] slide represents a molecule distribution."""
end

# ╔═╡ 646e52c4-5c4c-460c-be3d-b2382745ac89
begin
	datat=datas[:,:,Nframe]
	imt=frame(datat,spxfovx,spxfovy,G,pwr,ei,sigma,dl,texp,dc,readout_n)
	heatmap(imt)
end

# ╔═╡ ec93b3d6-bdef-4cb9-9be7-08693174beaf
begin
	moleculen=[]
	for i in range(1,Nt)
		append!(moleculen,sum(datas[:,:,i]))
	end
		
end

# ╔═╡ 420756fd-4f72-4eda-90be-e81939137677
plot(moleculen,xlabel="Frame",ylabel="Number of molecules in FOV")

# ╔═╡ 1c37e944-f18a-440d-ac01-bf8ecdb384c9
begin
	datat0=datas[:,:,1]
	imt0=frame(datat0,spxfovx,spxfovy,G,pwr,ei,sigma,dl,texp,dc,readout_n)
	threshim=Gray.(imt0) .>thresh
end

# ╔═╡ 1fa4b498-0c7a-404e-9760-d526a0a3d40a
begin
tt=Float64.(hcat(collect.(Tuple.(findall(x->x==1,threshim)))...))
clusters=dbscan(tt,1,min_cluster_size = 5)
end

# ╔═╡ e7a683e9-7545-4f57-835b-ea4a96320596

begin
	clusteri=clusters[clustn].core_indices
	clust=[tt[:,i] for i in clusteri]
	clustm=Int64.(transpose.(reduce(vcat,transpose.(clust))))
end

# ╔═╡ 05fe65a4-ebff-4f02-b717-a83ecdd64024
begin
	heatmap(threshim)
	scatter!([clustm[1,:][1]],[clustm[1,:][2]], markershape=:circle,label="ROI")
end

# ╔═╡ bf5b51d6-009b-4866-b489-2edd3d9a05a6
begin
	psincl=[imt0[row[1],row[2]] for row in eachrow(clustm)]
	scatter(clustm[:,1],clustm[:,2],marker_z=psincl)
end

# ╔═╡ 7e07dafb-2ff2-4a35-b9b9-e7bc84c1ccfa
begin
	psROI=[]
	signalROI=[]
	for i in range(1,Nt)
		datatt=datas[:,:,i]
		imtt=frame(datatt,spxfovx,spxfovy,G,pwr,ei,sigma,dl,texp,dc,readout_n)
		psinclt=[datatt[row[1],row[2]] for row in eachrow(clustm)]
		siginclt=[imtt[row[1],row[2]] for row in eachrow(clustm)]
		append!(psROI,sum(psinclt))
		append!(signalROI,sum(siginclt))
	end
end

# ╔═╡ 61700541-c9c2-41d0-beb1-3dd25351313c
plot(signalROI,xlabel="Frame",ylabel="Signal in ROI")

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
# ╠═cc82232e-fc86-413f-8b73-8044b4de6260
# ╟─790af81b-45ce-4fae-ab0b-7fb52b87a2af
# ╟─018997ec-0c72-4677-b6d1-e8dafa1c2811
# ╟─0782a2c9-7d5a-49bc-8ddb-71c57238de8b
# ╟─55352380-ff6d-42f6-863a-e5de8cced158
# ╟─003fab4d-5181-4b93-a889-7594598e3ccf
# ╠═abe2369a-8a9a-4e0e-9707-4184e6d6f339
# ╠═00cfc948-31e4-4d17-bf44-0a25c911cf68
# ╠═18fc7474-c7b8-4092-9dfa-92f07ce3cb2a
# ╟─39ec5a5f-8e06-4143-8016-39fc9020aaa1
# ╠═84e73dfd-a9a3-4bd4-93ce-2fc4a4e78f1b
# ╠═b32cd1dd-a93a-41e6-a9cd-92bb9295f3ad
# ╟─766d17b5-499a-4a28-8628-3031b023a851
# ╟─a24aab57-2145-40de-bfd0-b768fef6d01d
# ╠═f7f4ba8b-e1a6-4632-a785-767b2bfb646f
# ╠═d9e3e919-f1b6-4b0b-9e07-04d5855f3587
# ╟─3388900d-6fa1-4c0b-86b0-9ca541ad39ed
# ╠═50c17b83-57f5-48dc-aaa1-22e16f627022
# ╠═c4248750-7259-49ea-9aae-11ed7e3bf089
# ╠═245205d3-ef94-46cf-ba0a-48810c201610
# ╠═245c0ff7-22b2-47aa-969d-5cf4eca9831e
# ╠═da45d069-b08d-4221-9e7d-32371090e94d
# ╠═a45829e4-ac16-4e5d-9ce2-6007d62f37e8
# ╠═907729d2-6b55-4afb-a3fd-5a450dd57a63
# ╟─adebe61f-6355-46f3-b6d8-86fd50484366
# ╟─2e1ff7fa-1f63-4d52-87e2-b0d61fc1f85d
# ╠═c1af380d-8f25-4e3b-b5b3-5ac78f0a9913
# ╟─090f9903-8eb8-4846-909e-2cd2ffeb536e
# ╟─f778f467-001d-4484-a20f-fd85e10558b1
# ╟─58163cae-5a87-4de9-8747-f47dc0398389
# ╠═b3c0f88a-074f-4ca6-8642-16801816cc36
# ╟─8c9160ca-0048-4de6-8a39-ad00fe44f680
# ╠═f07e7613-1402-4abc-8b06-dc98df2b8c98
# ╠═bd7a2018-d882-4891-8f15-d90e92cb1077
# ╠═d32c1cc1-57b1-4d8e-a816-cf118e3353f1
# ╠═1cc91044-ff75-4cae-9513-fa99b19a4697
# ╠═cb384e7d-c506-4206-82f6-a34ac837fe66
# ╠═c063e645-f73d-41d0-9faa-253972d73f1c
# ╠═1d5d2334-7973-48c3-a9c4-e0617379ec2f
# ╠═40645b82-f470-4cab-bb63-d775e289d8b0
# ╟─04e8c605-f597-46af-88ed-c3c154bee3e2
# ╠═cfa76021-fdd0-4719-8c40-0b4cac126239
# ╠═76820384-74f6-46d9-9b70-83d93f1856ed
# ╟─7d7fdb32-f273-42f4-ab16-a9ab3693fb97
# ╟─5de5c815-f011-4dd0-9c18-6f769a7fd36a
# ╟─b95fed67-5021-4070-abef-d5150c6a82aa
# ╟─8e6aac57-72cf-40d7-90de-12f5717e3420
# ╠═631fd633-0519-42b0-89a7-eae72bba519a
# ╠═3d7cd04d-89b8-4c58-9507-2b2477155504
# ╟─7ae42d0e-190e-40bf-872e-29524b5a0dca
# ╟─445860e3-7dcf-4f0a-bd87-049e181a7c0d
# ╟─531b42a1-a373-4e8a-9ed9-6228f07ecfbe
# ╟─9f89d2bf-74c3-4a1d-ab7f-64bfde1f4dc5
# ╟─656f048c-d0bb-4eb1-b2ea-b35ec1839da6
# ╟─49ae9e64-7858-4bab-b89a-e6fb9028e168
# ╟─d5f267b2-6c6a-438c-97dd-75a05080cf87
# ╟─6c97d42b-f0e1-4e8d-a9cd-af6d79c8db72
# ╟─a4528adf-f185-4bce-b695-12b0e1cecfab
# ╟─48ad5b24-bb6e-4628-8df0-0e00a4dc336f
# ╟─1b13988e-8e3b-475e-8e08-8a92e4a638e4
# ╟─d5dfad8c-6b6c-4d08-93fc-7971c7bec339
# ╟─646e52c4-5c4c-460c-be3d-b2382745ac89
# ╟─ec93b3d6-bdef-4cb9-9be7-08693174beaf
# ╟─420756fd-4f72-4eda-90be-e81939137677
# ╠═092ff59e-3aae-4dde-a778-973eb83acba1
# ╠═095c751f-b5de-42ed-b3cc-e9d518fe0bc8
# ╠═1c37e944-f18a-440d-ac01-bf8ecdb384c9
# ╠═1fa4b498-0c7a-404e-9760-d526a0a3d40a
# ╠═f6aa8eee-4537-4b0e-a6b0-6efd2b27c9f0
# ╠═e7a683e9-7545-4f57-835b-ea4a96320596
# ╠═05fe65a4-ebff-4f02-b717-a83ecdd64024
# ╠═bf5b51d6-009b-4866-b489-2edd3d9a05a6
# ╠═7e07dafb-2ff2-4a35-b9b9-e7bc84c1ccfa
# ╠═61700541-c9c2-41d0-beb1-3dd25351313c
# ╠═c2289697-ab87-4f11-81e6-bee3439ca3cb
# ╠═564e24e7-015f-4b5e-86a5-107909120b7f
# ╠═32caa85b-f46a-4dc4-b3d8-d756bbd7e1f0
# ╠═1e44b00e-70f5-4c46-a99e-939dbf5bf7f6
# ╠═ed8d60e5-fe91-4048-945f-692b791bbdae
# ╠═7eada1cb-8890-4470-a9e2-f42620c9aa45
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
