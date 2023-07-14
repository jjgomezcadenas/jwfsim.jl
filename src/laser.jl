using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
#using Interpolations

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    fs, ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W


	"""
	struct Fov

Represent a field of view

# Fields
- `d::Unitful.Length`  : diameter of Fov
- `z::Unitful.Length`  : thickness
- `a::Unitful.Area`    : area (computed)
- `v::Unitful.Volume`  : volume (computed)

"""
struct Fov
    d::Unitful.Length
    z::Unitful.Length
	a::Unitful.Area
    v::Unitful.Volume

	function Fov(d,z)
		a = π * (d/2.)^2
		v = a * z
		new(d,z,a,v)
	end
end


abstract type Laser end

"""
Simple representation of a continous laser

# Fields
- `λ::typeof(1.0nm)`  : Laser wavelength
- `P::typeof(1.0mW)`  : Power

"""
struct CLaser <: Laser
	λ::Unitful.Length
	P::Unitful.Power
end


"""
Simple representation of a pulsed laser

# Fields
- `f::typeof(1.0nm)`  : Laser frequency
- `τ::typeof(1.0mW)`  : time

"""
struct PLaser <: Laser
	λ::Unitful.Length
	P::Unitful.Power
    f::typeof(1.0*MHz)
    τ::typeof(1.0*fs)
end

    
"""
	Representation of a Gaussian laser defined by a laser, 
    the location of the waist (z0) and the waist radius (w0)  
	
	# Fields
	- `laser::Laser`           : A Laser  
	- `w0::Unitful.Length`     : radius of the waist  
	- `zr::Unitful.Length`     : z for which the intensity of the beam is reduced by half
	- `I0::typeof(1.0mW/cm^2)` : Intensity of the beam at (0,0) (computed from w0)
	- `ρ0::typeof(1.0Hz/cm^2)` : Photon density at r = 0, z= 0
	- `θ0::Real `              : divergence of the beam  (computed from w0)
	
	"""
	struct GaussianLaser 
		laser::Laser
		w0::Unitful.Length
		z0::Unitful.Length
		I0::typeof(1.0mW/cm^2)
		ρ0::typeof(1.0Hz/cm^2)
		θ0::Real 

		function glaser_(λ::Unitful.Length, P::Unitful.Power, w0::Unitful.Length)
			zr = uconvert(mm, π * w0^2/λ)
			I0 = uconvert(mW/cm^2, 2.0 * P/(π * w0^2))
			ρ0 = uconvert(Hz/cm^2, n_photons(λ, 2.0 * P)/( π * w0^2))
			θ0 = w0/zr
			return zr, I0, ρ0, θ0
		end

        function GaussianLaser(laser::Laser, w0::Unitful.Length)
			z0, I0, ρ0, θ0 = glaser_(laser.λ, laser.P, w0)
			new(laser, w0, z0, I0, ρ0, θ0)
		end

        function GaussianLaser(λ::Unitful.Length, P::Unitful.Power, 
			                   w0::Unitful.Length)
			z0, I0, ρ0, θ0 = glaser_(λ, P, w0)
			new(CLaser(λ, P), w0, z0, I0, ρ0, θ0)
		end
	end


"""
Rate of photons (number of photons per unit time) produced by a laser
# Fields
- `laser::Laser`     : Laser
- `t::Unitful.Time`  : Time in which target is illuminated

"""
function n_photons(laser::Laser)
	uconvert(Hz, laser.P / photon_energy(laser.λ))
end


"""
Rate of photons (number of photons per unit time) corresponding to a wavelength
λ and a power P
# Fields
- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power

"""
function n_photons(λ::Unitful.Length, p::Unitful.Power)
	uconvert(Hz, p / photon_energy(λ))
end

"""
	photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area)

number of photons per unit time per unit area

# Fields

- `λ::Unitful.Length` : photon wavelength
- `p::Unitful.Power`  : Power
- `a::Unitful.Area`   : Area

"""

function photon_density(λ::Unitful.Length, p::Unitful.Power, a::Unitful.Area)
	return n_photons(λ, p)/ a
end


"""
	photon_density(l::Laser, fov::Fov)

number of photons per unit time per unit area, in a Fov illuminated by a laser

# Fields

- `laser::Laser` : Laser
- `fov::Fov`     : Field of view

"""
function photon_density(laser::Laser, fov::Fov)
	return n_photons(laser) / fov.a
end


"""
Given wavelength of photon return its energy.
# Fields
- `λ::Unitful.Length`  : Photon wavelength

"""
function photon_energy(λ::Unitful.Length)
	uconvert(eV, λ, Spectral())
end


function energy_per_pulse(lsr::PLaser)
	lsr.P * lsr.f
end

"""
Returns the diffractive limit for a laser with wavelength λ 
filling the entrance pupil of an objective with NA.    

# Fields

- `λ::Float64`  : λ of laser beam
- `NA::Float64` : NA of the objective, used to focus the laser.

"""
diffractive_limit(λ::Float64, NA::Float64; ff=1.22) = ff * λ / (2 * NA)

diffractive_limit(λ::Unitful.Length, NA::Float64; ff=1.22) = ff * λ / (2 * NA)


"""
Depth of penetration assuming a gaussian laser that fills fully the
objective (e.g, diffractive limit)

"""
function diffractive_limit_zr(λ::Unitful.Length, NA::Float64)
	λ/(π * NA^2)
end


"""
Returns the beam width: ``W(z) = W_0 \\sqrt{1 + (z/z_0)^2}``

# Fields
	- `gl::GaussianLaser`       : A gaussian Laser 

"""
function W(gl::GaussianLaser)
	Wz(z::Unitful.Length) = gl.w0 * sqrt(1.0 + (z/gl.z0)^2)
	return Wz
end


"""
Returns the beam radius of curvature: ``R(z) = z ( 1 + (z_0/z)^2) ``

# Fields
- `gl::GaussianLaser`       : A gaussian Laser 

"""
function R(gl::GaussianLaser)
	Rz(z::Unitful.Length) = z * (1.0 + (gl.z0/z)^2)
	return Rz
end


"""
Returns the beam Intensity: ``I(\\rho, z) = I_0 ( W_0 / W(z))^2 \\exp{-2 \\rho^2/W^2(z)}  ``

# Fields
- `gl::GaussianLaser`       : A gaussian Laser 

"""
function I(gl::GaussianLaser)
	Wz = W(gl)
	Irz(ρ::Unitful.Length, z::Unitful.Length) = 2.0 * gl.laser.P / (π * Wz(z))^2 * exp(-2.0 * ρ^2/Wz(z)^2)
	#Irz(ρ::Unitful.Length, z::Unitful.Length) = gl.I0 * (gl.w0 / Wz(z))^2 * exp(-2.0 * ρ^2/Wz(z)^2)
	return Irz
end


"""
Returns 2*w0 that is the diameter of the beam waist, also called spot size   
 
# Fields

- `gl::GaussianLaser` : A gaussian laser 

"""
spot_size(gl::GaussianLaser)  = 2 * gl.w0


"""
Returns 2θ that is the beam angular divergence    
 
# Fields

- `gl::GaussianLaser` : A gaussian laser 

"""
angular_divergence(gl::GaussianLaser)  = 2 * gl.θ0


"""
Returns 2z0   
 
# Fields

- `gl::GaussianLaser` : A gaussian laser 

"""

depth_of_focus(gl::GaussianLaser)  = 2 * gl.z0


"""
Returns the waist of a gaussian laser with parameter w0 resulting 
from focusing a paralell laser beam (depth of focus of beam
much longer than focal length of lens) with a lens of
focal length f.

# Fields

- `glaser::GaussianLaser` : A  gaussian laser 
- `f::Unitful.length`     : focal length of focusing lens. 

Returns the beam radius (beam waist) when a 
a paralell laser beam (depth of focus of beam
much longer than focal length of lens) diameter D
and wavelength  λ is focused by a lens of focal distance f

# Fields
- `D::Unitful.Length`   : initial diameter of beam
- `f::Unitful.Length`   : focal length of focusing lens.
- `λ::Unitful.Length`   : wavelength of laser  

"""
w0f(glaser::GaussianLaser, obj::Objective) = 2.0 * glaser.laser.λ * obj.f /(π * glaser.w0) 
	
w0f(λ::Unitful.Length, D::Unitful.Length, 
    f::Unitful.Length) = 2.0 * λ * f/(π * D)


# Wide field setup

"""
Returns the beam radius (beam waist) in a wide field setup characterized by
a focusing lens of focal length f and an objective with focusing length fMO.
The incident beam of initial diameter d0 is expanded a factor G and 
focused by the focusing lens at a distance f. 
An infinity corrected microscope is then located at a distance FMO, 
to produce a paralel beam o waist w0wf

# Fields
- `d0::Unitful.Length`  : initial diameter of beam
- `f::Unitful.Length`   : focal length of focusing lens.
- `fMO::Unitful.Length` : focal length of Objective.
- `G::Float64`          : Beam expansion factor 
-`σ::Float64            : 1.22 if defining w0 in terms of RMS, 1.05 if WHM  

"""
function w0wf(d0::Unitful.Length, f::Unitful.Length, fMO::Unitful.Length, 
              G::Float64, σ::Float64=1.22) 
	σ * π * d0 * G * fMO/(4.0 * f) 
end
                     
                     