using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
#using Interpolations

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W


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
Given wavelength of photon return its energy.
# Fields
- `λ::Unitful.Length`  : Photon wavelength

"""
function photon_energy(λ::Unitful.Length)
	uconvert(eV, λ, Spectral())
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
	Irz(ρ::Unitful.Length, z::Unitful.Length) = gl.I0 * (gl.w0 / Wz(z))^2 * exp(-2.0 * ρ^2/Wz(z)^2)
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

