using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
using Interpolations

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
	A, N, mol, mmol, V, L, M

	import PhysicalConstants.CODATA2018: N_A


GM = 1e-50 * cm^2*cm^2*s

"""
Abstract type representing a Molecular Sensor
"""
abstract type MSensor end

"""	
Represent a fluorescent molecule

# Fields
- `λ0::Units (nm)`                : minimum λ 
- `λl::Units (nm)`                : maximum λ
- `ϵ::Units(``cm^{-1} M^{-1}``)`  : molar extinction coefficient
- `Q::Float64`                    : Quantum efficiency
- σ ::typeof(cm^4*s)              : TPA cross section
- β ::typeof(1.0cm/W)             : TPA absorption coefficient 


"""
struct Fluorophore <: MSensor
    λ0::Unitful.Length 
	λl::Unitful.Length
    ϵ::Vector{typeof(1.0/(cm*M))}
    Q::Float64
    σ ::typeof(1.0*cm^4*s)
    
    function Fluorophore(λ0::Unitful.Length, λl::Unitful.Length, 
                         ϵ::Vector{typeof(1.0/(cm*M))}, Q::Float64)
        new(λ0, λl, ϵ, Q, 1.0*GM)
    end 

    function Fluorophore(λ0::Unitful.Length, λl::Unitful.Length, 
                         ϵ::Vector{typeof(1.0/(cm*M))}, Q::Float64,
                         σ ::typeof(1.0*cm^4*s))
        new(λ0, λl, ϵ, Q, σ)
    end

end


"""
Returns a function defining the cross section (in cm^2/mol)
of the fluorophore. 

The function interpolates function icam to the fluorophore data
describing the molar extinction coefficient, the computes the cross section
transforming to proper units and multiplying by the quantum efficiency, 
"""
function xsecfl(fl::Fluorophore, istep::Real)
    function icam(λ0::Real, λl::Real, step::Real, vϵ::Vector{Float64})
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
        wl=λ0:step:λl
        li = LinearInterpolation(wl, vϵ)
        return gfpdf_(li, λ0, λl)
    end

	function xs(λ::Unitful.Length)
		return log(10) * uconvert(cm^2/mol, fl.Q * feps(λ)) / N_A
	end

	feps = icam(fl.λ0/nm, fl.λl/nm, istep, fl.ϵ ./(cm^-1*M^-1))
	return xs
end


"""

Simplified representation of a fluorescent molecule

# Fields
- `ϵ::Units(``cm^{-1} M^{-1}``)` : molar extinction coefficient
- `Q::Float64`                   : Quantum efficiency
"""
struct SFluorophore <: MSensor
    ϵ::typeof(1.0/(cm*M))
    Q::Float64
	σ::typeof(1.0cm^2)
	function Fluorophore(ex, en, ϵ, Q)
		σ = log(10) * uconvert(cm^2/mol, ϵ) / N_A
		new(ϵ, Q, σ)
	end
end


"""

Simplified representation of a phosporescent molecule

# Fields
- `ϵ::Units(``cm^{-1} M^{-1}``)` : molar extinction coefficient
- `Q::Float64`                   : Quantum efficiency
- `λ::Unitful.Length`            : Lifetime of triplet state
- `ϵλ::Float64`                  : Fraction of triplet
"""
struct SPhospho <: MSensor
    ϵ::typeof(1.0/(cm*M))
    Q::Float64
    λ::Unitful.Length
    ϵλ::Float64
	σ::typeof(1.0cm^2)
    
	function Fluorophore(ϵ, Q, λ, ϵλ)
		σ = log(10) * uconvert(cm^2/mol, ϵ) / N_A
		new(ϵ, Q, λ, ϵλ, σ)
	end
end


