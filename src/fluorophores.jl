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


"""	
Represent a fluorescent molecule

# Fields
- `λ0::Units (nm)`           : minimum λ 
- `λl::Units (nm)`           : maximum λ
- `ϵ::Units(``cm^{-1} M^{-1}``)` : molar extinction coefficient
- `Q::Float64`                   : Quantum efficiency

"""
struct Fluorophore
    λ0::Unitful.Length 
	λl::Unitful.Length
    ϵ::Vector{typeof(1.0/(cm*M))}
    Q::Float64
end


"""
Returns a function defining the cross section (in cm^2/mol)
of the fluorophore. 

The function interpolates function icam to the fluorophore data
describing the molar extinction coefficient, the computes the cross section
transforming to proper units and multiplying by the quantum efficiency, 
"""
function xsecfl(fl::Fluorophore)
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

	function xs(λ::Unitful.Length)
		return log(10) * uconvert(cm^2/mol, fl.Q * feps(λ)) / N_A
	end

	feps = icam(fl.λ0/nm, fl.λl/nm, fl.ϵ ./(cm^-1*M^-1))
	return xs
end

