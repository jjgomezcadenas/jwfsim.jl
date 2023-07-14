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
Simple representation of a microscope objective

# Fields
- `NA::Float64`          : Numerical aperture
- `M::Float64`           : Magnification
- `f::Unitful.Length`    : focal length 
- `d::Unitful.Length`    : entrance pupil diameter
- `wd::Unitful.Length`   : working distance
- `T::Float64`           : Transmittance

"""
struct Objective
    NA::Float64 
    M::Float64 
	f::Unitful.Length
	d::Unitful.Length
    wd::Unitful.Length
    T::Float64 

    function Objective(NA::Float64,  M::Float64, f::Unitful.Length, d::Unitful.Length,
                       wd::Unitful.Length, T::Float64)
		new(NA, M, f, d, wd,T)
	end

	function Objective(NA::Float64,  M::Float64)
		new(NA, M, -1.0mm, -1.0mm, -1.0mm, 1.0)
	end

	function Objective(f::Unitful.Length, d::Unitful.Length, M::Float64)
		ff = f / d   # https://www.eckop.com/resources/optics/numerical-aperture-and-f-number/
		NA = 1.0/(2.0*ff)
		new(NA, M, f,d, -1.0mm, 1.0)
	end
end

    """
    Defines a CCD
    
    # Fields
    
        -`QE`         : quantum efficiency
        -`dc`         : dark current (in rms e-)
        -`readout_n`  : readout noise (in rms e-)
        -`ntpxx`      : number of physical pixels in x
        -`ntpxy`      : number of physicalpixels in y
        -`stpxx`      : physical pixel size in x
        -`stpxx`      : physical pixel size in y
        -`bining`     : defines the software bins
        -`ssensx`     : physical size CCD along x
        -`ssensy`     : physical size CCD along y
        -`npxx`       : number of sotware pixels in x
        -`npxx`       : number of sotware pixels in y
        -`readout_nt` : readout noise per bin (in rms e-)
        -`dct`        : dark current per bin
    """
    struct CCD
        QE::Float64
        ntpxx::Integer 
        ntpxy::Integer 
        stpxx::Unitful.Length
        stpxy::Unitful.Length
        readout_n::Float64
        dc::typeof(1.0s^-1)
        binning::Integer
        npxx::Integer 
        npxy::Integer
        spxx::Unitful.Length
        spxy::Unitful.Length
        readout_nt::Float64
        dct::typeof(1.0s^-1)
        ssensx::Unitful.Length
        ssensy::Unitful.Length
        function CCD(QE, ntpxx,ntpxy, stpxx, stpxy, readout_n, dc, binning)
            npxx       = Int(ntpxx/sqrt(binning))
            npxy       = Int(ntpxy/sqrt(binning))
            spxx       = stpxx * sqrt(binning)
            spxy       = stpxy * sqrt(binning)
            ssensx     = stpxx * ntpxx
            ssensy     = stpxy * ntpxy
            readout_nt = readout_n * binning
            dct        = dc * binning
            new(QE, ntpxx, ntpxy, stpxx, stpxy, readout_n, dc, binning, 
                npxx, npxy, spxx, spxy, readout_nt, dct, ssensx, ssensy)
        end
    end



"""
Return the efficiency of a generic CCD as a function of wavelength.

# Fields

- `lmin::Float64=350.0` : Minimum wavelength for which efficiency is defined
- `lmax::Float64=350.0` : Maximum wavelength for which efficiency is defined

"""
function ϵccd(lmin::Unitful.Length=350.0nm, 
              lmax::Unitful.Length=1000.0nm)
	function eff(l::Unitful.Length)
        ll    = uconvert(NoUnits, l/nm)
        llmin = uconvert(NoUnits, lmin/nm)
        llmax = uconvert(NoUnits, lmax/nm)
		if ll < llmin || ll > llmax
			return 0.
		else
			wl = 350.:50.:1000.
			ϵ = [0.3, 0.4,0.65,0.78,0.82,0.82,0.8,0.72,0.62,0.5,0.37,
			  0.24,0.12,0.07]
			e = CubicSplineInterpolation(wl, ϵ)
			return e(ll)
		end
	end
	return eff
end


"""
Return the efficiency of a a counting photon APD as a function of wavelength.

# Fields

- `lmin::Float64=400.0` : Minimum wavelength for which efficiency is defined
- `lmax::Float64=1000.0` : Maximum wavelength for which efficiency is defined

"""
function ϵapd(lmin::Unitful.Length=400.0nm, 
              lmax::Unitful.Length=1000.0nm)
	function eff(l::Unitful.Length)
        ll    = uconvert(NoUnits, l/nm)
        llmin = uconvert(NoUnits, lmin/nm)
        llmax = uconvert(NoUnits, lmax/nm)
		if ll < llmin || ll > llmax
			return 0.
		else
			wl = 400.0:50.0:1000.0
			ϵ = [0.1, 0.3,0.4,0.6,0.65,0.7,0.7,0.65,0.62,0.5,0.4,
			  0.25,0.1]
			e = CubicSplineInterpolation(wl, ϵ)
			return e(ll)
		end
	end
	return eff
end


"""
Compute the fraction of photons that make it through an iris
of diameter D located at a distance d from the emission point.

# Fields

- `d::Float64`   : distance between emission point and iris
- `D::Float64`   : Diameter of iris

"""
function geometrical_acceptance(d::Float64, D::Float64)
	return 0.5(1. - d/sqrt(d^2 + (D/2.)^2))
end


"""
Compute the geometrical transmission of an objective (depends only of NA).

# Fields

- `objective::Objective` : Objective

"""
function geometrical_transmission(objective::Objective)
	A = objective.NA
	A <=1 ? (1 - sqrt(1 - A^2)) /2 : 0.5
end


"""
Compute the transmission of an objective (includes average transmission).

# Fields

- `objective::Objective` : Objective

"""
transmission(obj::Objective) = geometrical_transmission(obj) * obj.T^2


"""
Return the tranmission of the MUE12900 as function of wavelength.

# Fields

- `lmin::Float64=400.0` : Minimum wavelength for which efficiency is defined
- `lmax::Float64=1000.0` : Maximum wavelength for which efficiency is defined

"""
function tmue12900(lmin::Unitful.Length=300.0nm, 
                   lmax::Unitful.Length=850.0nm)
	function eff(l::Unitful.Length)
        ll    = uconvert(NoUnits, l/nm)
        llmin = uconvert(NoUnits, lmin/nm)
        llmax = uconvert(NoUnits, lmax/nm)
		if ll < llmin || ll > llmax
			return 0.
		elseif ll <= 440 && ll >= llmin
            wl = 300.0:20.0:440.0
            ϵ = [0.0, 0.03, 0.22, 0.64, 0.78, 0.82 ,0.84, 0.84]
			e = CubicSplineInterpolation(wl, ϵ)
            return e(ll)
        elseif ll > 440 && ll < 450
            return 0.85
        else
            wl = 450.0:50.0:850.0
            ϵ = [0.86, 0.87, 0.87, 0.86, 0.84, 0.80 ,0.84, 0.8, 0.60]
			e = CubicSplineInterpolation(wl, ϵ)
			return e(ll)
		end
	end
	return eff
end


"""
Return the tranmission of the MUE12900 as function of wavelength.

# Fields

- `lmin::Float64=400.0` : Minimum wavelength for which efficiency is defined
- `lmax::Float64=1000.0` : Maximum wavelength for which efficiency is defined

"""
function tLMM40VUV(lmin::Unitful.Length=300.0nm, 
                   lmax::Unitful.Length=850.0nm)
	function eff(l::Unitful.Length)
        ll    = uconvert(NoUnits, l/nm)
        llmin = uconvert(NoUnits, lmin/nm)
        llmax = uconvert(NoUnits, lmax/nm)
		if ll < llmin || ll > llmax
			return 0.
        else
            return 0.85
		end
	end
	return eff
end

### define some elements of the setup:

# EPL 375 laser

function epl375(f=1.0MHz)
    λ = 370.0nm
    fr = 20.0MHz
    pr = 150.0μW 
    P = pr/(fr/f)
    tau = 75.0ps
    PLaser(λ, P, f, tau)
end

#PQ LDHLDH375B

function pqLDH375B(f=1.0MHz)
    λ = 376.0nm
    fr = 80.0MHz
    pr = 7.4mW 
    P = pr/(fr/f)
    tau = 56.0ps
    PLaser(λ, P, f, tau)
end
