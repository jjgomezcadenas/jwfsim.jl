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



# Wide field setup

"""
Returns the beam radius (beam waist) when a beam of diameter D
and wavelength  λ is focused by a lens of focal distance f

# Fields
- `D::Unitful.Length`   : initial diameter of beam
- `f::Unitful.Length`   : focal length of focusing lens.
- `λ::Unitful.Length`   : wavelength of laser  

"""
w0f(λ::Unitful.Length, D::Unitful.Length, f::Unitful.Length) = 2.0 * λ * f/(π * D)


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