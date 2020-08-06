import numpy 
import matplotlib 
from matplotlib import pylab
matplotlib.rcParams.update({'font.size':18, 'font.family':'DejaVu Sans'})
from scipy import signal



def Faraday2Lambda(lam2, phi_range, pol_lam2):

    """
    Computes Faraday Spectra using RM-Synthesis 
    as defined by Brentjens and de Bruyn (2005)

    """

    N = len(lam2)
    l20 = lam2.mean()
    fdata = numpy.zeros([len(phi_range)], dtype=complex)
    for k, phi in enumerate(phi_range):
        fdata[k] = pow(N, -1) * numpy.sum(pol_lam2 * 
                numpy.exp(-2j * (lam2-l20) * phi))    
    return fdata


def RMClean(lam2, phi_range, fspectrum, 
          niter=500, gain=0.1):

    fwhm = (3.8/ abs(lam2[0]-lam2[-1]))
    sigma = (fwhm/2.35482)
    Gauss = numpy.exp(-0.5 * (phi_range/sigma)**2) 

    # I am padding here to avoid edge effects.
    
    pad = abs(phi_range[-1]) * 2
    dpad = abs(phi_range[0]-phi_range[1])
    phi_pad = numpy.arange(-pad, pad, dpad)
    dshift = int(pad/(2.0 * dpad))

    rmsf_fixed = Faraday2Lambda(lam2, phi_pad, 1) 
    components = numpy.zeros([len(phi_range)], dtype=complex)

    for n in range(niter):
        temp = numpy.zeros([len(phi_range)], dtype=complex)
        f_amp = numpy.absolute(fspectrum)
        ind = numpy.where(f_amp == f_amp.max())[0]
        f_comp = fspectrum[ind[0]] * gain
        temp[ind[0]] = f_comp
        components += temp         
    
        dirac = numpy.zeros(len(phi_range))
        dirac[ind[0]] = 1
        rmsf = signal.convolve(rmsf_fixed, dirac, mode='same')
        rmsf = rmsf[dshift:-dshift+1]
 
        fspectrum -= f_comp * rmsf

    Fres = fspectrum
    fclean = signal.convolve(components, Gauss, mode='same') + Fres

    return fclean, components


def plot_data(lam, plam, phi, fphi):
    
    
    fig, (ax, ay) = pylab.subplots(1, 2, figsize=(12, 5))
    ax.plot(lam, abs(plam), 'b.')
    ax.set_xlabel('Wavelength [m$^2$]')
    ax.set_ylabel('Fractional Polarization')
          
    ay.plot(phi, abs(fphi), 'b')
    ay.set_xlabel('Faraday Depth [rad m$^{-2}$]')
    ay.set_ylabel('Faraday Spectrum')
    pylab.tight_layout()
    pylab.show()



def main(lam, pdata, phi_max=5000, phi_step=10,
         niter=1000, gain=0.1, plot=False):


    phi_range =  numpy.arange(-phi_max, phi_max+phi_step, phi_step)
    # this ensures that the middle value is zero. 
    fdirty = Faraday2Lambda(lam, phi_range, pdata)
    
    
    fclean, fcomp = RMClean(lam, phi_range, fdirty, 
                    niter=niter, gain=gain)

    if plot:
        plot_data(x, pdata, phi_range, fclean)


    return phi_range, fclean, fcomp


