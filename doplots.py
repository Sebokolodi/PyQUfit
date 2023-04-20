import rmsynthesis as rmsyn
import matplotlib
matplotlib.rcParams.update({'font.size':16})
import numpy
import corner
import pylab



def remove_ambiguity(PA, RM):

    """Corrects for npi-ambiquity"""

    difference = abs(numpy.diff(PA))
    ind = numpy.where(difference >= 1.0)[0] + 1
    pa_split = numpy.split(PA, ind)
    n = 0

    while len(ind) > 0:
        
        if RM > 0:
            pa_split[-2] = pa_split[-2] + numpy.pi/2.0
        
        if RM < 0:
            pa_split[-2] = pa_split[-2] - numpy.pi/2.0
            
        pa = numpy.hstack(pa_split)
        difference = abs(numpy.diff(pa))
        ind = numpy.where(difference >= 1.0)[0] + 1
        pa_split = numpy.split(pa, ind)
        n += 1
        if n > 3000:
            print("Couldn't solve for npi-amb. Returning the original PA." )
            return PA
            break   
    pa  = numpy.hstack(pa_split) 

    return pa



def do_plots(x, qdata, udata, sigmap, sigmaPA, 
                qmodel, umodel, nparams, phi_range,
                weighted_posteriors, parameter_names,
                output, redchisq, residual_data,
                residual_mean, residual_std, fdata, fmodel
                ):

    chains = weighted_posteriors.get_equal_weighted_posterior()
    
    figure = corner.corner(chains[:,:nparams], labels=parameter_names, 
        show_titles=True, use_math_text=True, max_n_ticks=3,
        title_fmt='.2f', quantiles=[0.16, 0.5, 0.84], color='k', 
        truth_color='bo')
    figure.savefig(output +'-TRIANGLE.png')

    pdata = qdata + 1j * udata   
    pmodel = qmodel + 1j * umodel
      
    angle_data = numpy.arctan2(udata, qdata)
    angle_model = numpy.arctan2(umodel, qmodel)
    
    ind = numpy.where(abs(fdata) == max(abs(fdata)))[0]
    rm = phi_range[ind]
    angle_data = remove_ambiguity(angle_data, rm)
    angle_model = remove_ambiguity(angle_model, rm)
    res_std1 = residual_mean - residual_std
    res_std2 = residual_mean + residual_std
 
    fig, ((aw, ax), (ay, az)) = pylab.subplots(2, 2, figsize=(15, 15))

    aw.errorbar(x, numpy.absolute(pdata), yerr=sigmap, fmt='bo',
         ms=5, ecolor='yellow', mfc='white', alpha=1, label='data')
    aw.plot(x, numpy.absolute(pmodel), 'r--', lw=2, label='model')
    aw.set_ylabel('Fractional Polarization')
    aw.set_xlabel('$\lambda^2$ [m$^{2}$]')
    aw.legend(loc='best')
            
    ax.errorbar(x, angle_data, yerr=sigmaPA, fmt='bo',
        ecolor='yellow', ms=2, mfc='white', alpha=1, label='$\chi^2_v$=%.2f'%(redchisq))
    ax.plot(x, angle_model, 'r*', ms=2, mfc='white')
    ax.set_ylabel('Polarization Angle [radians]')
    ax.set_xlabel('$\lambda^2$ [m$^{2}$]')
    ax.legend(loc='best')

    ay.plot(phi_range, numpy.absolute(fdata), 'b', lw=1)
    ay.plot(phi_range, numpy.absolute(fmodel), 'r', lw=1)
    ay.set_ylabel('Faraday Spectrum')
    ay.set_xlabel('Faraday Depth [rad m$^{-2}$]')
    #ay.legend(loc='best')

    az.plot(x, qdata, 'mo', label='q', ms=5, mfc='white')
    az.plot(x, qmodel, 'm--', lw=2, label='q model')
    az.plot(x, udata, 'go', label='u', ms=5, mfc='white')
    az.plot(x, umodel, 'g--',  lw=2, label='u model')
    az.set_ylabel('Fractional Polarization')
    az.set_xlabel('$\lambda^2$ [m$^{2}$]')
    az.legend(loc='best')

    pylab.tight_layout()
    pylab.savefig(output+ '-MODEL.png')
    #pylab.show()


