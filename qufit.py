#! /usr/bin/env python

import pymultinest
import os, sys
import define_model
import define_prior
import argparse
import numpy
import corner
import pylab
import json
import time
import matplotlib
from scipy import stats
matplotlib.rcParams.update({'font.size':16})



#TODO: include RM synthesis?
#TODO: try also fitting in Faraday depth space. 
#TODO: run from fits images. User must provide pixel locations.


def LogLike(cube, ndim, nparams):

    qmod, umod  = define_model.model(x, cube)
    chisq = 0.5 * ((q - qmod)/sigmaq)**2.0
    chisu = 0.5 * ((u - umod)/sigmau)**2.0
    prefactor = numpy.log(numpy.pi * sigmaq * sigmau)
    
    return numpy.sum(prefactor + chisq + chisu) * -1


def return_sigma(q, u, pabs, Q, U, I, Pabs, nQ, nU, nI):

    """compute the sigma from the noise """

    sq = abs(q * ( (nQ/Q)**2 + (nI/I)**2 )**0.5)
    su = abs(u * ( (nU/U)**2 + (nI/I)**2 )**0.5)
    np = ( (Q/Pabs)**2 * nQ**2 + (U/Pabs)**2 * nU**2 )**0.5
    sp = pabs * ((np/Pabs)**2 + (nI/I)**2 )**0.5
    sPA = (((u * sq)**2 + (q * su)**2)/(4.0* pabs**4 ))**0.5

    return sq, su, sp, sPA


def Faraday2Lambda():

    return


def check_length_data(f, Q, U, I, sigQ, sigU, sigI):

    """ checks if the number of files provided consistent """

    freq  = [item for item in f[0].split(',')]
    Qdata = [item for item in Q[0].split(',')]
    Udata = [item for item in U[0].split(',')]
    Idata = [item for item in I[0].split(',')] 
    Qnoise = [item for item in sigQ[0].split(',')]
    Unoise = [item for item in sigU[0].split(',')]
    Inoise = [item for item in sigI[0].split(',')] 

    # check the lengths of data in the files are consistent.
    ndimQ, ndimU, ndimI = len(Qdata), len(Udata), len(Idata)
    ndimNoiseQ, ndimNoiseU, ndimNoiseI = len(Qnoise), len(Unoise), len(Inoise)

    if ndimNoiseQ != ndimNoiseU or ndimNoiseQ != ndimNoiseI:
        sys.exit('The number of noise files are not equal.')
    
    if ndimQ != ndimU or ndimQ != ndimI:
        sys.exit('The number of Q/U/I files are not consistent.')

    if ndimNoiseQ != ndimQ:
        if ndimNoiseQ == 1:
            Qnoise = Qnoise * ndimQ
            Unoise = Unoise * ndimQ
            Inoise = Inoise * ndimQ
    if len(freq) != ndimQ:
        if len(freq) == 1:
            freq = freq * ndimQ

    return freq, Qdata, Udata, Idata, Qnoise, Unoise, Inoise



def do_plots(x, qdata, udata, sigmap, sigmaPA, qmodel, umodel, nparams, 
                weighted_posteriors, parameter_names,
                output, redchisq, residual_data, pval):

    chains = weighted_posteriors.get_equal_weighted_posterior()
    
    figure = corner.corner(chains[:,:nparams], labels=parameter_names, 
        show_titles=True, use_math_text=True, max_n_ticks=3,
        title_fmt='.2f', quantiles=[0.16, 0.5, 0.84], color='k', 
        truth_color='bo')
    figure.savefig(output +'TRIANGLE.png')

    
    pdata = (qdata**2 + udata**2)**0.5   
    pmodel = (qmodel**2 + umodel**2)**0.5    
    angle_data = numpy.arctan2(udata, qdata)
    angle_model = numpy.arctan2(umodel, qmodel)

    distrib = numpy.arange(-30, 30, 0.1)
    normal = (1.0/(2*numpy.pi)**0.5) * numpy.exp(-0.5 * distrib**2)
 
    fig, (ax, ay, az) = pylab.subplots(1, 3, figsize=(15, 5))
    ax.set_title('$\chi^2_v$=%.2f, p-val=%.3f'%(redchisq, pval))
    ax.errorbar(x, pdata, yerr=sigmap, fmt='bo', ms=2, ecolor='yellow')
    ax.plot(x, pmodel, 'r', lw=2)
    ax.set_ylabel('Fractional Polarization')
    ax.set_xlabel('$\lambda^2$ [m$^{-2}$]')
            
    ay.errorbar(x, angle_data, yerr=sigmaPA, fmt='b.', ecolor='yellow', ms=2)
    ay.plot(x, angle_model, 'r.', ms=2)
    ay.set_ylabel('Polarization Angle [radians]')
    ay.set_xlabel('$\lambda^2$ [m$^{-2}$]')

    az.hist(residual_data, histtype='step', density=True, color='b', lw=2)
    az.plot(distrib, normal, color='black', lw=2)
    az.set_xlabel('Normalized Residuals')
    az.set_ylabel('Probability Distribution')

    pylab.tight_layout()
    pylab.savefig(output+ 'model.png')
    pylab.show()



if __name__=='__main__':

#def main_text():
    parser = argparse.ArgumentParser(description='Performs QU-fitting ' 
             'to the data.')
    add = parser.add_argument

    add('-json', '--json', dest='json_file', help='a json file with multinest settings.')
    add('-f', '--freq', dest='freq', action='append', help='Frequency file (text)')
    add('-q', '--qdata', dest='qdata', action='append', help='Stokes q')
    add('-u', '--udata', dest='udata', action='append', help='Stokes u')
    add('-i', '--idata', dest='idata', action='append', help='Stokes i')  
    add('-nq', '--noiseq', dest='noiseq', action='append', help='Noise in q') 
    add('-nu', '--noiseu', dest='noiseu', action='append', help='Noise in u') 
    add('-ni', '--noisei', dest='noisei', action='append', help='Noise in i')
    add('-ec', '--eclip', dest='error_clip', help='Error to clip in fractional q and u.' 
        ' default is 10 percent', type=float, default=0.1)
    add('-nparam', '--nparam', dest='num_parameters', help='Number of parameters.', type=int) 
    add('-plot', '--make-plot', dest='make_plots', help='Make plots', type=bool, default=False) 
    add('-pref', '--prefix', dest='prefix', help='Prefix for output files.')
    add('-out', '--outdir', dest='outdir', help='Output directory for output files. '
        'if not there already it will be created. If not specified, all plots will be '
        'saved in "FIT". ')

    args = parser.parse_args()

    with open(args.json_file, 'r') as settings:
        settings = json.load(settings)
    kwargs = settings['multinest']        

    outdir = args.outdir or 'FIT'
    if not os.path.exists(outdir): os.mkdir(outdir)

    prefix = args.prefix or 'FIT-QU' 
    output =  os.path.join(outdir, prefix)
    write_stats = open(output + '-fitparameters.txt', 'w')

    nparams = args.num_parameters
    clip_error = args.error_clip

    freq, Qdata, Udata, Idata, Qnoise, Unoise, Inoise = check_length_data(
          args.freq, args.qdata, args.udata, args.idata, args.noiseq,
          args.noiseu, args.noisei)

    residuals = lambda data, model, sigmas: (data-model)/sigmas


    for i, (ftxt, qtxt, utxt, itxt, nq, nu, ni) in \
           enumerate(zip(freq, Qdata, Udata, Idata, Qnoise, 
               Unoise, Inoise)):

        prefix = args.prefix or 'FIT-QU'
        prefix =  prefix + '-%d-'%i
        output =  os.path.join(outdir, prefix)

        Q =  numpy.loadtxt(qtxt)
        U =  numpy.loadtxt(utxt)
        I =  numpy.loadtxt(itxt)
        noiseq =  numpy.loadtxt(nq)
        noiseu =  numpy.loadtxt(nu)
        noisei =  numpy.loadtxt(ni)
        freq =  numpy.loadtxt(ftxt)

        q, u = Q/I, U/I
        pabs = pow( pow(q, 2) + pow(u, 2), 0.5)
        Pabs = pow( pow(Q, 2) + pow(U, 2), 0.5)
        x = (3.0e8/freq)**2

        # compute sigma
        sigmaq, sigmau, sigmap, sigmaPA = return_sigma(
               q, u, pabs, Q, U, I, Pabs, noiseq,
               noiseu, noisei) 

        # remove Nan from the data.
        ind = ~numpy.isnan(pabs*sigmaq*sigmau)
        x = x[ind]
        q = q[ind]
        u = u[ind]
        sigmaq = sigmaq[ind]
        sigmau = sigmau[ind]
        sigmap = sigmap[ind]
        sigmaPA = sigmaPA[ind]

        # remove channels with error in fractional pol > error_clip
        ind_clip = numpy.where( (sigmaq + sigmau)/2.0 > clip_error)
        x = numpy.delete(x, ind_clip)
        q = numpy.delete(q, ind_clip)
        u = numpy.delete(u, ind_clip)
        sigmap = numpy.delete(sigmap, ind_clip)
        sigmaPA = numpy.delete(sigmaPA, ind_clip)
        sigmaq = numpy.delete(sigmaq, ind_clip)
        sigmau = numpy.delete(sigmau, ind_clip)


        # define parameters
        parameter_names =  ['p$_{%d}$'%k for k in range(nparams)]

        # run multinest    
        start = time.time()
        pymultinest.run(LogLike, define_prior.prior, nparams,
                   outputfiles_basename=output, **kwargs)
        a = pymultinest.analyse.Analyzer(nparams, 
                   outputfiles_basename=output)
        end = time.time()
        Dtime = end - start

        #TODO: return reduced chi-squared, and a histogram and p-values.

        # get important fitting results.
        best_fits = numpy.array([float(param) for param
                in a.get_best_fit()['parameters']])

        qmodel, umodel = define_model.model(x, best_fits)
      
        residualq = residuals(q, qmodel, sigmaq)
        residualu = residuals(u, umodel, sigmau)
        Ndata = len(x)
        Nnorm = pow((Ndata - nparams), -1)
         
        redchisq_q = Nnorm * sum(pow(residualq, 2))
        redchisq_u = Nnorm * sum(pow(residualq, 2))
        redchisq = 0.5 * (redchisq_q + redchisq_u)

        residual_data = numpy.hstack((residualq, residualu))
        histbins, histresidual = numpy.histogram(residual_data, density=True)
        statistics, pvalue = stats.kstest(histresidual, 'norm')


        errors_best_fits = numpy.array([a.get_stats()['marginals'][k]['sigma']
                    for k in range(nparams)])
        loglike = a.get_best_fit()['log_likelihood']
        output_stats = a.get_mode_stats()
        logZ = output_stats['global evidence']
        logZerr = output_stats['global evidence error']
        AIC = 2 * nparams - 2 * loglike
        BIC = nparams * numpy.log(Ndata) - 2 * loglike

        store_output = numpy.hstack((best_fits, errors_best_fits, loglike, 
              logZ, logZerr, redchisq, pvalue, AIC, BIC, Dtime))

        # write the output
        for k in range(len(store_output)):
            write_stats.write('%.4f \t'%store_output[k])
        write_stats.write('\n')

        # making plots.
        if args.make_plots:
            do_plots(x, q, u, sigmap, sigmaPA, qmodel, 
                umodel, nparams, a, parameter_names, 
                output, redchisq, residual_data, pvalue)

    write_stats.close()

#main_text()

