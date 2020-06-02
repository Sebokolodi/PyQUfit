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


def check_data():

    return 


def Faraday2Lambda():

    return


def Lambda2Faraday():

    return



def do_plots(x, qdata, udata, sigmaq, sigmau, nparams, 
                weighted_posteriors, best_fits, parameter_names,
                output):

    chains = weighted_posteriors.get_equal_weighted_posterior()
    
    figure = corner.corner(chains[:,:nparams], labels=parameter_names, 
        show_titles=True, use_math_text=True, max_n_ticks=3, title_fmt='.3f')
    figure.savefig(output +'TRIANGLE.png')

    qmodel, umodel = define_model.model(x, best_fits)
           
    fig, (ax, ay) = pylab.subplots(1, 2, figsize=(10, 5))
    ax.errorbar(x, qdata, yerr=sigmaq, fmt='co', alpha=0.4, ecolor='c', ms=2) 
    ax.plot(x, qmodel, 'r', lw=2)
    ax.set_ylabel('fractional q')
    ax.set_xlabel('$\lambda^2$ [m$^{-2}$]')
            
    ay.errorbar(x, udata, yerr=sigmau, fmt='co', alpha=0.4, ecolor='c', ms=2) 
    ay.plot(x, umodel, 'r', lw=2)
    ay.set_ylabel('fractional u')
    ay.set_xlabel('$\lambda^2$ [m$^{-2}$]')
    pylab.tight_layout()
    pylab.savefig(output+ 'model.png')



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

    for i, (ftxt, qtxt, utxt, itxt, nq, nu, ni) in \
           enumerate(zip(freq, Qdata, Udata, Idata, Qnoise, 
               Unoise, Inoise)):

        #i += 2001
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
        p = q + 1j * u
        x = (3.0e8/freq)**2

        sigmaq = abs(q * ( (noiseq/Q)**2 + (noisei/I)**2 )**0.5)
        sigmau = abs(u * ( (noiseu/U)**2 + (noisei/I)**2 )**0.5)

        # remove Nan from the data.
        ind = ~numpy.isnan(numpy.absolute(p)*sigmaq*sigmau)
        x = x[ind]
        q = q[ind]
        u = u[ind]
        sigmaq = sigmaq[ind]
        sigmau = sigmau[ind]

        # remove channels with error in fractional pol > error_clip
        ind_clip = numpy.where( (sigmaq + sigmau)/2.0 > clip_error)
        x = numpy.delete(x, ind_clip)
        q = numpy.delete(q, ind_clip)
        u = numpy.delete(u, ind_clip)
        sigmaq = numpy.delete(sigmaq, ind_clip)
        sigmau = numpy.delete(sigmau, ind_clip)

        parameter_names =  ['param%d'%k for k in range(nparams)]
        # run multinest    
        start = time.time()
        pymultinest.run(LogLike, define_prior.prior, nparams,
                   outputfiles_basename=output, **kwargs)
        a = pymultinest.analyse.Analyzer(nparams, 
                   outputfiles_basename=output)
        end = time.time()
        Dtime = end-start

        # get important fitting results.
        
        best_fits = numpy.array([float(param) for param
                in a.get_best_fit()['parameters']])

        errors_best_fits = numpy.array([a.get_stats()['marginals'][k]['sigma']
                    for k in range(nparams)])
        loglike = a.get_best_fit()['log_likelihood']
        stats = a.get_mode_stats()
        logZ = stats['global evidence']
        logZerr = stats['global evidence error']
        AIC = 2 * nparams - 2 * loglike
        BIC = nparams * numpy.log(len(x)) - 2 * loglike

        store_output = numpy.hstack((best_fits, errors_best_fits, loglike, 
              logZ, logZerr, AIC, BIC, Dtime))

        # write the output
        for k in range(nparams * 2 + 6):
            write_stats.write('%.4f \t'%store_output[k])
        write_stats.write('\n')

        # making plots.
        if args.make_plots:
            do_plots(x, q, u, sigmaq, sigmau, nparams, 
                a, best_fits, parameter_names,
                output)

    write_stats.close()

#main_text()

