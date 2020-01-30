#! /usr/bin/env python

import pymultinest
import os, sys
import model
import prior
import argparse
import numpy
import corner
import pylab
import json
import time


def LogLike(cube, ndim, nparams):

    qmod, umod  = model.model(x, cube)
    chisq = 0.5 * ((q - qmod)/sigmaq)**2.0
    chisu = 0.5 * ((u - umod)/sigmau)**2.0
    prefactor = numpy.log(numpy.pi * sigmaq * sigmau)
    
    return numpy.sum(prefactor + chisq + chisu) * -1



if __name__=='__main__':
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

    Qdata = [item for item in args.qdata[0].split(',')]
    Udata = [item for item in args.udata[0].split(',')]
    Idata = [item for item in args.idata[0].split(',')] 
    noiseQ = [item for item in args.noiseq[0].split(',')]
    noiseU = [item for item in args.noiseu[0].split(',')]
    noiseI = [item for item in args.noisei[0].split(',')] 
    frequency = [item for item in args.freq[0].split(',')]

    nparams = args.num_parameters
    outdir = args.outdir or 'FIT'
    if not os.path.exists(outdir): os.mkdir(outdir)

    ndimQ, ndimU, ndimI = len(Qdata), len(Udata), len(Idata)
    ndimNoiseQ, ndimNoiseU, ndimNoiseI = len(noiseQ), len(noiseU), len(noiseI)

    if ndimNoiseQ != ndimNoiseU or ndimNoiseQ != ndimNoiseI:
        sys.exit('The number of noise files are not equal.')
    
    if ndimQ != ndimU or ndimQ != ndimI:
        sys.exit('The number of Q/U/I files are not consistent.')

    if ndimNoiseQ != ndimQ:
        if ndimNoiseQ == 1:
            noiseQ = noiseQ * ndimQ
            noiseU = noiseU * ndimQ
            noiseI = noiseI * ndimQ
    if len(frequency) != ndimQ:
        if len(frequency) == 1:
            frequency = frequency * ndimQ
        
    write_stats = open(outdir + '/LoS-Stats.txt', 'w')

    for i, (ftxt, qtxt, utxt, itxt, nq, nu, ni) in \
           enumerate(zip(frequency, Qdata, Udata, Idata, noiseQ, 
               noiseU, noiseI)):

        Q =  numpy.loadtxt(qtxt)
        U =  numpy.loadtxt(utxt)
        I =  numpy.loadtxt(itxt)
        noiseq =  numpy.loadtxt(nq)
        noiseu =  numpy.loadtxt(nu)
        noisei =  numpy.loadtxt(ni)
        freq =  numpy.loadtxt(ftxt)

    
        prefix = args.prefix or 'FIT-QU' 
        prefix =  prefix + '-%d'%i
        output =  os.path.join(outdir, prefix)

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


        parameter_names =  ['param%d'%i for i in range(nparams)]
    
        start = time.time()
        pymultinest.run(LogLike, prior.prior, nparams,
                   outputfiles_basename=output, **kwargs)
        a = pymultinest.analyse.Analyzer(nparams, 
                   outputfiles_basename=output)
        end = time.time()

        Dtime = end-start
       
        best_fits = numpy.array([float(param) for param
                in a.get_best_fit()['parameters']])

        logZ = a.get_best_fit()['log_likelihood']
 
        array = numpy.hstack((best_fits, logZ, Dtime))

        for i in range(nparams + 2):
            write_stats.write('%.4f \t'%array[i])
        write_stats.write('\n')
          
        if args.make_plots:
       
            chains = a.get_equal_weighted_posterior()
            figure = corner.corner(chains[:,:nparams], labels=parameter_names, 
                   show_titles=True, use_math_text=True, max_n_ticks=3, 
                   title_fmt='.3f')
            figure.savefig(output +'-TRIANGLE.PNG')

            qmodel, umodel = model.model(x, best_fits)
            fig, (ax, ay) = pylab.subplots(1, 2, figsize=(10, 5))
            ax.errorbar(x, q, yerr=sigmaq, fmt='bo', alpha=0.4, ecolor='y', ms=2) 
            ax.plot(x, qmodel, 'r', lw=2)
            ax.set_ylabel('fractional q')
            ax.set_xlabel('$\lambda^2$ [m$^{-2}$]')
            
            ay.errorbar(x, u, yerr=sigmau, fmt='bo', alpha=0.4, ecolor='y', ms=2) 
            ay.plot(x, umodel, 'r', lw=2)
            ay.set_ylabel('fractional u')
            ay.set_xlabel('$\lambda^2$ [m$^{-2}$]')
            pylab.tight_layout()
            pylab.savefig(output+ '.PNG')
            #pylab.show()

    write_stats.close()

     


        
            
   
         

             

    

        













