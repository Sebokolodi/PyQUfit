#! /usr/bin/env python

import pymultinest
import os, sys
import define_model
import define_prior
import argparse
import numpy
import json
import time
import doplots 
import rmsynthesis as rmsyn
import glob


#TODO: try also fitting in Faraday depth space. 
#TODO: run from fits images. User must provide pixel locations.



def LogLike(cube, ndim, nparams):

    qmod, umod  = define_model.model(x, cube)
    chisq = 0.5 * ((q - qmod)/sigmaq)**2.0
    chisu = 0.5 * ((u - umod)/sigmau)**2.0
    prefactor = numpy.log(numpy.pi * sigmaq * sigmau)
    
    return numpy.sum(prefactor + chisq + chisu) * -1


def check_number_of_files(freq_files, data_files, noise_files):

    """ checks if the number of files provided consistent.
    If files for the noise/frequency is not equal to files 
    of the data, then use first noise/frequency file all the data,
    and ignore the rest.

    
    """

    freq  = [item for item in freq_files[0].split(',')]
    input_data = [item for item in data_files[0].split(',')]
    input_noise = [item for item in noise_files[0].split(',')]


    if  len(input_noise) != len(input_data):
        input_noise = [input_noise[0]] * len(input_data)
    if len(freq) != len(input_data):
        freq = [freq[0]] * len(input_data)

    return freq, input_data, input_noise



def return_sigma(q, u, pabs, Q, U, I, Pabs, nQ, nU, nI):

    """compute the sigma from the noise """

    sq = abs(q * ( (nQ/Q)**2 + (nI/I)**2 )**0.5)
    su = abs(u * ( (nU/U)**2 + (nI/I)**2 )**0.5)
    np = ( (Q/Pabs)**2 * nQ**2 + (U/Pabs)**2 * nU**2 )**0.5
    sp = pabs * ((np/Pabs)**2 + (nI/I)**2 )**0.5
    sPA = (((u * sq)**2 + (q * su)**2)/(4.0* pabs**4 ))**0.5

    return sq, su, sp, sPA


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


def extrap(x, xp, yp):
    """
    Wrapper to allow numy.interp to linearly extrapolate at function ends.
    
    np.interp function with linear extrapolation
    http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate
    -give-a-an-extrapolated-result-beyond-the-input-ran
    """
    y = numpy.interp(x, xp, yp)
    y = numpy.where(x < xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = numpy.where(x > xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]),
                 y)
    return y
     
     

if __name__=='__main__':

#def main_text():
    parser = argparse.ArgumentParser(description='Performs QU-fitting ' 
             'to the data.')
    add = parser.add_argument

    add('-json', '--json', dest='json_file', help='a json file with multinest settings.')
    add('-indata', '--input-data', dest='data', action='append',
       help='Text file should contain Stokes Q, U, and I with shape (N, 3) where N is data length.')  
    add('-mod', dest='model', help='A model to call', type=str, default='m1')
    add('-rm', '--error-clip', dest='error_clip', help='Error to clip in fractional q and u.' 
        ' default is 10 percent', type=float, default=0.8)
    add('-nparam', '--numparam', dest='num_parameters', help='Number of parameters.', type=int) 
    add('-phi', '--faraday-depth', dest='faraday_depth', help='Maximum Faraday depth in rad/m/m',
       type=float, default=5000) 
    add('-dphi', '--sample-depth', dest='sample_depth', help='Faraday depth sampling in rad/m/m',
        type=float, default=10) 
    add('-niter', '--rmclean-niter', dest='rmsyn_niter', 
        help='RMclean number of iteratation. Useful for plotting. Default=1000',
        type=float, default=10000) 
    add('-gain', '--rmclean-gain', dest='rmsyn_gain', 
        help='RMclean gain factor. Useful when plotting. Default=0.1',
        type=float, default=0.1) 
    add('-plot', '--make-plot', dest='make_plots', help='Make plots', type=bool, default=False) 
    add('-pref', '--prefix', dest='prefix', help='Prefix for output files.')
    add('-outdir', '--outdir', dest='outdir', help='Output directory for output files. '
        'if not there already it will be created. If not specified, all plots will be '
        'saved in "FIT". ')

    args = parser.parse_args()

    model_type = args.model
    
    
    
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

    #freq_files, data_files, noise_files = check_number_of_files(
    #      args.freq, args.data, args.noise)

    input_data = args.data
    
    residuals = lambda data, model, sigmas: (data-model)/sigmas



    #for i, (freqtxt, datatxt, noisetxt) in enumerate(zip(freq_files,
    #       data_files, noise_files)):
        #i += 201
    
    
    prefix = args.prefix or 'FIT-QU'
    #prefix =  prefix + '-%d-'%i
    output =  os.path.join(outdir, prefix)
    print(input_data[0])
    #try:
    input_data_LoS =  numpy.loadtxt(input_data[0])
    #except:
    #     sys.exit('>>> Something wrong with your input data files.')

        #try:
        #    input_noise_LoS =  numpy.loadtxt(noisetxt)
        #except:
        #    sys.exit('>>> Cannot open the input noise files.')

        #try:
        #    input_freq_LoS =  numpy.loadtxt(freqtxt)
        #except:
        #    sys.exit('>>> something wrong with your input frequency file.')
        
    input_freq_LoS = input_data_LoS[:, 0]
    Q = input_data_LoS[:, 2]
    U = input_data_LoS[:, 3]
    I = input_data_LoS[:, 1]
    noiseq = input_data_LoS[:, 5]
    noiseu = input_data_LoS[:, 6]
    noisei = input_data_LoS[:, 4]
    x = (3.0e8/input_freq_LoS)**2

    if len(Q) != len(noiseq) or len(Q) != len(x):
         sys.exit('>>> Length of input data, noise, and frequency not equal.')        

    q, u = Q/I, U/I
    pabs = pow( pow(q, 2) + pow(u, 2), 0.5)
    Pabs = pow( pow(Q, 2) + pow(U, 2), 0.5)


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

    residual_data = 0.5 * (residualq + residualu) 

    errors_best_fits = numpy.array([a.get_stats()['marginals'][k]['sigma']
                for k in range(nparams)])
    loglike = a.get_best_fit()['log_likelihood']
    output_stats = a.get_mode_stats()
    logZ = output_stats['global evidence']
    logZerr = output_stats['global evidence error']
    AIC = 2 * nparams - 2 * loglike
    BIC = nparams * numpy.log(Ndata) - 2 * loglike
    res_mean = residual_data.mean()
    res_std = residual_data.std()
    
    pdata = q + 1j * u   
    pmodel = qmodel + 1j * umodel
    
    phi_range = numpy.arange(-args.faraday_depth, 
                args.faraday_depth, args.sample_depth)
                    
    
    FWHM = 3.8/(x[-1]-x[0])
    sigma = FWHM/(2.35482)
          
          
    phi_range, fcleaned_data, fcomp = rmsyn.main(lam=x, pdata=pdata, phi_max=args.faraday_depth, phi_step=args.sample_depth,
         niter=args.rmsyn_niter, gain=args.rmsyn_gain, plot=False)
    
    phi_range, fcleaned_model, fcomp = rmsyn.main(lam=x, pdata=pmodel, phi_max=args.faraday_depth, phi_step=args.sample_depth,
         niter=args.rmsyn_niter, gain=args.rmsyn_gain, plot=False)
              
   
    store_output = numpy.hstack((best_fits, errors_best_fits, loglike, 
              logZ, logZerr, redchisq, res_mean, res_std, AIC, BIC, Dtime))

    # write the output
    for k in range(len(store_output)):
        write_stats.write('%.4f \t'%store_output[k])
    write_stats.write('\n')

    # making plots.
    if args.make_plots:
        doplots.do_plots(x=x, qdata=q, udata=u, sigmap=sigmap, sigmaPA=sigmaPA,
                qmodel=qmodel, umodel=umodel, nparams=nparams,  
                weighted_posteriors=a, parameter_names=parameter_names,
                output=output, redchisq=redchisq, residual_data=residual_data, 
                residual_mean=res_mean, fdata=fcleaned_data, fmodel=fcleaned_model,
                residual_std=res_std, phi_range=phi_range)
    write_stats.close()
    
    remove_file = ['.points', '.dat', 'summary.txt', '.ptprob', '.iterinfo']
    for remove in remove_file: 
        files = glob.glob(output + '*%s'%remove)
        for f in files:
            os.system('rm -rf %s'%f)

    print('Execution time %.2f seconds'%Dtime)
