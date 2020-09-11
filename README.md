# PyQUfit
A software program designed for fitting a user-defined model to Stokes Q and U astronomical data. This employs nested sampling algorithm via PyMultinest [reference here] to evaluate the evidence and parameters from the posterior probabilities. 

Installation:


Package requirements:
pyMultinest 
Numpy
Matplotlib


Executing the code:

1. Required Files (Not optional):   
  1.1 Stokes Q, U and I data in separate text files.      
  1.2 Frequency text file.      
  1.3 Noise in Q, U, and I images.      

2. Required information     
  2.1 user-defined model in a script 'define_model.py'        
  2.2 user-defined priors intervals in a script 'define_priors.py'      


How to run this:

./qufit.py -h 

Performs QU-fitting to the data.

optional arguments:
  -h, --help            show this help message and exit
  -json JSON_FILE, --json JSON_FILE
                        a json file with multinest settings.
  -f FREQ, --freq FREQ  Frequency file (text)
  -indata DATA, --input-data DATA
                        Text file should contain Stokes Q, U, and I with shape
                        (N, 3) where N is data length.
  -noise NOISE, --noise NOISE
                        Text file should contain noise in Stokes Q, U and I
                        maps. See indata.
  -rm ERROR_CLIP, --error-clip ERROR_CLIP
                        Error to clip in fractional q and u. default is 10
                        percent
  -nparam NUM_PARAMETERS, --numparam NUM_PARAMETERS
                        Number of parameters.
  -phi FARADAY_DEPTH, --faraday-depth FARADAY_DEPTH
                        Maximum Faraday depth in rad/m/m
  -dphi SAMPLE_DEPTH, --sample-depth SAMPLE_DEPTH
                        Faraday depth sampling in rad/m/m
  -niter RMSYN_NITER, --rmclean-niter RMSYN_NITER
                        RMclean number of iteratation. Useful for plotting.
                        Default=1000
  -gain RMSYN_GAIN, --rmclean-gain RMSYN_GAIN
                        RMclean gain factor. Useful when plotting. Default=0.1
  -plot MAKE_PLOTS, --make-plot MAKE_PLOTS
                        Make plots
  -pref PREFIX, --prefix PREFIX
                        Prefix for output files.
  -outdir OUTDIR, --outdir OUTDIR
                        Output directory for output files. if not there
                        already it will be created. If not specified, all
                        plots will be saved in "FIT".



So on a terminal you can simply do:

./qufit.py -indata data.txt -f frequency.txt  -noise noisedata.txt -nparm 4 -json config.json -phi_range 5000 -dphi 100 -niter 1000 -gain 0.1 -plot True -pref test -outdir TEST 

You can also provide a list of lines of sight on the terminal, separated by a comma. However, in case this is tedious, we encourage using run-all.py, in this you read all the data files into the script and pass it to qufit.py. All you need in this case would be 

ipython run-all.py
