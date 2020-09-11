# PyQUfit
A software program designed for fitting a user-defined model to Stokes Q and U astronomical data. This employs nested sampling algorithm via PyMultinest [reference here] to evaluate the evidence and parameters from the posterior probabilities. 

Installation:

At this point, you simply clone this github page. Not yet in pypi.



Package requirements:

1. pyMultinest   
2. Numpy   
3. Matplotlib  


The files:

1. Required Files (Not optional):   
  1.1 Stokes Q, U and I data in separate text files.      
  1.2 Frequency text file.      
  1.3 Noise in Q, U, and I images.      

2. Required information     
  2.1 user-defined model in a script 'define_model.py'        
  2.2 user-defined priors intervals in a script 'define_priors.py'      


How to find options to parse:  
 
    ./qufit.py -h 


So on a terminal you can simply do:

    ./qufit.py -indata data.txt -f frequency.txt  -noise noisedata.txt -nparm 4 -json config.json -phi_range 5000 -dphi 100 -niter 1000 -gain 0.1 -plot True -pref test -outdir TEST 
  
Where -indata is the input data (text file) containting Stokes Q, U, and I with shape (N, 3) where N is data length. -f is text file containing the frequencies, -noise is the noise estimates for each plane (should be same dimension as the indata), -json is Jason file, this is needed -- it configures the multinest inputs, -phi_range is the maximum Faraday depth, and dphi is the sampling size of the Faraday depths in rad/m/m. Niter is number of iterations to use for iteration when performing RM-synthesis clean, and gain is similar to any cleaning (RM-synthesis is  performed mainly for plotting purposes TODO: should make optional). If you want plots output, indicate with True. The code plots Fractional pol vs wavelength squared, position angle vs. wavelength-squared, and Faraday Spectrum. 


NOTE: at this point, this is done mainly for text files. The reason being is that QU-fitting through multinest takes long especially for the data sizes that led to this code (~ 600000 pixels each with 2000  planes). Usually multinest takes roughly 20 seconds for a single pixel, with 4 parameter model and 1200 frequency planes, thus, for 600000 pixels this amount to ~3000 hours (125 days)!!! 



You can also provide a list of lines of sight on the terminal, separated by a comma. However, in case this is tedious, we encourage using run-all.py, in this you read all the data files into the script and pass it to qufit.py. All you need in this case would be 

    ipython run-all.py


Output:

1. Multinest output   
2. Plot showing fractional vs. lam^2, PA vs lam^2, Faraday spectrum vs Faraday depth, and residuals.    
3. Text file with fitted parameters, errors, loglikelihood, evidence, error in evidence, reduced Chi-squared, AIC, BIC, and excution time.    
