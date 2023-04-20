# PyQUfit
A software program designed for fitting a user-defined model to Stokes Q and U astronomical data. This employs nested sampling algorithm via PyMultinest [reference here] to evaluate the evidence and parameters from the posterior probabilities. 

Installation:

At this point, you simply clone this github page. Not yet on pypi.


Package requirements:

1. pyMultinest   
2. Numpy   
3. Matplotlib  


The files:

1. Required a text file with freq, Stokes Q, U, I, Error in Stokes Q, U and I.     

2. Required information  
   For now, the user need to specify using True/False, the model they want to fit. I plan to automate this at some point.   
   A user needs to specify this inside the following files:
   i.'define_model.py'        
   ii.'define_priors.py'
   Script (i) defines the models to use. The user can easily update this with their model of interest. 
   Script (ii) defines the priors intervals to use for each model. If a user changes the model in script (i), they will
   also need to change it in script (ii).


How to find options to parse:  
 
    ./qufit.py -h 


So on a terminal you can simply do:

    ./qufit.py -indata data.txt -nparm 4 -json config.json -phi_range 5000 -dphi 100 -niter 1000 -gain 0.1 -plot True -pref test -outdir TEST 
  
Where -indata is the input data (text file) containting freq, Stokes Q, U, and I, Error in Q, U and I. With shape (N, 7) where N is data length.
-json is Jason file, this is needed -- it configures the multinest inputs, -phi_range is the maximum Faraday depth, and dphi is the sampling size of the Faraday depths in rad/m/m. Niter is number of iterations to use for iteration when performing RM-synthesis clean, and gain is similar to any cleaning (RM-synthesis is  performed mainly for plotting purposes so it is optional). If you want plots output, indicate with True. The code plots Fractional pol vs wavelength squared, position angle vs. wavelength-squared, and Faraday Spectrum. 


NOTE: at this point, this is done mainly for text files. The reason being is that QU-fitting through multinest takes long especially for the data sizes that led to this code (~ 600000 pixels each with 2000  planes). Usually multinest takes roughly 20 seconds for a single pixel, with 4 parameter model and 1200 frequency planes, thus, for 600000 pixels this amount to ~3000 hours (125 days)!!! 



You can also provide a list of lines of sight on the terminal, separated by a comma. However, in case this is tedious, we encourage using run-all.py, in this you read all the data files into the script and pass it to qufit.py. All you need in this case would be 

    ipython run-all.py


Output:

1. Multinest output   
2. Plot showing fractional vs. lam^2, PA vs lam^2, Faraday spectrum vs Faraday depth, and residuals.    
3. Text file with fitted parameters, errors, loglikelihood, evidence, error in evidence, reduced Chi-squared, AIC, BIC, and excution time.    
