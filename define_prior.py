import numpy
from priors_distribution import Priors

pri = Priors()


m1 = False
m11 = True
m2 = False
m4 = False

if m1:
    def prior(cube, ndim, nparams):

        """ Define your priors here. You will need to ensure that this is consistent with 
        define model script. For example, cube[0] for model 4 corresponds to fractional 
        polarisation. The default prior distribution is uniform ('U'). You can change the 
        intervals of the parameters (mininum and maximum). 
         """

        cube[0] = pri.GeneralPrior(cube[0], 'U', 0, 1)                  
        cube[1] = pri.GeneralPrior(cube[1], 'U', -numpy.pi/2.0, numpy.pi/2.0)
        cube[2] = pri.GeneralPrior(cube[2], 'U', -2000, 2000.0)  

if m11:
    def prior(cube, ndim, nparams):

        """ Define your priors here. You will need to ensure that this is consistent with 
        define model script. For example, cube[0] for model 4 corresponds to fractional 
        polarisation. The default prior distribution is uniform ('U'). You can change the 
        intervals of the parameters (mininum and maximum). 
         """

        cube[0] = pri.GeneralPrior(cube[0], 'U', 0, 1)                  
        cube[1] = pri.GeneralPrior(cube[1], 'U', -numpy.pi/2.0, numpy.pi/2.0)
        cube[2] = pri.GeneralPrior(cube[2], 'U', -2000, 2000.0)                 
 
        cube[3] = pri.GeneralPrior(cube[3], 'U', 0, 1)                  
        cube[4] = pri.GeneralPrior(cube[4], 'U', -numpy.pi/2.0, numpy.pi/2.0)
        cube[5] = pri.GeneralPrior(cube[5], 'U', -2000, 2000.0)   
         
if m2:
    def prior(cube, ndim, nparams):

        """ Define your priors here. You will need to ensure that this is consistent with 
        define model script. For example, cube[0] for model 4 corresponds to fractional 
        polarisation. The default prior distribution is uniform ('U'). You can change the 
        intervals of the parameters (mininum and maximum). 
         """

        cube[0] = pri.GeneralPrior(cube[0], 'U', 0, 1)                  
        cube[1] = pri.GeneralPrior(cube[1], 'U', -numpy.pi/2.0, numpy.pi/2.0)
        cube[2] = pri.GeneralPrior(cube[2], 'U', -2000, 2000.0)
        cube[3] = pri.GeneralPrior(cube[3], 'U', 0, 400)                  

    
if m4:
    def prior(cube, ndim, nparams):

        """ Define your priors here. You will need to ensure that this is consistent with 
        define model script. For example, cube[0] for model 4 corresponds to fractional 
        polarisation. The default prior distribution is uniform ('U'). You can change the 
        intervals of the parameters (mininum and maximum). 
         """

        cube[0] = pri.GeneralPrior(cube[0], 'U', 0, 1)                  
        cube[1] = pri.GeneralPrior(cube[1], 'U', -numpy.pi/2.0, numpy.pi/2.0)
        cube[2] = pri.GeneralPrior(cube[2], 'U', -2000, 2000.0)
        cube[3] = pri.GeneralPrior(cube[3], 'U', 0, 400)                  
 
  
        cube[4] = pri.GeneralPrior(cube[4], 'U', 0, 1)                  
        cube[5] = pri.GeneralPrior(cube[5], 'U', -numpy.pi/2.0, numpy.pi/2.0)
        cube[6] = pri.GeneralPrior(cube[6], 'U', -2000, 2000.0)
        cube[7] = pri.GeneralPrior(cube[7], 'U', 0, 400)  


