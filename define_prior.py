import numpy
from priors_distribution import Priors

pri = Priors()


def prior(cube, ndim, nparams):

    cube[0] = pri.GeneralPrior(cube[0], 'U', 0, 1)                  
    cube[1] = pri.GeneralPrior(cube[1], 'U', -numpy.pi/2.0, numpy.pi/2.0)
    cube[2] = pri.GeneralPrior(cube[2], 'U', -6500.0, 6500.0)
    cube[3] = pri.GeneralPrior(cube[3], 'U', 0, 2500.0) 
    cube[4] = pri.GeneralPrior(cube[4], 'U', -2500.0, 2500.0) 
 

    




