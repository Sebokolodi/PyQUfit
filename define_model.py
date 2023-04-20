import numpy

m1 = False
m11 = True
m2 = False
m4 = False


if m1:
    def model(x, cube): # m1
        """
        p = p1 * exp(2j * PA_1) exp(2j RM_1 lam2)
        """
        P = cube[0] * numpy.exp(2j * cube[1] + 2j * cube[2] * x) 
   
        return P.real, P.imag


if m11:
    def model(x, cube): # m11
        """
         p = p1 * exp(2j * PA_1) exp(2j RM_1 lam2) + 
             p2 * exp(2j * PA_2) exp(2j RM_2 lam2)
        """
        P = cube[0] * numpy.exp(2j * cube[1] + 2j * cube[2] * x) + \
            cube[3] * numpy.exp(2j * cube[4] + 2j * cube[5] * x)
   
        return P.real, P.imag
    
if m2:
    def model(x, cube): #m2
        """   
        p = p0 * exp(2j * param) exp(2j RMlam2 - 2 sigma^2 lam2^2)
        """

        P = cube[0] * numpy.exp(2j * cube[1]) * numpy.exp(2j * cube[2] * x -\
            2 * cube[3]**2 *x**2)

        return P.real, P.imag
        
if m4:
    def model(x, cube): # model 2 # m4

        """
        p = p01 * exp(2j * param1) exp(2j RM1*lam2 - 2 sigma_1^2 lam2^2)+
        p = p02 * exp(2j * param2) exp(2j RM2*lam2 - 2 sigma_2^2 lam2^2)
        """

        P = cube[0] * numpy.exp(2j * cube[1]) * numpy.exp(2j * cube[2] * x -\
           2*cube[3]**2 *x**2) + cube[4] * numpy.exp(2j * cube[5]) * \
           numpy.exp(2j * cube[6] * x - 2*cube[7]**2 *x**2)

        return P.real, P.imag
        
        

    


