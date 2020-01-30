import numpy

def model(x, cube):

    P = cube[0] * numpy.exp(2j * cube[1]) * ( numpy.sinc(cube[4] * x) *
          numpy.exp(2j * cube[2] * x - 2 * cube[3]**2 * x**2 ))

    return P.real, P.imag




