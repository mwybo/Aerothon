"""
A module containing a generic Secant solver.
"""
from __future__ import division # let 5/2 = 2.5 rather than 2

#
# Try to import psycho to improve performance
#
try:
    import psyco
    psyco.full()
except ImportError:
#    print "You do not have Psycho installed. Psycho can drastically improve performance"
    pass

import numpy as npy
from scalar.scalar import _scalar
import copy

################################################################################
class SecantError(Exception):
    """ Error handler for non convergence of the secant solver """
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg


################################################################################
def SecantSolver(x0, x1, y, tol, itmax, xmin = None, xmax = None):
    """
    Finds a single the root of the 1-D equation y(x)
    using the Secant method.

    Inputs:
        x0     - First starting point
        x1     - Second starting point (typically x0 + dx)
        y      - Function of the form y = y(x)
        tol    - Tolerance to converge to
        itmax  - Maximum number of iterations
        xmin   - Minimum allowed x value
        xmax   - Maximum allowed x value

    Example:

    a = 2.
    b = 3.
    c = 1.
    d = 1.
    def y(x):
        return a*x**3 - b*x**2 + d*x - c
    
    >>> SecantSolver(x0 = 0., x1 = 0.1, y = y, tol = 1.e-5, itmax = 100)
    1.3981609724830093
    """
    xkm = copy.deepcopy(x0)
    xk  = copy.deepcopy(x1)

    ykm = y(x0)
    yk  = y(x1)

#    if isinstance(ykm, Unum):
#        if npy.isnan(ykm._value) or npy.isinf(ykm._value):
#            ykm = y(x0)
#            raise SecantError(message)
#    else:
#        if npy.isnan(ykm) or npy.isinf(ykm):
#            ykm = y(x0)
#            raise SecantError(message)


    xkp = x0  # Use the initial guess as the starting point
    ykp = ykm #
    it = 0 
    while abs(ykp) > tol and it < itmax:
            it += 1

            xkp = xk - (xk - xkm)/(yk - ykm)*yk
            
            if xmin is not None and xkp < xmin:
                xkp = xmin
                
            if xmax is not None and xkp > xmax:
                xkp = xmax

#            if isinstance(xkp, Unum):
#                if npy.isnan(xkp._value) or npy.isinf(xkp._value):
#                    raise SecantError(message)
            
            ykp = y(xkp)

#            if isinstance(ykp, Unum):
#                if npy.isnan(ykp._value) or npy.isinf(ykp._value):
#                    raise SecantError(message)
#            else:
#                if npy.isnan(ykp) or npy.isinf(ykp):
#                    raise SecantError(message)

            xkm = xk
            ykm = yk

            xk = xkp
            yk = ykp
            
            #print 'Secant', it, xkp, ykp

    if abs(ykp) > tol:
        message = 'SecantSolver did not converge to tol=' + str(tol) + ' in ' \
                    + str(itmax) + ' iterations.\n' + \
                    'Last iteration: x= ' + str(xkp) + ' y(x)=' + str(ykp)
        raise SecantError(message)
    
    message = 'SecantSolver did not compute a valid solution.\n' + \
              'Last iteration: x= ' + str(xkp) +' y(x)=' + str(ykp)

    if isinstance(xkp, _scalar):
        if npy.isnan(xkp.num) or npy.isinf(xkp.num):
            raise SecantError(message)
    else:
        if npy.isnan(xkp) or npy.isinf(xkp):
            raise SecantError(message)
    
    return xkp    

################################################################################
class BisectionError(Exception):
    """ Error handler for non convergence of the bisection solver """
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg

################################################################################
def BisectionSolver(x_lower, x_upper, y, tol, itmax):
    """
    Finds a single the root of the 1-D equation y(x)
    using the bisection method.

    Inputs:
        x_lower - Lower left bound of the equation
        x_upper - Upper right bound of the equation
        y       - Function of the form y = y(x)
        tol     - Tolerance to converge to
        itmax   - Maximum number of iterations

    Example:

    a = 2.
    b = 3.
    c = 1.
    d = 1.
    def y(x):
        return a*x**3 - b*x**2 + d*x - c
    
    >>> BisectionSolver(x_lower = 0., x_upper = 2, y = y, tol = 1.e-5, itmax = 100)
    1.398162841796875
    """
    xl = copy.deepcopy(x_lower)
    xr  = copy.deepcopy(x_upper)

    yl = y(xl)
    yr = y(xr)
    
    if yl*yr > 0:
        raise BisectionError("The function must have different signs at the end points for the Bisection solver.")

    it = 0
    yp = 1 + tol
    while abs(yp) > tol and it < itmax:
            it += 1
            xp = (xl + xr)/2; yp = y(xp)
            
            if yp*yl < 0:
                xr = xp; yr = yp
            else:
                xl = xp; yl = yp
            
            #print 'Bisection', it, xp, yp

    if abs(yp) > tol:
        raise SecantError('BisectionSolver did not converge to tol=' + str(tol) + ' in '
                    + str(itmax) + ' iterations.\n' +
                    'Last iteration value y(x)=' + str(yp))
    
    return xp    
################################################################################
if __name__ == '__main__':
    import doctest
    a = 2.
    b = 3.
    c = 1.
    d = 1.
    def y(x):
        return a*x**3 - b*x**2 + d*x - c

    doctest.testmod()

    x = SecantSolver(x0 = 0., x1 = 0.1, y = y, tol = 1.e-5, itmax = 100)
    print 'SecantSolver y(x) = 0 -> x = ', x
    print 'The answer should be 1.39816'
    
    x = BisectionSolver(x_lower = 0., x_upper = 2, y = y, tol = 1.e-5, itmax = 100)
    print 'BisectionSolver y(x) = 0 -> x = ', x
    print 'The answer should be 1.39816'
    
