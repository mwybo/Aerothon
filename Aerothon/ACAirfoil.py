from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import ARCDEG, RAD
from scalar.units import AsUnit
import numpy as npy
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline, UnivariateSpline
import pylab as pyl
import glob
import cmath as math
from RootSolvers import SecantSolver
from scipy import integrate
from scipy.optimize import brentq
from ACBase import ACBase

import os

################################################################################
class ACAirfoilError(Exception):
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg


################################################################################
class ACAirfoil:
    """
    A class for reading in a XFOIL airfoil file

    Attributes
        name     - Name of the airfoil
        TC       - Thickness to chord ratio
        x        - x coordinates
        y        - y coordinates
        Inverted - Determines if the airfoil should be inverted
    """

#===============================================================================
    def __init__(self, airfoildir, fext = '.txt', nalpha = 20):
        """
        Read in a XFOIL airfoil file and polar data

        Inputs:
            airfoildir - Directory containing coordiantes and polars for this airfoil
            nalpha     - Number of angles of attack to use to construct the interpolation
        """ 
        
        file = __file__
        if file.find(' ') == -1:
            #
            # Using an absolute path, the airfoil directory is one folder up 
            #
            rpath = __file__
            rpath = rpath[:rpath.rfind(os.sep)] #Remove this file name
            rpath = rpath[:rpath.rfind(os.sep)] #Remove the Aerothon folder
            rpath += os.sep
        
        else:
            #
            # Find the relative position of the airfoil directory (AVL can't use absolute path due to spaces in windows)
            #
            file = file[:file.rfind('Aerothon')]
            trunk = file[file[:-2].rfind(os.sep):] #Find the trunk of the repository (not everyone checks it out properly)
            
            path = os.path.abspath(os.path.curdir) #Find the relative path to the airfoils
            path = path[path.rfind(trunk)+len(trunk):]
            rpath = ''
            for i in range(path.count(os.sep)+1):
                rpath = os.pardir + os.sep + rpath

        AirfoilDir = rpath + 'Airfoils' + os.sep
        
        #
        # Add the relative path to the airfoil file names
        #
        afbase = AirfoilDir + airfoildir + os.sep + airfoildir
        airfoil    = afbase + '.dat'
        airfoilinv = afbase + '_inv.dat'
        polardir = AirfoilDir + airfoildir
        
        self.name        = airfoildir
        self.filename    = airfoil
        self.filenameinv = airfoilinv
        self.Inverted    = False
        self._ReadAirfoil(airfoil)
        self._CreatePolar(polardir, fext, nalpha)
        self._Area() # Compute the area of the airfoil

#===============================================================================
    def GetRealitveFilename(self):
        """
        Returns a relative path to the airfoil file name
        """
        #
        # Find the relative position of the airfoil directory (AVL can't use absolute path due to spaces in windows)
        #
        file = __file__
        file = file[:file.rfind('Aerothon')]
        trunk = file[file[:-2].rfind(os.sep):] #Find the trunk of the repository (not everyone checks it out properly)
        
        path = os.path.abspath(os.path.curdir) #Find the relative path to the airfoils
        path = path[path.rfind(trunk)+len(trunk):]
        rpath = ''
        for i in range(path.count(os.sep)+1):
            rpath = os.pardir + os.sep + rpath

        AirfoilDir = rpath + 'Airfoils' + os.sep
                
        return AirfoilDir + self.name + os.sep + self.name + '.dat'


#===============================================================================
    def _ReadAirfoil(self, airfoil):
        """
        Read in a XFOIL airfoil file

        Inputs:
            airfoil - File name with airfoil coordinates
        """
        #
        # Read in the x, y coordinates
        #
        try:
            f = open(airfoil,'r')
            self.airfoilname = f.readline()
            x = []
            y = []
            for line in f:
                if line == '/n':
                    continue
                if line.find('\t') is not -1:
                    xl = line.lstrip().partition('\t')
                else:
                    xl = line.lstrip().partition(' ')
                x.append(float(xl[0]))
                y.append(float(xl[2]))
        except IOError:
            print "While reading '", airfoil, "'"
            raise

        f.close()
        #
        # Clean up the name a bit
        #
        self.name = self.name.rstrip('\n')
        self.name = self.name.rstrip(' ')

        #
        # Split the coordinates into top and bottom
        #
        le = 0
        while x[le+1] < x[le]: le += 1

        #
        # Top coordinates
        #
        xtop = x[0:le+1]
        ytop = y[0:le+1]
                
        xtop.reverse()
        ytop.reverse()

        #
        # Bottom coordinates
        #
        xbot = x[le:]
        ybot = y[le:]
        
        def FixEdges(x,y):
            if x[0] > 0:
                x = [0] + x
                y = [0] + y
            if x[-1] < 1:
                x = x + [1]
                y = y + [y[-1]]
            return x, y
                
        xtop, ytop = FixEdges(xtop, ytop)
        xbot, ybot = FixEdges(xbot, ybot)

        top = interp1d(xtop, ytop, kind='linear')
        bot = interp1d(xbot, ybot, kind='linear')

        npt = 200
        TC = 0
        for i in range(1,npt):
            xi = float(i)/float(npt)
            TC = max(TC, top(xi) - bot(xi))

        self.x = npy.array(x)
        self.y = npy.array(y)

        self.TC = TC
        self.top = top
        self.bot = bot

#===============================================================================
    def _ReadPolar(self, filename):
        """
        Read in a XFOIL polar file

        Inputs:
            filename - File name with polar data

        Output:
            polar - File name with polar data

        Designed for XFOIL and XFLR5 polars
        """


        try:

            f = open(filename,'r')
            line = f.readline()
            #
            # Find the header line before the data
            #
            while line.find('-------') == -1:
                # Grab the Reynolds number from the file
                if line.find('Re =') != -1:
                    line = line[line.find('Re =')+4:]
                    line = line[:line.find('Ncrit')]
                    line = line.replace(' ','')
                    Re = float(line)

                line = f.readline()
            #
            # Read in the polar
            #
            alpha = []
            Cl = []
            Cd = []
            Cdp = []
            Cm = []
            Clmax = -99
            Clmin = 99
            AlphaClmax = 0
            AlphaClmin = 0
            for line in f:
                vals = line.lstrip()
                vals = vals.split('  ')
                if len(vals) == 1:
                    break #Just in case there is a newline at the end of the file

                alpha.append(float(vals[0]))
                Cl.append(float(vals[1]))
                Cd.append(float(vals[2]))
                Cdp.append(float(vals[3]))
                Cm.append(float(vals[4]))
                
                if Cl[-1] > Clmax:
                    Clmax = Cl[-1]
                    AlphaClmax = alpha[-1]
                    
                if Cl[-1] < Clmin:
                    Clmin = Cl[-1]
                    AlphaClmin = alpha[-1]

        except IOError, (errno, strerror):
            print "While reading '", self.filename, "'"
            raise
        except ValueError:
            print "Data in the airfoil ", self.filename, " is corrupt..."
            f.close()
            raise

        f.close()

        #print 'alpha ', alpha
        #print 'Cl', Cl
        #print 'Cd', Cd
        #print 'Cdp', Cdp
        #print 'Cm', Cm

        class PolarClass:
            #
            # Add the methods so the polars can be sorted
            #
            def __lt__(self, other):
                return self.Re < other.Re
            def __gt__(self, other):
                return self.Re > other.Re

        polar = PolarClass()

        alpha = npy.array(alpha)
        Cl = npy.array(Cl)
        Cd = npy.array(Cd)
        Cm = npy.array(Cm)
        Cdp = npy.array(Cdp)

        #
        # Save off all the inAsUnition about this polar
        #
        polar.Re = Re
        polar.amax = npy.max(alpha)
        polar.amin = npy.min(alpha)

        polar.Clmax = npy.max(Cl)
        polar.Clmin = npy.min(Cl)

        polar.AlphaClmax = AlphaClmax
        polar.AlphaClmin = AlphaClmin
        
        polar.alpha = alpha
        polar.Cl = Cl
        polar.Cd = Cd
        polar.Cm = Cm
        

        #
        # Perform a linear Least Squares Curvefit through the polar
        # to get an approximate average slope. The range is limited to
        # the range of alphas between Clmin and Clmax
        #
        def PolarSlope(C):
 #           a = (polar.amax + polar.amin)/2.
            imin = npy.argmin(C)
            imax = npy.argmax(C)
            imin, imax = min(imin, imax), max(imin, imax)
            ai = alpha[imin:imax]
            Ci = C[imin:imax]
            #
            # Do the linear Least Squares Curvefit through the entire polar
            # and return it's derivative
            #
            Cp = npy.polyfit(ai, Ci, 1)
            return npy.polyder(Cp)

        polar.dCl_da = PolarSlope(Cl)
        polar.dCd_da = PolarSlope(Cd)
        polar.dCm_da = PolarSlope(Cm)

        return polar

#===============================================================================
    def _CreatePolar(self, polardir, fext, nalpha):
        """
        Read in a XFOIL polar file

        Inputs:
            polardir - A directory with polars for the airfoil
            fext     - File extension (.dat of .txt or ...)
            naplha   - Number of angles of attack used to construct the interpolation data
        """

        #
        # Read all polar files in the directory
        #
        pdir = polardir + '/*.' + fext.lstrip('.')
        pfiles = glob.glob(pdir)
        
        if len(pfiles) == 0:
            raise ACAirfoilError("Could not find any airfoil polars in " + pdir)
        
        self.polars = []
        polars = self.polars
        for pf in pfiles:
            pfile = pf
            polars.append(self._ReadPolar(pfile))

        #
        # Make sure the polars are sorted by Reynolds number
        #
        polars.sort()

        #
        # Get the bounds for the angles of attach
        #
        amax = -90
        amin = 90
        for polar in polars:
            amax = max(amax, polar.amax)
            amin = min(amin, polar.amin)
        #
        # Save off the range of angles attack to
        # enable the AlphaRange method
        #
        self.AlphaMax = amax
        self.AlphaMin = amin

        #
        # Create the angles of attack
        #
        alpha = self.AlphaRange() / ARCDEG

        #
        # Construct the 2D interpolation arrays
        #
        Cl     = npy.empty((len(polars),len(alpha)))
        Cd     = npy.empty((len(polars),len(alpha)))
        Cm     = npy.empty((len(polars),len(alpha)))

        Rey         = npy.empty(len(polars))
        dCl_da      = npy.empty(len(polars))
        dCd_da      = npy.empty(len(polars))
        dCm_da      = npy.empty(len(polars))
        Clmax       = npy.empty(len(polars))
        Clmin       = npy.empty(len(polars))
        AlphaClmax  = npy.empty(len(polars))
        AlphaClmin  = npy.empty(len(polars))

        for j in range(len(polars)):
            polar = polars[j]
            Re           = polar.Re
            Rey[j]       = Re
            dCl_da[j]    = polar.dCl_da
            dCd_da[j]    = polar.dCd_da
            dCm_da[j]    = polar.dCm_da
            Clmax[j]     = polar.Clmax
            Clmin[j]     = polar.Clmin
            AlphaClmax[j] = polar.AlphaClmax
            AlphaClmin[j] = polar.AlphaClmin
            for i in range(len(alpha)):
                #
                # Use npy.interp here because it "extrapolates" the end values
                #
                Cl[j,i] = npy.interp(alpha[i], polar.alpha, polar.Cl)
                Cd[j,i] = npy.interp(alpha[i], polar.alpha, polar.Cd)
                Cm[j,i] = npy.interp(alpha[i], polar.alpha, polar.Cm)

        #
        # Save off the interpolation data
        #
        self.Clint = RectBivariateSpline(Rey, alpha, Cl, kx = 1, ky = 1)
        self.Cdint = RectBivariateSpline(Rey, alpha, Cd, kx = 1, ky = 1)
        self.Cmint = RectBivariateSpline(Rey, alpha, Cm, kx = 1, ky = 1)

        #
        # Save off some other data
        #
        self.Rey      = Rey
        self.AlphaMax = amax
        self.AlphaMin = amin
        self.ReMin    = min(Rey)
        self.ReMax    = max(Rey)

        self.Clmaxint = interp1d(Rey, Clmax, bounds_error = True)
        self.Clminint = interp1d(Rey, Clmin, bounds_error = True)

        self.AlphaClmaxint = interp1d(Rey, AlphaClmax, bounds_error = True)
        self.AlphaClminint = interp1d(Rey, AlphaClmin, bounds_error = True)

        self.dCl_daint = interp1d(Rey, dCl_da, bounds_error = True)
        self.dCd_daint = interp1d(Rey, dCd_da, bounds_error = True)
        self.dCm_daint = interp1d(Rey, dCm_da, bounds_error = True)

#===============================================================================
    def _CheckRe(self, Re):
        "Checks that a Reynolds number is within the interpolation range"
        
        errmsg = "Reynolds number " + str(Re) + " is outside the range of " + self.name + '.' + \
                 "\nRe min and max are: " + str(self.ReMin) + ", " + str(self.ReMax)
        try:
            if any(Re > self.ReMax) or  any(Re < self.ReMin):
                raise ACAirfoilError(errmsg)
        except:
            if Re > self.ReMax or Re < self.ReMin:
                raise ACAirfoilError(errmsg)

#===============================================================================
    def _Area(self):
        "Computes the non-dimensional area of the airfoil"
        self.Area = integrate.romberg(function=self.Thickness, a=min(self.x), b=max(self.x), divmax = 100, tol=1.e-2)

#===============================================================================
    def GetFilename(self):
        "Returns the filename of the airfoil"
        if self.Inverted:
            return self.filenameinv
        else:
            return self.filename

#===============================================================================
    def Lower(self, x = None):
        "Gives the lower y coordinate at x"
        
        if x is None:
            if self.Inverted:
                return -max(self.y)
            else:
                return min(self.y)
                
        if self.Inverted:
            return -self.top(x)
        else:
            return self.bot(x)

#===============================================================================
    def Upper(self, x = None):
        "Gives the upper y coordinate at x"
        
        if x is None:
            if self.Inverted:
                return -min(self.y)
            else:
                return max(self.y)

        if self.Inverted:
            return -self.bot(x)
        else:
            return self.top(x)

#===============================================================================
    def Thickness(self, x = None):
        "Gives thickness at x"
        if x is None:
            return self.TC
        else:
            return self.top(x) - self.bot(x)

#===============================================================================
    def AlphaZeroLift(self, Re):
        """
        Provides the angle of attack for zero lift

        Inputs:
            Re - Reynolds number
        """
        self._CheckRe(Re)
        def y(alpha):
            return self.Clint(Re, alpha)[0,0]

        #
        # Zero lift typically occurs near alpha = 0
        #
        x0 = 0.
        x1 = 1.
        try:
#            a0 = SecantSolver(x0 = x0, x1 = x1, y = y, tol = 1.e-5, itmax = 100) * ARCDEG
            a0 = brentq(f = y, a = -30, b = 30, rtol = 1.e-5, maxiter = 500) * ARCDEG
        except:
            a0 = 0*ARCDEG
        
        if self.Inverted:
            return a0
        else:
            return -a0


#===============================================================================
    def Clmax(self, Re):
        """
        Provides an interpolated 2D airfoil Clmax value

        Inputs:
            Re - Reynolds number
        """
        self._CheckRe(Re)
        if self.Inverted:
            return -self.Clminint(float(Re))
        else:
            return self.Clmaxint(float(Re))

#===============================================================================
    def Clmin(self, Re):
        """
        Provides an interpolated 2D airfoil Clmin value

        Inputs:
            Re - Reynolds number
        """
        self._CheckRe(Re)
        if self.Inverted:
            return -self.Clmaxint(Re)[0]
        else:
            return self.Clminint(Re)[0]
        
#===============================================================================
    def AlphaClmax(self, Re):
        """
        Provides an interpolated angle of attack for the 2D airfoil Clmax value

        Inputs:
            Re - Reynolds number
        """
        self._CheckRe(Re)
        if self.Inverted:
            return -self.AlphaClminint(float(Re))*ARCDEG
        else:
            return self.AlphaClmaxint(float(Re))*ARCDEG

#===============================================================================
    def AlphaClmin(self, Re):
        """
        Provides an interpolated angle of attack for the 2D airfoil Clmin value

        Inputs:
            Re - Reynolds number
        """
        self._CheckRe(Re)
        if self.Inverted:
            return -self.AlphaClmaxint(Re)*ARCDEG
        else:
            return self.AlphaClminint(Re)*ARCDEG
                
#===============================================================================
    def _GetAlpha(self, alpha2d):
        """
        Retrieves the list of angles of attack for a polar properly sorted for a 
        polar.
        """
        if self.Inverted:
            try:
                a = npy.array([-alpha2d[i] / ARCDEG for i in range(-1,-len(alpha2d)-1,-1)]) * ARCDEG
            except:
                a = -alpha2d
        else:
            a = alpha2d
            
        return a if hasattr(a,'__len__') else [a]

#===============================================================================
    def _GetCoef(self, alpha2d, Re, Cint, Cname):
        """
        Retrieves the interpolated polar.
        
        Inputs:
            alpha2d - A range of 2D angels of attack
            Cint    - An interpolation of a ceffocient (Cl, Cd, or Cm)
            Re      - Reynolds number for the interpolation
            Cname   - The name of the coefficient used in an error message
        """
        self._CheckRe(Re)
        alpha = self._GetAlpha(alpha2d)
        
        if not hasattr(Re,'__len__'):
            Re = npy.array([Re])
        
        try:
            C = []
            if self.Inverted:
                for a in alpha:
                    C.append( -Cint(Re, a / ARCDEG)[:,0] )
                C = npy.array(C)
                C = npy.array([[C[i,j] for i in range(-1,-len(C[:,j])-1,-1)] for j in range(len(Re))])
            else:
                for a in alpha:
                    C.append( Cint(Re, a / ARCDEG)[:,0] )
                C = npy.transpose(npy.array(C))
                
            C = C[0] if len(C) == 1 else C
            return C[0] if len(C) == 1 else C
        
        except AssertionError:
            raise ACAirfoilError(self.name + " failed to interpolate "+ Cname + " for AoA=" + str(a) + ", Re=" + str(Re))

#===============================================================================
    def Cl(self, alpha2d, Re):
        """
        Provides an interpolated 2D airfoil Cl value

        Inputs:
            alpha2d - Angle of attack
            Re      - Reynolds number
        """
        return self._GetCoef(alpha2d, Re, self.Clint, 'Cl')

#===============================================================================
    def Cd(self, alpha2d, Re):
        """
        Provides an interpolated 2D airfoil Cd value

        Inputs:
            alpha2d - Angle of attack
            Re      - Reynolds number
        """
        Cd = self._GetCoef(alpha2d, Re, self.Cdint, 'Cd')
        #
        # Drag must be flipped back... oh well
        #
        if self.Inverted:
            return -Cd
        else:
            return Cd

#===============================================================================
    def Cm(self, alpha2d, Re):
        """
        Provides an interpolated 2D airfoil Cm value

        Inputs:
            alpha2d - Angle of attack
            Re      - Reynolds number
        """
        return self._GetCoef(alpha2d, Re, self.Cmint, 'Cm')

#===============================================================================
    def dCl_da(self, Re, alpha2d = None):
        """
        Provides an interpolated 2D airfoil Cl value.
        If an angle of attack is given the derivative is evaluated
        at that angle of attack, otherwise a derivative of the linear
        least squares curvefit between the min and max is given.

        Inputs:
            Re      - Reynolds number
            alpha2d - Angle of attack
        """
        self._CheckRe(Re)
        if alpha2d is None:
            return self.dCl_daint(Re) / ARCDEG
        else:
            return self._da(alpha2d, Re, self.Clint)

#===============================================================================
    def dCd_da(self, Re, alpha2d = None):
        """
        Provides an interpolated 2D airfoil Cd value.
        If an angle of attack is given the derivative is evaluated
        at that angle of attack, otherwise a derivative of the linear
        least squares curvefit between the min and max is given.

        Inputs:
            Re      - Reynolds number
            alpha2d - Angle of attack
        """
        self._CheckRe(Re)
        if alpha2d is None:
            return self.dCd_daint(Re) / ARCDEG
        else:
            return self._da(alpha2d, Re, self.Cdint)
        
#===============================================================================
    def dCm_da(self, Re, alpha2d = None):
        """
        Provides an interpolated 2D airfoil Cm value.
        If an angle of attack is given the derivative is evaluated
        at that angle of attack, otherwise a derivative of the linear
        least squares curvefit between the min and max is given.

        Inputs:
            Re      - Reynolds number
            alpha2d - Angle of attack
        """
        self._CheckRe(Re)
        if alpha2d is None:
            return self.dCm_daint(Re) / ARCDEG
        else:
            return self._da(alpha2d, Re, self.Cmint)

#===============================================================================
    def _da(self, alpha2d, Re, interp):
        """
        Calculates a alpha derivatives for the RectBivariateSpline interp.
        This is a hack since parder has not been added to the
        interfaces of fitpack yet and derivatives cannot be computed
        with a BivariateSpline

        Inputs:
            alpha2d - Angle of attack
            Re      - Reynolds number
            interp  - RectBivariateSpline with interpolation data
        """

        if self.Inverted:
            a0 = -alpha2d[1] / ARCDEG
            a1 = -alpha2d[0] / ARCDEG
        else:
            a0 = alpha2d[0] / ARCDEG
            a1 = alpha2d[1] / ARCDEG
            
        #
        # Return the derivative
        #
        return (interp(Re, a1) - interp(Re, a0))[0][0] / (a1-a0) / ARCDEG

#===============================================================================
    def AlphaRange(self, nalpha = 30):
        """
        Creates an evenly spaced range of angles of attack
        for plotting purposes
        """
        #
        # Create a range of angles of attack
        #
        amax = self.AlphaMax
        amin = self.AlphaMin
        
        if self.Inverted:
            return npy.linspace(-amax, -amin, nalpha) * ARCDEG
        else:
            return npy.linspace(amin, amax, nalpha) * ARCDEG

#===============================================================================
    def XYAtAngle(self, alpha2d = 0 * ARCDEG):
        """
        Returns the airfoils coordinates rotated about the leading edge to the
        desired angle of attack

        Inputs:
            alpha2d - Angle of attack of the airfoil
        """
        #
        # Create a range of angles of attack
        #
        if self.Inverted:
            a = -alpha2d / RAD
        else:
            a = alpha2d / RAD
            
        xnew =  self.x * math.cos(a).real + self.y * math.sin(a).real
        ynew = -self.x * math.sin(a).real + self.y * math.cos(a).real

        if self.Inverted:
            return xnew, -ynew
        else:
            return xnew, ynew

#===============================================================================
    def PlotPolar(self, fig = 1, nalpha = 30, Re = None):
        """
        Plot the interpolated polars as a 3D surface

        Inputs:
            nalpha - Number of angles of attack
        """
        from mpl_toolkits.mplot3d.axes3d import Axes3D
        import matplotlib.pyplot as plt
        #
        # Create the angles of attack
        #
        alpha = self.AlphaRange(nalpha)

        #
        # Create the reynolds numbers
        #
        if Re is None:
            Rey = self.Rey #npy.linspace(self.Rey[0], self.Rey[-1], 20)
        else:
            Rey = Re
        
        Cl = self.Cl(alpha, Rey)
        Cd = self.Cd(alpha, Rey)
        Cm = self.Cm(alpha, Rey)

        alpha = alpha / ARCDEG

        figure = plt.figure(fig)
        ax221 = figure.add_subplot(2, 2, 1, projection='3d')
        ax222 = figure.add_subplot(2, 2, 2, projection='3d')
        ax223 = figure.add_subplot(2, 2, 3, projection='3d')
        ax224 = figure.add_subplot(2, 2, 4, projection='3d')

        alpha, Rey = npy.meshgrid(alpha, Rey)

        colormap = pyl.cm.hot

        ax221.plot_surface(alpha,Rey,Cl,rstride=1,cstride=1,cmap=colormap)
        ax222.plot_surface(alpha,Rey,Cd,rstride=1,cstride=1,cmap=colormap)
        ax223.plot_surface(alpha,Rey,Cm,rstride=1,cstride=1,cmap=colormap)
        ax224.plot_surface(Cd,Rey,Cl,rstride=1,cstride=1,cmap=colormap)
                        
#        for polar in self.polars:
#            Res = [polar.Re]*len(polar.alpha)
#            ax221.plot(polar.alpha,Res,polar.Cl,'r')
#            ax222.plot(polar.alpha,Res,polar.Cd,'r')
#            ax223.plot(polar.alpha,Res,polar.Cm,'r')
#            ax224.plot(polar.Cd   ,Res,polar.Cl,'r')
            
        
        ax221.set_xlabel(r'$\alpha[^o]$'); ax221.set_ylabel('Re'); ax221.set_zlabel(r'$C_l$')
        ax222.set_xlabel(r'$\alpha[^o]$'); ax222.set_ylabel('Re'); ax222.set_zlabel(r'$C_d$')
        ax223.set_xlabel(r'$\alpha[^o]$'); ax223.set_ylabel('Re'); ax223.set_zlabel(r'$C_m$')
        ax224.set_xlabel(r'$C_d$');        ax224.set_ylabel('Re'); ax224.set_zlabel(r'$C_l$')

        Remin = npy.min(Rey)
        Remax = npy.max(Rey)
        
        if len(Rey) > 1:
            title = self.name + '  ' + str(Remin) + ' < Re < ' + str(Remax)
        else:
            title = self.name + '  ' + ' Re ' + str(Remax)
            
        pyl.annotate(title, xy=(.025, .975),
                     xycoords='figure fraction',
                     horizontalalignment='left', verticalalignment='top',
                     fontsize=20)

#===============================================================================
    def Plot2DPolar(self, fig = 1, nalpha = 30, Re = None):
        """
        Plot the interpolated polars

        Inputs:
            nalpha - Number of angles of attack
        """
        #
        # Create the angles of attack
        #
        alpha = self.AlphaRange(nalpha)

        #
        # Create the reynolds numbers
        #
        if Re is None:
            Rey = npy.linspace(self.Rey[0], self.Rey[-1], 20)
        else:
            try:
                a = len(Re)
                Rey = Re
            except:         
                Rey = [Re]
        
        Cl = self.Cl(alpha, Rey)
        Cd = self.Cd(alpha, Rey)
        Cm = self.Cm(alpha, Rey)

        alpha = alpha / ARCDEG

        figure = pyl.figure(fig)
        ax221 = figure.add_subplot(321)
        ax222 = figure.add_subplot(322)
        ax223 = figure.add_subplot(323)
        ax224 = figure.add_subplot(324)
        ax225 = figure.add_subplot(325)
        
        if len(Rey) == 1:
            ax221.plot(alpha,Cl)
            ax222.plot(alpha,Cd)
            ax223.plot(alpha,Cm)
            ax224.plot(Cd,Cl)
            ax225.plot(Cl,Cl/Cd)
        else:
            for i in xrange(len(Rey)):
                ax221.plot(alpha,Cl[i,:])
                ax222.plot(alpha,Cd[i,:])
                ax223.plot(alpha,Cm[i,:])
                ax224.plot(Cd[i,:],Cl[i,:])
                ax225.plot(Cl[i,:],Cl[i,:]/Cd[i,:])
                                        
        ax221.set_xlabel(r'$\alpha[^o]$'); ax221.set_ylabel(r'$C_l$')
        ax222.set_xlabel(r'$\alpha[^o]$'); ax222.set_ylabel(r'$C_d$')
        ax223.set_xlabel(r'$\alpha[^o]$'); ax223.set_ylabel(r'$C_m$')
        ax224.set_xlabel(r'$C_d$');        ax224.set_ylabel(r'$C_l$')
        ax225.set_xlabel(r'$C_l$');        ax225.set_ylabel(r'$C_l/C_d$')

        Remin = npy.min(Rey)
        Remax = npy.max(Rey)
        
        if len(Rey) > 1:
            title = self.name + '  ' + str(Remin) + ' < Re < ' + str(Remax)
        else:
            title = self.name + '  ' + ' Re ' + str(Remax)
            
        pyl.annotate(title, xy=(.025, .975),
                     xycoords='figure fraction',
                     horizontalalignment='left', verticalalignment='top',
                     fontsize=20)

#===============================================================================
    def PlotReportPolars(self, fig = 1, nalpha = 30, Re = None):
        """
        Plot the interpolated polars

        Inputs:
            nalpha - Number of angles of attack
        """
        #
        # Create the angles of attack
        #
        alpha = self.AlphaRange(nalpha)

        #
        # Create the reynolds numbers
        #
        if Re is None:
            Rey = npy.linspace(self.Rey[0], self.Rey[-1], 20)
        else:
            try:
                a = len(Re)
                Rey = Re
            except:         
                Rey = [Re]
        
        Cl = self.Cl(alpha, Rey)
        Cd = self.Cd(alpha, Rey)
        Cm = self.Cm(alpha, Rey)

        alpha = alpha / ARCDEG

        figure = pyl.figure(fig)
        ax221 = figure.add_subplot(311)
        ax222 = figure.add_subplot(312)
        ax223 = figure.add_subplot(313)
        
        if len(Rey) == 1:
            ax221.plot(alpha,Cl)
            ax222.plot(alpha,Cd)
            ax223.plot(alpha,Cm)
        else:
            for i in xrange(len(Rey)):
                ax221.plot(alpha,Cl[i,:])
                ax222.plot(alpha,Cd[i,:])
                ax223.plot(alpha,Cm[i,:])
                                        
        ax221.set_ylabel(r'$C_l$')
        pyl.setp( ax221.get_xticklabels(), visible=False)
        ax222.set_ylabel(r'$C_d$')
        pyl.setp( ax222.get_xticklabels(), visible=False)
        ax223.set_ylabel(r'$C_m$')
        ax223.set_xlabel(r'$\alpha[^o]$'); 
                
#===============================================================================
    def PlotAirfoil(self, fig = 1, Alpha2d = 0*ARCDEG, subfig = None, yo = None, clr = 'b'):
        """
        Plots the arifoil at the given angle of attack

        Inputs:
            fig     - Figure to plot to
            Alpha2d - Angle of attack to plot the polar at
            subfig  - Subplot number
            yo      - The y origin of the airfoil (good for plotting multiple airfoils)
            clr     - The color the airfoil is plotted with
        """
        
        pyl.figure(fig)
        if subfig is not None:
            pyl.subplot(subfig)
        pyl.axis('equal')
        x, y = self.XYAtAngle(Alpha2d)
        
        if yo is not None:
            for i in range(len(y)):
                y[i] = y[i] + yo
        
        pyl.plot(x, y, clr)

################################################################################
if __name__ == '__main__':
    #
    # Read in an arifoil
    #
#    af = ACAirfoil('S1223')
#    af = ACAirfoil('NACA0012')
    af = ACAirfoil('e423')

    Re      = 235643
    alpha2d = (5 * ARCDEG, 6 *ARCDEG)
    a0      = af.AlphaZeroLift(Re)
    af.Inverted = False

    print 'Alpha Range      : ', af.AlphaRange()
    print 'Name             : ', af.name
    print 'Thicness to chord: ', af.TC
    print 'dCl_da           : ', af.dCl_da(Re, alpha2d)
    print 'dCm_da           : ', af.dCm_da(Re, alpha2d)
    print 'AlphaZeroLift    : ', AsUnit(a0 , "deg")
    print 'Cl(a0)           : ', af.Cl(a0,Re)
    print 'Area             : ', af.Area

    af.Plot2DPolar()
    af.PlotPolar(fig=2)
#    print af.AlphaRange()
    af.PlotAirfoil(3, 0*ARCDEG)
    pyl.show()


