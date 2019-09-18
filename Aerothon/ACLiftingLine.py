"""
A class for analyzing lifting surfaces with lifting line theory

A good reference can be found at
http://www.desktop.aero/appliedaero/potential3d/liftingline.html

"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACBase
import numpy as npy
import pylab as pyl
from scipy.interpolate import interp1d
from scalar.units import IN, ARCDEG, RAD, LBF, FT, SEC
from scalar.units import AsUnit
from math import sqrt, sin, cos, atan, pi
import AeroUtil
from ACBase import Angle, Velocity, Length, Unitless

################################################################################
class ACLiftingLine(ACBase):
    """
    A class for lifting line theory analasys
    """
#===============================================================================
    def __init__(self, LiftSurf):
        super(ACLiftingLine, self).__init__()
        
        UnitList = self.UnitList
        NoneList = self.NoneList
        
        LiftSurf.Parents.append( self )
        self.__dict__['LiftSurf'] = LiftSurf
        
        self.__dict__['Alpha3d'] = 0*ARCDEG   ; UnitList['Alpha3d']  = Angle
        self.__dict__['nSpan']   = 51         ; UnitList['nSpan']    = Unitless
        self.NoneList['Vinf']    = None       ; UnitList['Vinf']     = Velocity
        self.NoneList['Alt']     = None       ; UnitList['Alt']      = Length
        self.__dict__['Omega']   = 0.1        ; UnitList['Omega']    = Unitless
        self.__dict__['itmax']   = 100        ; UnitList['itmax']    = Unitless
        self.__dict__['tol']     = 1e-7       ; UnitList['tol']      = Unitless
        
        self.param.nSpan = 0
       
#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACLiftingLine, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'Vinf':
            return self._CalcVinf()
        elif key == 'Alt':
            return self._CalcAlt()
        #
        # If no calculation exist return None
        #
        return None
         
#===============================================================================
    def _CalcVinf(self):
        """
        Computes the reference velocity as the lift off velocity
        """
        return self.LiftSurf.GetV_LO()
            
#===============================================================================
    def _CalcAlt(self):
        """
        Computes the reference velocity as the lift off altitude
        """
        return self.LiftSurf.GetAlt_LO()
             
#===============================================================================
    def Re(self, y):
        """
        Computes the local Reynolds number
        
        Inputs:
            y - Spanwise location
        """
        Alt = self.Alt
        V   = self.Vinf
        rho = AeroUtil.DensAlt(Alt)
        mu  = AeroUtil.Mu(AeroUtil.TempAlt(Alt))
        c   = self.LiftSurf.Chord(y)

        return rho * V * c / mu
        
#===============================================================================
    def InitEllipticCirculation(self):
        """
        Initializes the circulation for an elliptic lift distribution
        """
        y = 0*IN
        Re = self.Re(y)
        Cl = self.LiftSurf.Cl(self.Alpha3d, y, Re)
        Cr = self.LiftSurf.Chord(y)
        S  = self.LiftSurf.S
        b  = self.LiftSurf.b
        Vinf = self.Vinf
        
        Gam0 = 2*Vinf*S*Cl/(b*pi) /(FT**2/SEC)
        Gamma = npy.empty(self.nSpan)
        #y = self.param.y = b/2*npy.cos( npy.linspace(pi, 0,self.nSpan) )
        if self.LiftSurf.FullWing:
            y = self.param.y = npy.linspace(-b/2/IN, b/2/IN,self.nSpan)*IN
        else:
            y = self.param.y = npy.linspace(0, b/IN,self.nSpan)*IN
            
        for i in xrange(self.nSpan):
            Gamma[i] = Gam0*self.LiftSurf.EllipChord(y[i])/Cr

        self.param.Gamma   = Gamma * (FT**2/SEC)
        self.param.Alpha_i = npy.zeros(self.nSpan)*ARCDEG

#===============================================================================
    def Alpha_i(self, y):
        """
        Calculates the induced angle of attach
        
        Inputs:
            y - Spanwise location
        """
        if self.dirty: self.Refresh()
        
        return self.param.Alpha_i(y/IN)*ARCDEG

#===============================================================================
    def Cl(self, y):
        """
        Calculates the local lift coefficient
        
        Inputs:
            y - Spanwise location
        """
        if self.dirty: self.Refresh()
        
        return self.param.Cl(y/IN)

#===============================================================================
    def CL(self):
        """
        Calculates the total lift coefficient
        """
        if self.dirty: self.Refresh()
        
        return self.param.CL

#===============================================================================
    def CDi(self):
        """
        Calculates the total induced drag
        """
        if self.dirty: self.Refresh()
        
        return self.param.CDi

#===============================================================================
    def CDv(self):
        """
        Calculates the total viscous drag
        """
        if self.dirty: self.Refresh()
        
        return self.param.CDv

#===============================================================================
    def Cd(self, y):
        """
        Calculates the local airfoil drag
        
        Inputs:
            y - Spanwise location
        """
        if self.dirty: self.Refresh()
        
        return self.param.Cd(y/IN)

#===============================================================================
    def CD(self):
        """
        Calculates the total drag
        """
        if self.dirty: self.Refresh()
        
        return self.CDv() + self.CDi()

#===============================================================================
    def OswaldEff(self):
        """
        Calculates the Oswald efficiency
        """
        if self.dirty: self.Refresh()
        
        return self.CL()**2/(pi*self.LiftSurf.AR*self.CDi())

#===============================================================================
    def _Alpha_i(self, y):
        """
        Calculates the induced angle of attach for the current circulation
        
        Inputs:
            y - Spanwise location
        """
        
        y0 = y
        Vinf = self.Vinf
        Gamma = self.param.Gamma
        y     = self.param.y
        nSpan = self.param.nSpan
        
        #dy = (y[1:nSpan] - y[0:nSpan-1])
        dgdy = (Gamma[1:nSpan] - Gamma[0:nSpan-1]) #/dy
        ymid = y0 - (y[1:nSpan] + y[0:nSpan-1])/2
        
        integral = npy.sum( dgdy/ymid /(FT/SEC) ) *FT/SEC #*dy )
        
        return atan( -1/(4*pi*Vinf)*integral )*RAD

#===============================================================================
    def _CL(self):
        """
        Calculates the total lift of the wing
        """
        Vinf  = self.Vinf
        S     = self.LiftSurf.S
        y     = self.param.y
        Gamma = self.param.Gamma
        nSpan = self.nSpan
        
        dy   = (y[1:nSpan] - y[0:nSpan-1])
        Gavg = (Gamma[1:nSpan] + Gamma[0:nSpan-1])/2
        
        Integral = npy.sum( Gavg*dy / (FT**3/SEC) ) *FT**3/SEC
        
        return 2/(Vinf*S)*Integral
        
#===============================================================================
    def _CDi(self, Alpha_i):
        """
        Calculates the total induced drag of the wing
        
        Inputs:
            Alpha_i - Induced angle of attack
        """
        Vinf  = self.Vinf
        S     = self.LiftSurf.S
        y     = self.param.y
        Gamma = self.param.Gamma
        nSpan = self.nSpan
        
        Alpha_i = Alpha_i/RAD
        
        dy   = (y[1:nSpan] - y[0:nSpan-1])
        AiGavg = (Gamma[1:nSpan]*npy.sin(Alpha_i[1:nSpan]) + Gamma[0:nSpan-1]*npy.sin(Alpha_i[0:nSpan-1]))/2
        
        Integral = npy.sum( AiGavg*dy / (FT**3/SEC) ) *FT**3/SEC
        
        return -2/(Vinf*S)*Integral

#===============================================================================
    def _Cd(self, Alpha_i, Re):
        """
        Calculates the total drag from the airfoil
        
        Inputs:
            Alpha_i - Induced angle of attack
        """
        Vinf  = self.Vinf
        b     = self.LiftSurf.b
        y     = self.param.y
        nSpan = self.nSpan
        Alpha3d = self.Alpha3d
        
        Cd = npy.empty(nSpan)
        for i in xrange(nSpan):
            Alpha_eff = Alpha3d + Alpha_i[i]
            Cd[i] = self.LiftSurf.Cd(Alpha_eff, y[i], Re[i])
                
        dy   = (y[1:nSpan] - y[0:nSpan-1])
        Cdavg = (Cd[1:nSpan] + Cd[0:nSpan-1])/2
        Integral = npy.sum( Cdavg*dy / FT ) *FT
        
        return 1/b*Integral, Cd
    
#===============================================================================
    def _UpdateGamma(self):
        """
        Calculates the circulation for the lifting surface
        """
        if self.param.nSpan != self.nSpan:
            self.InitEllipticCirculation()
            self.param.nSpan = self.nSpan
        
        Alpha3d = self.Alpha3d
        Vinf    = self.Vinf
        nSpan   = self.nSpan
        y       = self.param.y
        Gamma   = self.param.Gamma
        w       = self.Omega
        itmax   = self.itmax
        
        c = self.LiftSurf.Chord(y)
        Re = self.Re(y)
        
        def Update(start, end, inc):
            L2 = 0
            for i in xrange(start, end, inc):
                yi = y[i]
                Alpha_i = self._Alpha_i(yi)
                Alpha_eff = Alpha3d + Alpha_i
                
                Cl = self.LiftSurf.Cl(Alpha_eff, yi, Re[i])
                
                G_new = 0.5*Cl*c[i]*Vinf*w + (1-w)*Gamma[i]
                
                L2 += ( (G_new - Gamma[i])/(FT**2/SEC) )**2
                
                Gamma[i] = G_new
        
            return L2
        
        L2 = self.tol + 1
        it = 0
        
        while L2 > self.tol and it < itmax:
            it += 1
            L2 = Update(1, nSpan-1, 1)
            #if L2 > self.tol:
            #    L2 = Update( nSpan-2, 1, -1)
                
        #print 'Converged in : ', it
        
        Alpha_i = npy.empty(nSpan)
        for i in xrange(nSpan):
            Alpha_i[i] = self._Alpha_i(y[i]) / ARCDEG
            
        Alpha_i[0]       = Alpha_i[1]
        Alpha_i[nSpan-1] = Alpha_i[nSpan-2]
        
        Cl = 2*Gamma/(Vinf*c)
        
        self.param.Alpha_i = interp1d(y/IN, Alpha_i, kind='linear')
        self.param.Cl      = interp1d(y/IN, Cl, kind='linear')
        self.param.CL      = self._CL()
        self.param.CDi     = self._CDi(Alpha_i*ARCDEG)
        self.param.CDv, Cd = self._Cd(Alpha_i*ARCDEG, Re)
        self.param.Cd      = interp1d(y/IN, Cd, kind='linear')

        if it == itmax:
            print "WARNING: ", self.LiftSurf.name, " - Lifting line failed to converge to ", self.tol, " in ", itmax, " iterations."
            print "         Final L2 norm : ", L2, ' | oswald eff : ', self.OswaldEff()

#===============================================================================
    def Refresh(self):
        """
        Updates all the internal quantities.

        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACLiftingLine, self).Refresh()
        self.param.refreshing = True

        #
        # Compute the circulation of the wing
        #
        self._UpdateGamma()
                
        self.param.refreshing = False        

################################################################################
if __name__ == '__main__':
    from ACWing import ACMainWing
    
    # Create a wing and changing some parameters
    wing = ACMainWing(1)

#    wing.Lift_LO       = 27.00654 * LBF
    wing.V_max_climb   = 65 * FT/SEC
    wing.V_Stall       = 40 * FT/SEC
#    wing.V_Stall_Init   = 39.2 * FT/SEC
    wing.Alt_LO        = 600 * FT

    wing.S = 900 * IN**2
#    wing.AR = 16
    wing.b = 10. * FT
    wing.TR = [1, 1, 0.7]
    wing.Fb = [0.3, 0.8, 1]
    wing.Gam = [0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG]
    wing.CEdge = 'LE'
    wing.ConstUpper = True
    wing.Airfoil = 'S1223'
    wing.FullWing = True
    wing.ClSlopeAt = (0 * ARCDEG, 7 * ARCDEG)
    wing.CmSlopeAt = (0 * ARCDEG, 8 * ARCDEG)
#    wing.color = 'r'
    
#    wing.Draw()
    
    LLT = ACLiftingLine(wing)
    
    y = npy.linspace(-wing.b/2/IN, wing.b/2/IN, 61)*IN
    
    #LLT.nSpan = 11
    LLT.Alpha3d = 5*ARCDEG
    print LLT.CL()
    print LLT.CDi()
    print LLT.CDv()
    print LLT.CD()
    print LLT.CL()/LLT.CD()
    print LLT.OswaldEff()
    print 'Wing CD : ', wing.CD(5*ARCDEG)
    print AsUnit( wing.Chord(0*IN), 'in' )
    print AsUnit( wing.Chord(18*IN), 'in' )
    print AsUnit( wing.Chord(48*IN), 'in' )
    print AsUnit( wing.Chord(60*IN), 'in' )
    
    pyl.subplot(121)
    pyl.plot(y/IN, LLT.Alpha_i(y)/ARCDEG)
    pyl.xlabel("y (in)")
    pyl.ylabel("Induced Angle (deg)")

    pyl.subplot(122)
    pyl.plot(y/IN, LLT.Cl(y))
    pyl.xlabel("y (in)")
    pyl.ylabel("Cl")
    
    pyl.show()

