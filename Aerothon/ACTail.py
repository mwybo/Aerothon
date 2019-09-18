"""
A collection of classes for tail surfaces; vertical, horizontal, and V-tail
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACLiftSurf import ACLiftSurf
from ACBase import Unitless, Length
from scalar.units import IN, LBF, FT, SEC, ARCDEG, RAD
from scalar.units import AsUnit
import AeroUtil
import numpy as npy
import pylab as pyl
import cmath as math
from RootSolvers import SecantSolver
from scipy.optimize import brentq, brent
from Memoize import memoize

################################################################################
class ACTailSurf(ACLiftSurf):
    """
    A class for a tail surface

    Attributes:
        VC     - Volume coefficient
        L      - Length from the quarter chord to the CG of the aircraft
        eta    - Ratio of tail dynamic pressure to freestream dynamic pressure
    """

#===============================================================================
    def __init__(self, mainwing, index):
        """
        Inputs:
            mainwing  - The main wing used to calculate properties of the tail surface
            index     - Index of the tail surface for AVL
        """
        super(ACTailSurf, self).__init__(index)
        UnitList = self.UnitList

        self.NoneList['VC']  = None     ;   UnitList['VC']  = Unitless
        self.NoneList['L']   = None     ;   UnitList['L']   = Length
        self.__dict__['eta'] = 0.9      ;   UnitList['eta'] = Unitless

        self.SetWing(mainwing)
        self.param.AircraftXcg = 0*IN # The CG of the aircraft. It's up to the aircraft to update this

#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACTailSurf, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'VC':
            ans = self._CalcVC()
        elif key == 'L':
            ans = self._CalcL()
        return ans

#===============================================================================
    def _CalcArea(self):
        " Calculates the tail surface area. "

        VC   = self.VC
        L    = self.L
        Sw   = self.param.mainwing.S

        #
        # Choose the appropriate characteristic length
        # depending on if this is a vertical or horizontal surface
        #
        if self.Axis[0] == 1:
            Lw = self.param.mainwing.MAC()
        else:
            Lw = self.param.mainwing.b/2

        return VC * Lw * Sw/L

#===============================================================================
    def _CalcL(self):
        " Calculates the distance tail surface quarter chord is from the CG of the aircraft. "

        VC   = self.VC
        Sw   = self.param.mainwing.S
        S    = self.S

        #
        # Choose the appropriate characteristic length
        # depending on if this is a vertical or horizontal surface
        #
        if self.Axis[0] == 1:
            Lw = self.param.mainwing.MAC()
        else:
            Lw = self.param.mainwing.b/2

        return VC * Lw * Sw/S

#===============================================================================
    def _CalcVC(self):
        """
        Calculates the tail volume coefficient. 
        """

        L    = self.L                # Length from the aircraft Xcg to the quarter chord of the tail
        Sw   = self.param.mainwing.S # Wing area
        St   = self.S                # Tail area

        #
        # Choose the appropriate characteristic length
        # depending on if this is a vertical or horizontal surface
        #
        if self.Axis[0] == 1:
            Lw = self.param.mainwing.MAC()
        else:
            Lw = self.param.mainwing.b/2

        return L*St/(Lw * Sw)

#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.
        """
        
        #
        # Because special handling is required when computing S, don't call the base
        # class to check the S,b,AR equation
        #
        # super(ACMainWing, self)._CheckConsistent()
        
        NoneList = self.NoneList
        
        if NoneList.has_key('S'):
            #
            # S is calculated so we need three other variables to satisfy both equations for S
            # However, S should be calculated from Lift_LO and V_Stall for a main wing
            #
            self._CheckEquation(['VC','L','AR','b'], Need = 3, message = ' to compute S')
            self._CheckEquation(['VC','L'], Need = 2)
        else:
            #
            # Since we have S, just have to check that the two equations for S are valid
            #
            self._CheckEquation(['S','VC','L'])
            self._CheckEquation(['S','AR','b'])

#===============================================================================
    def SetAircraftXcg(self, AircraftXcg):
        " Sets the CG or the aircraft. "
        self.param.AircraftXcg = AircraftXcg
        #
        # Update the location before calling the base class
        #
        self.X[0] = self.L + self.param.AircraftXcg
        self.dirty = True
        
#===============================================================================
    def RefLen(self, y = None):
        """
        Returns the reference length used for calculating Reynolds number
        
        Inputs:
             y - Spanwise location           
        """
        
        return self.MAC() if y is None else self.Chord(y)
        
#===============================================================================
    def Re(self, V = None, y = None):
        " Calculates the Reynolds number for the wing. "
        Alt = self.param.mainwing.Alt_LO
        
        if V is None:
            V   = self.param.mainwing.GetV_LO()
        
        V  *= self.eta
        rho = AeroUtil.DensAlt(Alt)
        mu  = AeroUtil.Mu(AeroUtil.TempAlt(Alt))
        c   = self.RefLen(y)

        return rho * V * c / mu
    
#===============================================================================
    def SetWing(self, Wing):
        " Set's the wing associated with this tail surface "
        
        self.param.mainwing = Wing
        if Wing is not None:
            Wing.Parents.append( self )
 
#===============================================================================
    def Alpha2D(self, alpha3d, del_c = 0*ARCDEG, Re = None):
        """
        Calculates the 2D sideslip angle given a 3D angle of attack
        Inputs:
            alpha3d - 3D angle of attack
            del_c  - Control surface deflection
            Re     - Reynolds number
        """
        if self.dirty:
            self.Refresh()

        if Re is None:
            Re = self.Re()
            
        def y(alpha2d):
            alpha2d *= ARCDEG
            return (alpha3d - self.DownWash(alpha2d, del_c, Re) + self.i - alpha2d) / ARCDEG

        x0 = alpha3d / ARCDEG
        x1 = x0 + 0.1
         
        #return SecantSolver(x0 = x0, x1 = x1, y = y, tol = 1e-5, itmax = 100)*ARCDEG
        return brentq(f = y, a = -30, b = 30, rtol = 1.e-5, maxiter = 500)*ARCDEG

#===============================================================================
    def CL(self, alpha2d, del_c, Re = None):
        """
        Calculates CL accounting for control surfaces

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
            del_c    - Control surface deflection
        """
        if self.dirty:  self.Refresh()

        if Re is None:
            Re = self.Re()

        CL = ACLiftSurf.CL(self, alpha2d, Re)

        for control in self.Controls.__dict__.itervalues():
            CL += control.DCl(del_c)

        return CL

#===============================================================================
    def CD(self, alpha2d, del_c, Re = None):
        """
        Calculates CD accounting for control surfaces

        Inputs:
            alhpah2d - Angle of attack the wing
            del_c    - Control surface deflection
            Re       - Reynolds number
        """
        if self.dirty:  self.Refresh()

        if Re is None:
            Re = self.Re()

        S = self.S
        CD = self.af.Cd(alpha2d, Re) + self.CDi(alpha2d, del_c, Re)
        CDWinglet = self._WingletCD(self.Winglets, alpha2d) #Computes the winglet drag times winglet area
        CD += CDWinglet/S

        for control in self.Controls.__dict__.itervalues():
            CD += control.DCd(del_c)

        return CD
    
#===============================================================================
    def CDi(self, alpha2d, del_c, Re = None):
        """
        Calculates vortex induced drag of the tail surface

        Inputs:
            alhpah2d - Angle of attack the tail
            del_c    - Control surface deflection
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()
        
        #
        # The lifting surface calls self.af.Cl due to the presence of winglets.
        # however, we have to add the lift due to the control surface deflections
        # here.
        #
        CL = self.af.Cl(alpha2d, Re) * self.FWCF
                
        for control in self.Controls.__dict__.itervalues():
            CL += control.DCl(del_c)
        
        return CL**2/(math.pi * self.o_eff * self.AR)

#===============================================================================
    def CM(self, alpha2d, del_c, Re = None):
        """
        Calculates CM accounting for control surfaces

        Inputs:
            alhpah2d - Angle of attack the wing
            del_c    - Control surface deflection
            Re       - Reynolds number
        """
        if self.dirty:  self.Refresh()

        if Re is None:
            Re = self.Re()

        CM = ACLiftSurf.CM(self, alpha2d, Re)

        for control in self.Controls.__dict__.itervalues():
            CM += control.DCm(del_c)

        return CM
    
#===============================================================================
    def GetV_LO(self):
        """
        Retrieves the lift off velocity for the tail surface.
        """
        return self.param.mainwing.GetV_LO()
    
#===============================================================================
    def GetAlt_LO(self):
        """
        Retrieves the lift off altitude for the tail surface.
        """
        return self.param.mainwing.GetAlt_LO()

#===============================================================================
    def AlphaWing(self, alpha2d, del_c = 0*ARCDEG, Re = None):
        """
        Computes the three dimensional wing angle of attack for plotting purposes
        
        Inputs:
            alhpah2d - 2D angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        return alpha2d + self.DownWash(alpha2d, del_c, Re)

#===============================================================================
    def DownWash(self, alpha2dw, del_c, Re = None):
        """
        Calculates the shed vortex induced downwash angle assuming an elliptical loading.

        Ref: Abbott and Doenhoff pg 7, Prandtl wing theory

        Inputs:
            alhpah2dw - Angle of attack the wing
            del_c     - Control surface deflection
            Re        - Reynolds number
        """
        if self.dirty:  self.Refresh()

        if Re is None:
            Re = self.Re();
        #
        # Calling self.CL here will case an infinite loop for the Alpha2DHTail calculation
        # in the derived ACHorizontalTail class
        #
        return ACTailSurf.CL(self, alpha2dw, del_c, Re)/(math.pi * self.o_eff * self.AR) * RAD


#===============================================================================
    def Refresh(self):
        """
        Updates all the internal quantities.
        """
        self.param.refreshing = True

        #
        # Update the location before calling the base class
        #
        self.X[0] = self.L + self.param.AircraftXcg

        #
        # Now call the base class
        #
        super(ACTailSurf, self).Refresh()
        self.param.refreshing = True

        mainwing = self.param.mainwing
        #
        # Complete the refresh
        #
        V   = mainwing.GetV_LO() * self.eta
        Alt = mainwing.Alt_LO
        mac = self.MAC()

        Re = AeroUtil.Re(Alt, V, mac)
        self.param.AlphaZeroLift = self.af.AlphaZeroLift(Re)


        self.param.refreshing = False

################################################################################
class ACHorizontalTailSurf(ACTailSurf):
    """
    A class for a horizontal tail surface. This is required to account
    for the down wash of the main wing.

    Attributes:
        DWF - Main wing downwash factor (1.0 - 2.0)
    """

#===============================================================================
    def __init__(self, mainwing, index):
        """
        Inputs:
            mainwing  - The main wing used to calculate properties of the horizontal tail surface
            index     - Index of the horizontal tail surface for AVL
        """
        ACTailSurf.__init__(self, mainwing, index)
        UnitList = self.UnitList

        self.__dict__['DWF'] = 2.0  ;   UnitList['DWF'] = Unitless

#===============================================================================
    def _GetAlpha2DH(self, alpha2dw, del_e, Re, iht):
        """
        Calculates possible arrays of 2D horizontal angles of attack

        Inputs:
            alhpah2dw - Angle of attack the wing
            del_e     - Elevator deflection
            Re        - Reynolds number
            i         - Incidence angle
        """
        try:
            if len(alpha2dw) == len(del_e):
                alpha2dh = npy.empty(len(del_e))
                for j in range(len(del_e)):
                    alpha2dh[j] = self.Alpha2DHTail(alpha2dw[j], del_e[j], Re = Re, iht = iht) / ARCDEG
                alpha2dh *= ARCDEG
        except:
            alpha2dh = self.Alpha2DHTail(alpha2dw, del_e, Re = Re, iht = iht)
        
        return alpha2dh
    
#===============================================================================
    def CL(self, alpha2dw, del_e = 0*ARCDEG, Re = None, iht = None):
        """
        Calculates CL accounting for control surfaces and
        altered AoA due to the wing down wash

        Inputs:
            alhpah2dw - Angle of attack the wing
            del_e     - Elevator deflection
            Re        - Reynolds number of the tail surface
            i         - Incidence angle
        """
        if self.dirty:  self.Refresh()
        
        alpha2dh = self._GetAlpha2DH(alpha2dw, del_e, Re, iht)

        return ACTailSurf.CL(self, alpha2dh, del_e, Re)

#===============================================================================
    def CD(self, alpha2dw, del_e = 0*ARCDEG, Re = None, iht = None):
        """
        Calculates CD accounting for control surfaces and
        altered AoA due to the wing down wash

        Inputs:
            alhpah2dw - Angle of attack the wing
            del_e     - Elevator deflection
            Re        - Reynolds number
            iht       - Incidence angle
        """
        if self.dirty:  self.Refresh()

        alpha2dh = self._GetAlpha2DH(alpha2dw, del_e, Re, iht)

        return ACTailSurf.CD(self, alpha2dh, del_e, Re)

#===============================================================================
    def CDi(self, alpha2dw, del_e = 0*ARCDEG, Re = None, iht = None):
        """
        Calculates CDi accounting for control surfaces and
        altered AoA due to the wing down wash

        Inputs:
            alhpah2dw - Angle of attack the wing
            del_e     - Elevator deflection
            Re        - Reynolds number
            iht       - Incidence angle
        """
        if self.dirty:  self.Refresh()

        alpha2dh = self._GetAlpha2DH(alpha2dw, del_e, Re, iht)

        whtw = self._MainWingDownWash(alpha2dw)
        whth = ACTailSurf.DownWash(self, alpha2dh, del_e, Re)

        CL = ACTailSurf.CL(self, alpha2dh, del_e, Re)
        
        return CL*( whth + whtw )

#===============================================================================
    def CM(self, alpha2dw, del_e = 0*ARCDEG, Re = None, iht = None):
        """
        Calculates CM accounting for control surfaces and
        altered AoA due to the wing down wash

        Inputs:
            alhpah2dw - Angle of attack the wing
            del_e     - Elevator deflection
            Re        - Reynolds number
            iht       - Incidence angle
        """
        if self.dirty:  self.Refresh()

        alpha2dh = self._GetAlpha2DH(alpha2dw, del_e, Re, iht)

        return ACTailSurf.CM(self, alpha2dh, del_e, Re)

#===============================================================================
    def _MainWingDownWash(self, alpha2dw):
        """
        Calculates the down wash fom the wing

        Inputs:
            alhpah2dw - Angle of attack the wing
            del_e     - Elevator deflection
            Re        - Reynolds number
            iht       - Incidence angle
        """
        if self.dirty:  self.Refresh()

        return self.DWF * self.param.mainwing.DownWash(alpha2dw)

#===============================================================================
#    @memoize
    def Alpha2DHTail(self, alpha2dw, del_e, Re = None, iht = None):
        """
        Calculates the 2D angle of attack that the tail surface

        Solves the equation

        whtw = DWF*DWW(a2dw)
        a2dht = afus(a2dw) + iht - wht(del_c, a2dht) - whtw

        Where:
            DWF  - Down wash factor (1.0 - 2.0)
            DWW  - Down wash of the wing
            afis - Wing angle of attack relative to the fuselage
            wht  - Down wash of the horizontal tail

        Inputs:
            alhpah2dw    - Angle of attack the wing
            del_e        - Elevator deflection
            initalpha2dh - Initial guess for the tail surface angle of attack
            Re           - Reynolds number of the tail surface
            iht          - Incidence angle of the horizontal tail
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re   = self.Re()
            
        if iht is None:
            iht  = self.i

        try:
            nalpha = len(alpha2dw)            
        except:
            nalpha = 1

        alpha2dh = npy.empty(nalpha)

        #
        # Solve the equation using a secant method
        #
        alpha2dw = self.MakeList(alpha2dw)

        #last = alpha2dw[0] / ARCDEG
        for i in range(nalpha):
            afus = self.param.mainwing.AlphaFus(alpha2dw[i])
            whtw = self._MainWingDownWash(alpha2dw[i])

            # alpha2d == afus + iht - HTail.DownWash - whtw
            def y(alpha2d):
                alpha2d *= ARCDEG 
                return (afus + iht - ACTailSurf.DownWash(self, alpha2d, del_e, Re) - whtw - alpha2d) / ARCDEG

#            x0 = alpha2dw[i]
#            x1 = alpha2dw[i] - 1. * ARCDEG
#            last = SecantSolver(x0 = x0, x1 = x1, y = y, tol = 1.e-5, itmax = 100)
            try:
                last = brentq(f = y, a = -30, b = 30, rtol = 1.e-5, maxiter = 100)
            except ValueError:
                raise ValueError( "Failed to compute horizontal tail angle of attack. This is likely due to a HTail that is too small." )
            
            alpha2dh[i] = last * ARCDEG

        return alpha2dh if len(alpha2dh) > 1 else alpha2dh[0]
    
    
################################################################################
if __name__ == '__main__':
    
    from ACWing import ACMainWing
    wing = ACMainWing(1)

    wing.Lift_LO       = 26.799 * LBF
    wing.V_max_climb   = 65 * FT/SEC
    wing.V_Stall       = 41.27 * FT/SEC
#    wing.V_Stall_Init       = 39.2 * FT/SEC
    wing.Alt_LO        = 600 * FT

#    wing.S = 900 * IN**2
    wing.AR = 16
#    wing.b = 10. * FT
    wing.TR = [1, 0.7, 0.7]
    wing.Fb = [0.3, 0.8, 1]
    wing.Gam = [0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG]
    wing.CEdge = 'LE'
    wing.ConstUpper = True
    wing.Airfoil = 'S1223'
    wing.Symmetric = True
    wing.ClSlopeAt = (3 * ARCDEG, 4 * ARCDEG) 
#    wing.color = 'r'
    wing.LLFile = 'LLTPolar.txt'

    ail = wing.AddControl('ail')


#------------------------------------------------------------------------------    
    htail = ACHorizontalTailSurf(wing,2)
    htail.Airfoil = 'NACA0012'
#    htail.AR = 3
    htail.TR = 0.5
    htail.Gam = 0 * ARCDEG
    htail.Lam = 0 * ARCDEG
    htail.S   = 100 * IN**2
    htail.b   = 20 *IN
#    htail.L        = 30 * IN
    htail.VC = 0.85
    htail.SweepFc = 0.75 #Set the sweep to zero about the elevator hinge
    htail.DWF = 2.
    htail.i = 0. * ARCDEG
    htail.ClSlopeAt = None
    htail.FullWing = True
    
    htail.X = [10*IN, 0*IN, 10*IN]

    elev = htail.AddControl('elev')
    elev.Fb = 1.
    elev.Fc = 0.25
    elev.Ft = 0.

    vtail = ACTailSurf(wing, 3)
    vtail.Airfoil = 'NACA0012'
    vtail.VC = 0.06
    vtail.AR = 1.6
    vtail.TR = 0.7
    vtail.Symmetric = True
    vtail.Axis = (0, 1)
    vtail.L = 51.572 * IN
    vtail.SweepFc = 0.6 #Set the sweep to zero about the rudder hinge

    vtail.FullWing = True
    vtail.X = [10*IN, 10*IN, 10*IN]
    vtail.Symmetric = True
    
    rud = vtail.AddControl('rud')
    rud.Fb = 1.
    rud.Fc = 0.4
    rud.Ft = 0.
    
    
    print AsUnit( htail.S, "m**2")
    print 'Htail span: ', htail.b, 'area', AsUnit( htail.S, "in**2"), 'rChord', htail.Chord(0*IN)
    print 'VTail span: ', vtail.b, 'area', AsUnit( vtail.S, "in**2"), 'rChord', vtail.Chord(0*IN)

    del_c = 0 * ARCDEG
    alpha2dw = 5 * ARCDEG
    print 'Horizontal Tail angle', htail.Alpha2DHTail(alpha2dw, del_c, iht = 0.*ARCDEG)
    print 'Horizontal Tail Cl', htail.CL(alpha2dw, del_c)
    print 'Horizontal Tail DW', AsUnit( htail.DownWash(alpha2dw, del_c), "deg")
    print 'Horizontal Tail Re', htail.Re()
    
    
    htail.Draw2DAirfoilPolars(2)
    htail.Draw()
    vtail.Draw()
    pyl.show()
    