"""
Class for wing surface.

Key implementation is the area calculation. Here it is implemented as

S = MaxLift/(q_LO * CL_LO)
where
MaxLift   - The desired maximum lift at lift off
q         - Dynamic pressure at lift off
CL_LO     - CL at lift off
CLmax     - Maximum lift coefficient of the wing
V_LOstall - Min velocity factor to prevent stall (typically 1.2)

CL_LO = CLmax/V_LOstall**2

"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACLiftSurf import ACLiftSurf
from ACBase import ACBase, ACComponent, ACStructural, Unitless, Velocity, Force, Area, Length, Angle, g
from ACAVLWing import ACAVLWing
from scalar.units import FT, SEC, LBF, IN, ARCDEG, RAD, LBM
from scalar.units import AsUnit
import AeroUtil
import numpy as npy
import pylab as pyl
import cmath as math
from scipy.interpolate import RectBivariateSpline, UnivariateSpline, interp1d
from RootSolvers import SecantSolver
from scipy.optimize import brentq

################################################################################
class ACWingError(Exception):
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg

################################################################################
class ACMainWing(ACLiftSurf):
    """
    A class for a main wing

    Attributes
        V_max_climb  - Maximum climb velocity
        V_Stall      - Stall velocity
        Lift_LO      - Desired lift from the wing at lift of
        Alt_LO       - Lift off altitude
        V_LOstall    - Min velocity factor to prevent stall (typically 1.2)
        S_Reinit     - Estimated initial reynolds number for wing area calculation
        S_IterMax    - Maximum number of iterations to calculate wing area
        S_Init       - Initial guess for the area
        V_Stall_Init - Initial guess for the V_Stall
        
        DrawDetail   - If true, the elliptic chord distribution and quarter chord are drawn
    """

#===============================================================================
    def __init__(self, index):
        """
        Inputs:
            index - Index of the main wing for AVL
        """
        super(ACMainWing, self).__init__(index)
        UnitList = self.UnitList

        self.NoneList['V_Stall']     = None             ;   UnitList['V_Stall']     = Velocity
        self.NoneList['Lift_LO']     = None             ;   UnitList['Lift_LO']     = Force
        
        self.__dict__['V_max_climb'] = 0 * FT/SEC       ;   UnitList['V_max_climb'] = Velocity
        self.__dict__['Alt_LO']      = 1129 * FT        ;   UnitList['Alt_LO']      = Length
        self.__dict__['V_LOstall']   = 1.2              ;   UnitList['V_LOstall']   = Unitless
        
        self.__dict__['V_Stall_Init']= 40 * FT/SEC      ;   UnitList['V_Stall_Init']= Velocity
        self.__dict__['S_ReInit']    = 4e5              ;   UnitList['S_ReInit']    = Unitless
        self.__dict__['S_IterMax']   = 50               ;   UnitList['S_IterMax']   = Unitless
        self.__dict__['S_Init']      = 1700 * IN**2     ;   UnitList['S_Init']      = Area
        self.__dict__['DrawDetail']  = True
        
        #
        # Create the drawings for the elliptic chord distribution and quarter chord
        #
        self.Drawings.Ellips = ACComponent.ComponentDrawing()
        self.Drawings.Ellips.style = '--'

        self.Drawings.QuarterChord = ACComponent.ComponentDrawing()
        self.Drawings.QuarterChord.color = 'm'
        self.Drawings.QuarterChord.style = '--'

        #
        # Initialize internal variables
        #
        self.param.S = None
        self.param.VStall = None
        
        self.param.ComputingS = False #A flag to track if S is currently being computed

#===============================================================================
    def _UpdateArea(self):
        " Calculates the wing area. "

        #
        # If the wing is clean don't recalculate the area
        #
        if not self.NoneList.has_key('S'):
            self.param.S = self.S
            return
        
        if not self.NoneList.has_key('AR') and not self.NoneList.has_key('b'):
            self.param.S = super(ACMainWing,self)._CalcArea()
            return
        
        Lift_LO   = self.Lift_LO
        Alt       = self.Alt_LO
        V         = self.GetV_LO()
        rho       = AeroUtil.DensAlt(Alt)
        q_LO      = 0.5 * rho * V**2
        V_LOstall = self.V_LOstall
        FWCF      = self.FWCF
        Re        = self.S_ReInit
        mu        = AeroUtil.Mu(AeroUtil.TempAlt(Alt))
        self.__dict__['S'] = self.S_Init

        # Get an it
        CLmax   = self.CLmax(Re)
        CL_LO   = CLmax/V_LOstall**2
        self.S  = Lift_LO/(q_LO * CL_LO)

        #
        # Set that S is being calculated for the consistency check
        #
        self.__dict__['param'].ComputingS = True #To make sure self does not get dirty
        #
        # Because the area is dependent on Clmax which is dependent
        # on Reynolds number which is dependent on the MAC which
        # is dependent on the area, the area must be calculated
        # iteratively
        #
        it = 0
        Cinit = 0
        tol = 1e-2
        resid = 1 + tol
        S_old = self.S_Init
        while it < self.S_IterMax and abs(resid) > tol:
            it     += 1
            #
            # To improve performance, a hack is done here to limit the
            # update to only what is required. This is needed for the MAC equation.
            #
            self._CalcChordSpanLoc()  # Force an update of chord locations
            self._CalcChord()         # Force an update of chords
            self.param.MAC = None     # Force an update to the MAC
            self.__dict__['dirty'] = False

            mac     = self.MAC()
            Re      = rho * V * mac / mu
            
            CLmax   = self.CLmax(Re)
            CL_LO   = CLmax/V_LOstall**2
            S       = Lift_LO/(q_LO * CL_LO)
            self.S  = S

            resid   = S / (IN**2) - S_old / (IN**2)
            S_old   = S

#            print self.name, ':  CLmax = ', CLmax, 'Re = ',Re, 'S = ', AsUnit( S, 'in' ), 'resid =', abs(resid)

        self.param.S = self.S
        self.S = None
        self.S_Init = self.S
        self.S_ReInit = Re
        self.param.ComputingS = False #S has been computed

#===============================================================================
    def _CalcArea(self):
        """ Returns the computed area """
        if self.dirty: self.Refresh()
        return self.param.S
 
#===============================================================================
    def CL_max_climb(self, Re):
        """
        Calculates the lift coefficient for max climb

        Inputs:
            Re - Reynolds number
        """

        V_max_stall = self.V_max_climb / self.V_Stall
        return self.CLmax(Re) / V_max_stall**2

#===============================================================================
    def RefLen(self, y = None):
        """
        Returns the reference length used for calculating Reynolds number
        
        Inputs:
             y - Spanwise location           
        """
        
        return self.MAC() if y is None else self.Chord(y)
    
#===============================================================================
    def Re(self, y = None, V = None):
        """
        Calculates the Reynolds number for the wing.

        This is a relatively cheap calculation and is therefore left
        as a calculation rather than storing a self.param.Re
        to reduce the interdependencies.
        
        Inputs:
            y - Spanwise location
            V - Freestrea velocity
        """
        if V is None:
            V   = self.GetV_LO()

        Alt = self.Alt_LO
        rho = AeroUtil.DensAlt(Alt)
        mu  = AeroUtil.Mu(AeroUtil.TempAlt(Alt))
        c   = self.RefLen(y)

        return rho * V * c / mu

#===============================================================================
    def AlphaFus(self, alpha2d, Re = None):
        """
        Calculates the angle of attack relative to the fuselage for plotting purposes

        Inputs:
            alhpah2d - Angle of attack of the wing
            Re       - Reynolds number
        """
        if self.dirty:
            self.Refresh()

        return self.AlphaWing(alpha2d, Re) - self.i        

#===============================================================================
    def Alpha2D(self, AlphaFus, Re = None):
        """
        Calculates the 2D angle of attack given a fuselage angle of attack

        Inputs:
            AlphaFus - Angle of attack of the fuselage
            Re       - Reynolds number
        """
        if self.dirty:
            self.Refresh()

        if Re is None:
            Re = self.Re()
            
        def y(alpha2d):
            alpha2d *= ARCDEG
            return (a - self.DownWash(alpha2d, Re) + self.i - alpha2d) / (ARCDEG)

        AlphaFus = self.MakeList(AlphaFus)

        a2d = []
        for a in AlphaFus:
            x0 = a / (ARCDEG)
            x1 = x0 + 0.1
#            a2d.append(SecantSolver(x0 = x0, x1 = x1, y = y, tol = 1e-5, itmax = 100))
            a2d.append(brentq(f = y, a = -30, b = 30, rtol = 1.e-5, maxiter = 500))
        
        return a2d[0]*ARCDEG if len(a2d) == 1 else npy.array(a2d)*ARCDEG
    
#===============================================================================
    def _UpdateAlpha2d_LO(self):
        """
        Calculates the two dimensional angle of attack at Lift off
        """
        Re = self.Re()
        CL_LO = self.GetCL_LO()
        
        AlphaClmin = self.af.AlphaClmin(Re)/ARCDEG
        AlphaClmax = self.af.AlphaClmax(Re)/ARCDEG
        
        def y(alpha2d):
            return self.CL(alpha2d*ARCDEG, Re) - CL_LO
        
#        self.param.Alpha2d_LO = SecantSolver(x0 = 0, x1 = 1, y = y, tol = 1e-5, itmax = 100)*ARCDEG
        self.param.Alpha2d_LO = brentq(f = y, a = AlphaClmin, b = AlphaClmax, rtol = 1.e-5, maxiter = 500)*ARCDEG
        
#===============================================================================
    def _UpdateVStall(self):
        """
        Determines Velocity at which wing stalls
        """
        #
        # If the wing VStall then don't calculate it
        #
        if not self.NoneList.has_key('V_Stall'):
            self.param.VStall = self.V_Stall
            return

        Lift_LO   = self.Lift_LO
        Alt       = self.Alt_LO
        V_Stall   = self.V_Stall_Init
        rho       = AeroUtil.DensAlt(Alt)
        V_LOstall = self.V_LOstall
        S         = self.S
        mu        = AeroUtil.Mu(AeroUtil.TempAlt(Alt))

        #
        # To improve performance, a hack is done here to limit the
        # update to only what is required. This is needed for the MAC equation.
        #
        self._CalcChordSpanLoc()  # Force an update of chord locations
        self._CalcChord()         # Force an update of chords
        self.__dict__['dirty'] = False

        mac       = self.MAC()

        #
        # Because VStall is dependent on Clmax which is dependent
        # on Reynolds number which is dependent on the V_LO which
        # is dependent on VStall, the VStall must be calculated
        # iteratively
        #
        it = 0
        Cinit = 0
        tol = 1e-3
        resid = 1 + tol
        V_Stall_old = 0 * FT/SEC
        while it < self.S_IterMax and abs(resid) > tol:
            it     += 1
            V = V_Stall * V_LOstall
            Re = rho * V * mac / mu
            CLmax   = self.CLmax(Re)
            
            V_Stall = ((2*Lift_LO)/(rho*CLmax*S))**0.5
            
            resid = (V_Stall - V_Stall_old) / (FT/SEC)
            V_Stall_old = V_Stall
            
            #print 'CLmax = ', CLmax, 'Re = ',Re, 'V_Stall = ', AsUnit( V_Stall, 'ft/s' )

        self.param.VStall = V_Stall

        
#===============================================================================
    def _CalcVStall(self):
        """ Returns the computed V_Stall """
        if self.dirty: self.Refresh()
        return self.param.VStall
        
#===============================================================================
    def _CalcLift_LO(self):
        """
        Calculates the lift at lift off
        """
        S         = self.S
        Alt_LO    = self.Alt_LO
        V_LO      = self.GetV_LO()

        CL_LO     = self.GetCL_LO()
        rho       = AeroUtil.DensAlt(Alt_LO)
        q         = 0.5*rho*V_LO**2

        # Just multiply the CL at lift off with q and area
        Lift_LO = q * CL_LO * S

        return Lift_LO
    
#===============================================================================
    def GetCL_LO(self):
        """
        Returns the lift off CL.
        """
        CLmax     = self.CLmax()
        V_LOstall = self.V_LOstall

        return CLmax/V_LOstall**2
    
#===============================================================================
    def GetV_LO(self):
        """
        Determines Velocity required for Lift-Off
        """
        if self.dirty: self.Refresh()
        
        V_Stall = self.V_Stall
        V_LOstall = self.V_LOstall
        
        VLO = V_LOstall * V_Stall

        return VLO

#===============================================================================
    def GetAlt_LO(self):
        """
        Returns the altitude for Lift-Off. This method is needed so winglets work properly
        """
        if self.dirty: self.Refresh()

        return self.Alt_LO

#===============================================================================
    def GetAlpha2d_LO(self):
        """
        Returns the angle of attack at Lift-Off.
        """
        if self.dirty: self.Refresh()

        return self.param.Alpha2d_LO

#===============================================================================
    def GetAlpha3d_LO(self):
        """
        Returns the 3D angle of attack at Lift-Off.
        """
        if self.dirty: self.Refresh()

        return self.AlphaWing(self.param.Alpha2d_LO)
    
#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACMainWing,self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'Lift_LO':
            return self._CalcLift_LO()
        elif key == 'V_Stall':
            return self._CalcVStall()
        #
        # If no calculation exist return None
        #
        return None
    

#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.
        """
        super(ACMainWing, self)._CheckConsistent()
        
        NoneList = self.NoneList
        
        if self.param.ComputingS or NoneList.has_key('S'):
            #
            # S is calculated so we need three other variables to satisfy both equations for S
            # However, S should be calculated from Lift_LO and V_Stall for a main wing
            #
            self._CheckEquation(['Lift_LO','V_Stall','AR','b'], Need = 3, message = ' to compute S')
        else:
            #
            # Since we have S, just have to check that the two equations for S are valid
            #
            self._CheckEquation(['S','Lift_LO','V_Stall'])
            self._CheckEquation(['S','AR','b'])


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

        self.__dict__['dirty'] = False
        self.param.refreshing = True
        #
        # Due to the cyclic behavior of the area calculation, it must be updated before anything else.
        #
        self._CheckConsistent()
        self._CalcAirfoils() #The airfoil is also needed for the area
        self._UpdateArea()
        self._UpdateVStall()
                
        #
        # Now update the base class (This will also call UpdateDrawing)
        #
        super(ACMainWing, self).Refresh()
        
        #
        # Now that everything is up to date, update the Alpha2d calculation
        #
        self._UpdateAlpha2d_LO()

        self.param.refreshing = False
#===============================================================================
    def UpdateDrawing(self):
        """
        Updates all drawings of the wing
        """
        super(ACMainWing, self).UpdateDrawing()
        
        if self.DrawDetail:
            #
            # Also draw an elliptic chord distribution
            #
            Ellips = self.Drawings.Ellips
    
            b = self.b/2. if self.FullWing else self.b
            
            npel = 31 # Number of points to draw the ellips
            x = npy.empty(npel)
            y = npy.empty(npel)
    
            for i in range(npel):
                yl = b * float(i)/float(npel-1)
                xl = self.LE(yl) + self.EllipChord(yl)
                y[i] = (yl + self.X[1]) / IN
                x[i] = xl / IN
    
            Ellips.xt = x
            Ellips.yt = y
    
            #
            # and a quarter chord line
            #
            QuarterChord = self.Drawings.QuarterChord
            nspan = len(self.param.span)
            x = npy.empty(nspan)
            y = self.param.span
    
            for i in range(nspan):
                yl = y[i] * IN
                xl = self.LE(yl) + 0.25*self.Chord(yl)
                x[i] = xl / IN
    
            QuarterChord.xt = x
            QuarterChord.yt = y + self.X[1] / IN
        

################################################################################
class ACBiWing(ACStructural, ACAVLWing):
    """
    A class for a bi wing
    
    Attributes
        GapInterp  - An array of gap/span values used for the interpolation of Oswald efficiency
        OeffInterp - An array of Oswald efficiencies (Use AVL to find the o_eff for a few gap/span values)
        BWFCInterp - An array of biwing correction factors (use the same gap values as GapInterp
        Gap        - The gap/span value of the lower wing
        Stagger    - Percent Stagger of the LE of top wing to the LE of the bottom wing over average MAC
                          (0.0 is no stagger, + is top wing forward)
        b          - Optional span arguement used for equal span, limited design
        S          - Total area of the two wings
        AR         - The aspect ratio which is the average span squared divided by the total wing area
        Lift_Ratio - Ratio of lift between the two wings
        Lift_LO    - Desired lift from the wing at lift off
        FullWing   - True/False If true a full wing will be generated, otherwise only a half wing
        
        V_max_climb  - Maximum climb velocity
        V_Stall      - Stall velocity
        Alt_LO       - Lift off altitude
        V_LOstall    - Min velocity factor to prevent stall (typically 1.2)
 
        ClSlopeAt   - Angle of attack to calculate the Cl slope of the biwing
        CmSlopeAt   - Angle of attack to calculate the Cm slope of the biwing
       
        AVLexe      - The path to the AVL executable (default 'avl.exe')
    """
    
#===============================================================================
    def __init__(self, index, index2, index3):
        """
        Inputs:
            index  - Index of the Lower wing for AVL
            index2 - Index of the Upper wing for AVL
            index3 - Index of the optional end plate
        """
        super(ACBiWing,self).__init__()
        UnitList = self.UnitList
        
        self.__dict__['LowerWing'] = ACMainWing(index)
        self.__dict__['UpperWing'] = ACMainWing(index2)
        self.NoneList['EndPlate']  = None
        self.param.EndPlateIndex = index3
        
        self.__dict__['GapInterp']   = [0.1,0.2,0.3,0.4]  ;   UnitList['GapInterp']   = Unitless
        self.__dict__['OeffInterp']  = [1.,1.,1.,1.]      ;   UnitList['OeffInterp']  = Unitless
        self.__dict__['BWCFInterp']  = [1.,1.,1.,1.]      ;   UnitList['OeffInterp']  = Unitless
        self.__dict__['Gap']         = 0.1                ;   UnitList['Gap']         = Unitless
        self.__dict__['Stagger']     = 0.0                ;   UnitList['Stagger']     = Unitless
        self.NoneList['b']           = None               ;   UnitList['b']           = Length
        self.NoneList['S']           = None               ;   UnitList['S']           = Area
        self.NoneList['TotalS']      = None               ;   UnitList['TotalS']      = Area
        self.NoneList['AR']          = None               ;   UnitList['AR']          = Unitless
        self.__dict__['FullWing']    = False
 
        self.NoneList['V_Stall']     = None               ;   UnitList['V_Stall']     = Velocity
        self.NoneList['Lift_LO']     = None               ;   UnitList['Lift_LO']     = Force
        self.NoneList['Lift_Ratio']  = None               ;   UnitList['Lift_Ratio']  = Unitless
        
        self.__dict__['V_max_climb'] = 0 * FT/SEC         ;   UnitList['V_max_climb'] = Velocity
        self.__dict__['Alt_LO']      = 600 * FT           ;   UnitList['Alt_LO']      = Length
        self.__dict__['V_LOstall']   = 1.2                ;   UnitList['V_LOstall']   = Unitless
 
        self.NoneList['ClSlopeAt']  = None                ;   UnitList['ClSlopeAt']   = Angle
        self.NoneList['CmSlopeAt']  = None                ;   UnitList['CmSlopeAt']   = Angle
       
        self.__dict__['AVLexe']     = None
        
        #
        # If the lower of upper wing is dirtied, the biwing will be dirtied
        #
        self.LowerWing.Parents.append( self )
        self.UpperWing.Parents.append( self )

        #
        # Give some reasonable names to the wings
        #
        self.LowerWing.name = 'LowerWing'
        self.UpperWing.name = 'UpperWing'
        #
        # This is needed to access control surfaces on the wings
        #
        self.NoneList['Controls'] = None

#===============================================================================
    def CreateEndPlate(self):
        """
        Creates and end plate for the biwing
        """

        self.EndPlate = ACEndPlate(self.param.EndPlateIndex, self)
        #
        # If the endplate is dirtied, the biwing will be dirtied
        #
        self.EndPlate.name = 'End Plate'
        self.EndPlate.Parents.append( self )

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the bi wing and its parts to the Weight Table
        
        Input:
            PartName    - The name give to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WeightTable[PartName] = self
        
        self.UpperWing.AddToWeightTable('UpperWing', WeightTable[PartName])
        self.LowerWing.AddToWeightTable('LowerWing', WeightTable[PartName])
        
        if self.EndPlate is not None:
            self.EndPlate.AddToWeightTable('EndPlate', WeightTable[PartName])
        
#===============================================================================
    def CL(self,alpha2d, Re = None):
        """
        Calculates CL for biplane combination taken to be area average of two wing CLs
        """
        if self.dirty: self.Refresh()
        
        UWS  = self.UpperWing.S
        LWS  = self.LowerWing.S
        UWCL = self.UpperWing.CL(alpha2d, Re)
        LWCL = self.LowerWing.CL(alpha2d, Re)
        BWCF = self.BWCF()
        
        return (UWCL*UWS + LWCL*LWS) * BWCF / (UWS + LWS)

#===============================================================================
    def CL2(self,alpha2d, Re = None):
        """
        Calculates CL for biplane combination using inAsUnition from Altman et. al.
        Taken from AIAA 2009-1086-302, Altman, Kang. Pt II Pg 3
        """
        if self.dirty: self.Refresh()
        
        St = self.Stagger
        G  = self.Gap * self.b / self.MAC()
        
        dCL_da_UW = self.UpperWing.dCL_da()
        dCL_da_LW = self.LowerWing.dCL_da()
        
        dCL_da_predict = (-0.007 * St + 0.015)*G +(0.018 * St + 0.05)
    
        Corr_Fctr      = dCL_da_predict / ((dCL_da_UW + dCL_da_LW)/2) 
        
        UWS  = self.UpperWing.S
        LWS  = self.LowerWing.S
        UWCL = self.UpperWing.CL(alpha2d, Re)
        LWCL = self.LowerWing.CL(alpha2d, Re)
        
        return (UWCL*UWS + LWCL*LWS) * Corr_Fctr/ (UWS + LWS)

#===============================================================================
    def CDi(self, alpha2d, Re = None):
        """
        Calculates induced drag coefficient for biplane as the sum of the two wings
        """
        if self.dirty: self.Refresh()
        
        o_eff = self.O_eff()
        CL    = self.CL(alpha2d, Re)
        AR    = self._CalcAR()
        
        return CL**2/(math.pi * o_eff * AR)
    
#===============================================================================
    def CD(self, alpha2d, Re = None):
        """
        Calculates CD including the induced drag
        
        Inputs:
            alhpah2d - Angle of attack the biwing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()
        
        #
        # Compute an area weighted drag coefficient between the wing
        #
        UW = self.UpperWing
        LW = self.LowerWing
        EP = self.EndPlate
        
        if Re is None:
            Re = LW.Re()
         
        S    = UW.S + LW.S
        CD   = UW.af.Cd(alpha2d, Re) + LW.af.Cd(alpha2d, Re) + self.CDi(alpha2d, Re)
        if self.EndPlate is not None:
            EPS  = EP.S if not EP.Symmetric else 2*EP.S 
            CDEP = EP.af.Cd(0*ARCDEG, Re)*EPS #The angle of attack changes with sideslip for the endplates
        else:
            CDEP = 0
        
        return CD + CDEP / S

#===============================================================================
    def CM(self, alpha2d, Re = None):
        """
        Calculates CM as the sum of the wing CMs
        
        Inputs:
            alpha2d - Two dimensional wing angle of attack
            Re      - Reynolds number
        """
        if self.dirty: self.Refresh()

        UW = self.UpperWing
        LW = self.LowerWing

        #
        # Normalize the CM to the area of the two wings
        #
        return (UW.CM(alpha2d,Re)*UW.S + LW.CM(alpha2d,Re)*LW.S) / (UW.S + LW.S)
 
  #===============================================================================
    def dCL_da(self, alpha2d = None):
        """
        Calculates CL slope of the two wings
        
        Inputs:
            alpha2d - Two dimensional wing angle of attack
        """
        if alpha2d is None:
            alpha2d = self.ClSlopeAt

        return self._da(self.CL, alpha2d)
   
 #===============================================================================
    def dCM_da(self, alpha2d = None):
        """
        Calculates CM slope of the two wings
        
        Inputs:
            alpha2d - Two dimensional wing angle of attack
        """
        if alpha2d is None:
            alpha2d = self.CmSlopeAt
            
        return self._da(self.CM, alpha2d)
   
#===============================================================================
    def _da(self, Cfn, alph2d):
        """
        Calculates the slope of the function Cfn
        
        Inputs:
            Cfn     - The function from which a slope is desired
            alpha2d - The 2D wing angle of attack
        """
        if self.dirty:  self.Refresh()
        
        a1 = alph2d[0]
        a2 = alph2d[1]
        dela = a2 - a1
        #
        # Calculate the Cfn slope
        #
        dCfn_da = (Cfn(a2) - Cfn(a1))/dela

        return dCfn_da

#===============================================================================
    def CLmax(self, Re = None):
        """
        Calculates max CL with the finite wing correction factor

        Inputs:
            Re - Reynolds number
        """
        if self.dirty: self.Refresh()

        UWS  = self.UpperWing.S
        LWS  = self.LowerWing.S
        UWCLmax = self.UpperWing.CLmax()
        LWCLmax = self.LowerWing.CLmax()
        BWCF = self.BWCF()
        
        return (UWCLmax*UWS + LWCLmax*LWS) * BWCF / (UWS + LWS)
                       
#===============================================================================
#    def DownWash(self, alpha2d, Re=None):
#        """
#        Calculate Downwash as the sum of the downwash from both wings
#        """
#        if self.dirty: self.Refresh()
#        
#        UW = self.UpperWing
#        LW = self.LowerWing
#
#        return UW.DownWash(alpha2d, Re) + LW.DownWash(alpha2d, Re)

#===============================================================================
    def DownWash(self, alpha2d, Re = None):
        """
        Calculate Downwash as the sum of the downwash from both wings
        """
        if self.dirty: self.Refresh()
        
        o_eff = self.O_eff()
        CL    = self.CL(alpha2d, Re)
        AR    = self._CalcAR()
        
        return CL/(math.pi * o_eff * AR)*RAD

#===============================================================================
    def Xac(self, Re = None):
        """
        Calculate Aerodynamic Center of biplane wings
        """
        if self.dirty: self.Refresh()
        
        UW = self.UpperWing
        LW = self.LowerWing
        
        #
        # Make sure that the X[0] location is up to date.
        # setting Wing.X[0] = Xw does not dirty the BiWing
        #
        self._UpdateStagger()
                
        return (UW.Xac(Re) + LW.Xac(Re)) / 2
    
#===============================================================================
    def MAC(self):
        """
        Return Average MAC of wings
        """
        if self.dirty: self.Refresh()

        return (self.UpperWing.MAC() + self.LowerWing.MAC()) / 2

#===============================================================================
    def CG(self):
        """
        Computes center or gravity of the biwing
        """
        if self.dirty: self.Refresh()

        UW = self.UpperWing
        LW = self.LowerWing
        EP = self.EndPlate
        
        CG      = UW.CG()*UW.Weight + LW.CG()*LW.Weight
        TWeight = UW.Weight + LW.Weight
        if EP is not None:
            CG      += EP.CG()*EP.Weight
            TWeight += EP.Weight
            
        return CG/TWeight

#===============================================================================
    def MOI(self):
        """
        Computes the moments of inertia of the biwing
        """
        if self.dirty: self.Refresh()
        
        UW = self.UpperWing
        LW = self.LowerWing
        EP = self.EndPlate
        
        CG  = self.CG()
        MOI = npy.array([0,0,0]) * LBM*IN**2

        def AddMOI(PCG, PM, PMOI):
            dCG = PCG - CG
            MOI[0] += PMOI[0] + PM*(dCG[1]**2 + dCG[2]**2) #Ixx
            MOI[1] += PMOI[1] + PM*(dCG[0]**2 + dCG[2]**2) #Iyy
            MOI[2] += PMOI[2] + PM*(dCG[1]**2 + dCG[0]**2) #Izz
        
        AddMOI(UW.CG(), UW.Weight / g, UW.MOI())
        AddMOI(LW.CG(), LW.Weight / g, LW.MOI())
        if EP is not None:
            AddMOI(EP.CG(), EP.Weight / g, EP.MOI())
        
        return MOI

#===============================================================================
    def AlphaWing(self,alpha2d, Re=None):
        """
        Return the three dimensional angle of attack for the finite wing
        """
        if self.dirty: self.Refresh()
        
        return alpha2d + self.DownWash(alpha2d, Re)
        
#===============================================================================
    def AlphaFus(self,alpha2d, Re=None):
        """
        Return Relative angle of attack of wing to fuselage
        """
        if self.dirty: self.Refresh()
        
        return self.AlphaWing(alpha2d, Re) - (self.LowerWing.i + self.UpperWing.i)/2
    
#===============================================================================
    def Alpha2D(self, AlphaFus, Re = None):
        """
        Calculates the 2D angle of attack given a fuselage angle of attack

        Inputs:
            AlphaFus - Angle of attack of the fuselage
            Re       - Reynolds number
        """
        if self.dirty:
            self.Refresh()
            
        def y(alpha2d):
            alpha2d *= ARCDEG
            return (a - self.DownWash(alpha2d, Re) + self.LowerWing.i - alpha2d) / ARCDEG
        
        AlphaFus = self.MakeList(AlphaFus)
        
        a2d = []
        for a in AlphaFus:
            x0 = a / ARCDEG
            x1 = x0 + 0.1
#            a2d.append(SecantSolver(x0 = x0, x1 = x1, y = y, tol = 1e-5, itmax = 100))
            a2d.append( brentq(f = y, a = -30, b = 30, rtol = 1.e-5, maxiter = 500) )
         
        return a2d[0]*ARCDEG if len(a2d) == 1 else npy.array(a2d)*ARCDEG
    
#===============================================================================
    def _UpdateAlpha2d_LO(self):
        """
        Calculates the two dimensional angle of attack at Lift off
        """
        Re = self.LowerWing.Re()
        CL_LO = self.GetCL_LO()
        
        def y(alpha2d):
            return self.CL(alpha2d*ARCDEG, Re) - CL_LO
        
#        self.param.Alpha2d_LO = SecantSolver(x0 = 0, x1 = 1, y = y, tol = 1e-5, itmax = 100)*ARCDEG
        self.param.Alpha2d_LO = brentq(f = y, a = -30, b = 30, rtol = 1.e-5, maxiter = 500)*ARCDEG

#===============================================================================
    def GetAlpha2d_LO(self):
        """
        Returns the angle of attack at Lift-Off.
        """
        if self.dirty: self.Refresh()

        return self.param.Alpha2d_LO

#===============================================================================
    def GetAlpha3d_LO(self):
        """
        Returns the 3D angle of attack at Lift-Off.
        """
        if self.dirty: self.Refresh()

        return self.AlphaWing(self.param.Alpha2d_LO)
    
#===============================================================================
    def Upper(self, y, x = None):
        """
        Returns the upper side if the upper wing given a spanwise location

        Inputs:
            y - span location
            x - chord location (Min y value if None)
                Note that the airfoil max and min y values may
                not correspond to the same x location
        """
        if self.dirty: self.Refresh()
        
        self._UpdateStagger()
        
        return self.UpperWing.Upper(y,x)

#===============================================================================
    def Lower(self, y, x = None):
        """
        Returns the lower side given a spanwise location

        Inputs:
            y - span location
            x - chord location (Min y value if None)
                Note that the airfoil max and min y values may
                not correspond to the same x location
        """
        if self.dirty: self.Refresh()

        self._UpdateStagger()
        
        return self.LowerWing.Lower(y,x)
    
#===============================================================================
    def GetCL_LO(self):
        """
        Returns the lift off CL.
        """
        CLmax     = self.CLmax()
        V_LOstall = self.V_LOstall
        
        return CLmax/V_LOstall**2

#===============================================================================
    def Re(self, y = None):
        """
        Returns average reynolds number of the two wings.
        """
        if self.dirty: self.Refresh()
        
        return (self.LowerWing.Re(y) + self.UpperWing.Re(y))/2

#===============================================================================
    def GetV_LO(self):
        """
        Determines Velocity required for Lift-Off
        """
        if self.dirty: self.Refresh()
        
        V_Stall = self.V_Stall
        V_LOstall = self.V_LOstall
        
        VLO = V_LOstall * V_Stall

        return VLO
    
#===============================================================================
    def GetLift_LO(self):
        """
        Returns the Lift produced at Lift-off conditions
        """
        if self.dirty: self.Refresh()
        
        L_LO = self.UpperWing.Lift_LO + self.LowerWing.Lift_LO

        return L_LO

#===============================================================================
    def GetAlt_LO(self):
        """
        Returns the altitude for Lift-Off. This method is needed so winglets work properly
        """
        if self.dirty: self.Refresh()

        return self.Alt_LO

#===============================================================================
    def AlphaRange(self, nalpha = 30):
        if self.dirty: self.Refresh()
        
        return self.LowerWing.AlphaRange(nalpha)
    
#===============================================================================
#    def _CalcMunksNumber(self):
#        """
#        Calculation of influence coefficient as determined in NACA Technical Note 182
#            Pradtl - Induced Drag of Multiplanes
#            Pages 7 to 12
#        """
#        
#        UW = self.UpperWing
#        LW = self.LowerWing
#        G = self.Gap
#        
#        # b1 is larger span, b2 is smaller span
#        if (UW.b >= LW.b):
#            b1 = UW.b
#            b2 = LW.b
#        else:
#            b1 = LW.b
#            b2 = UW.b
#            
#        r = b2/b1
#        b  = (UW.b + LW.b)/2
##        Gb = G/b
#        Gb = G
#        
#        sigma1 = (1 - 0.66 * Gb)/(1.055 + 3.7 * Gb)
#        tau = (1 - r)/(1 + r)
#        s = 0.8 * sigma1 * (1 - sigma1) - 0.1
#        t = 0.56 / (sigma1 + s - 0.22)
#        sigma = sigma1 + s - (s**2 + (tau/t)**2)**0.5
#        kappa = (1 - sigma**2) / (r * (r + (1 / r) - (2 * sigma)))
#        
#        k = 1 / kappa**0.5
#        
#        return k, kappa

#===============================================================================
    def O_eff(self):
        """
        Calculates the oswald efficiency based upon the Gap/(lower wing span) ratio using an interpolation
        """
        if self.dirty: self.Refresh()
        
        gap = self.Gap
        
        return self.param.O_eff(gap)
        
#===============================================================================
    def _UpdateO_eff(self):
        """
        Updates the O_eff vs. Gap/(lower wing span) interpolation
        """
        self.param.O_eff = interp1d(self.GapInterp, self.OeffInterp, kind='linear')

#===============================================================================
    def BWCF(self):
        """
        Calculates the biwing correction factor based upon the Gap/(lower wing span)
           ratio using an interpolation
        """
        if self.dirty: self.Refresh()

        return self.param.BWCF(self.Gap)
        
#===============================================================================
    def _UpdateBWCF(self):
        """
        Updates the biwing correction factor vs. Gap/(lower wing span) interpolation
        """
        self.param.BWCF = interp1d(self.GapInterp, self.BWCFInterp, kind='linear')

#===============================================================================
    def _UpdateStagger(self):
        """
        Updates the stagger of the BiWing
        """
        #
        # This is really a 'Refreshing' operation, but it also needs to be called
        # from Xac()
        #
        unrefresh = False
        if not self.param.refreshing:
            unrefresh = True
            self.param.refreshing = True
        
        LW = self.LowerWing
        UW = self.UpperWing
        
        #
        # Lower Wing is reference wing, so update first.
        #
        for i in range(3):
            LW.X[i] = self.X[i]
                
        #
        # Update the position of the upper wing
        #
        UW.X[0] = self.X[0] - (self.Stagger * self.LowerWing.MAC())
        UW.X[2] = self.X[2] + (self.Gap * LW.b)
        
        if unrefresh:
            self.param.refreshing = False
        
#===============================================================================
    def _GetControls(self):
        """
        Creates a class of the control surfaces on the wings
        """
        class ControlsClass:
            pass
        Controls = ControlsClass()
        
        def AddControls(LiftSurf):
            for control in LiftSurf.Controls.__dict__.itervalues():
                Controls.__dict__[control.name] = control
            
        AddControls(self.LowerWing)
        AddControls(self.UpperWing)
        
        if self.EndPlate is not None:
            AddControls(self.EndPlate)
            
        return Controls
        
#===============================================================================
    def _CalcSpan(self):
        """
        Calculates span of the biwing as the average of the two spans
        """    
        if self.dirty: self.Refresh()
        
        UW = self.UpperWing
        LW = self.LowerWing
        
        return (UW.b + LW.b)/2

#===============================================================================
    def _CalcArea(self):
        """
        Calculates Effective Area of Biplane
            This is the area of the monowing with the equivalent coefficient of drag
        """
        if self.dirty: self.Refresh()
            
        UW = self.UpperWing
        LW = self.LowerWing
        
        return UW.S + LW.S
    
#===============================================================================
    def _CalcAR(self):
        """
        Calculates Effective Aspect Ratio of Biplane
            Taken to be the inverse of Munk's Area Ratio  (see Munk: General Biplane Theory)
        """
        if self.dirty: self.Refresh()
        
        UW = self.UpperWing
        LW = self.LowerWing
        b  = (UW.b + LW.b)/2
        
        return b**2 / (UW.S + LW.S)

#===============================================================================
    def _CalcWeight(self):
        """
        Calculates the weight of the biwing.
        """
        #
        # This is called from ACStructural.__getattr__
        #
        Weight = self.LowerWing.Weight
        Weight += self.UpperWing.Weight
        
        if self.EndPlate is not None:
            Weight += self.EndPlate.Weight #This will account for a symmetric endplate as well

        return Weight

#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACBiWing, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'S':
            return self._CalcArea()
        elif key == 'AR':
            return self._CalcAR()
        elif key == 'b':
            return self._CalcSpan()
        elif key == 'Controls':
            return self._GetControls()
        #
        # If no calculation exist return None
        #
        return None

#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.
        """        
#        self._CheckEquation(['AR'], Need = 0) #Must never specify an aspect ratio, it is always calculated
        
#===============================================================================
    def WriteAVLWing(self, filename):
        """
        Writes out an AVL input file for the lifting surface

        Inputs:
            filname - The file name of the AVL file to be written
        """
        if self.dirty: self.Refresh()
        
        fAVL = open(filename, 'w')
        
        mach = 0. #Mach number is just hard coded to zero for now...

        iYsym = 0  # This is an integer
        iZsym = 0  # This is an integer
        Zsym  = 0. # This is a floating point

        Sref = self.S / (IN**2)
        Cref = self.LowerWing.Chord(0*IN) / (IN)
        Bref = self.b / (IN)

        Xref = self.X[0] / (IN)
        Yref = 0.0 if self.Symmetric else self.X[1] / (IN)
        Zref = self.X[2] / (IN)

        CDp = 0.

        fAVL.write(self.name + '\n')
        fAVL.write(str(mach)  + '\t\t\t\t| Mach\n')
        fAVL.write(str(iYsym) + ' ' +  str(iZsym) + ' ' +  str(Zsym) + '\t\t\t| iYsym  iZsym  Zsym\n')
        fAVL.write(str(Sref)  + ' ' +  str(Cref)  + ' ' +  str(Bref) + '\t\t\t| Sref   Cref   Bref\n')
        fAVL.write(str(Xref)  + ' ' +  str(Yref)  + ' ' +  str(Zref) + '\t\t\t| Xref   Yref   Zref\n')
        fAVL.write(str(CDp)   +                                        '\t\t\t\t| CDp  (optional)\n')

        self._WriteAVLSurface(fAVL)

        fAVL.close()
#===============================================================================
    def _WriteAVLSurface(self, fAVL):
        """
        Writes out the bi-wing to an AVL file

        Inputs:
            fAVL - The open file to write to
        """
        self.LowerWing._WriteAVLSurface(fAVL)
        self.UpperWing._WriteAVLSurface(fAVL)
         
        if self.EndPlate is not None:
            self.EndPlate._WriteAVLSurface(fAVL, Sspace = -1.0) #Cluster the spacing at both ends of the endplate
        

#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """
        Draws this lifting surface and it's control surfaces

        Inputs:
            fig   - Integer number of figure to draw into
            top   - Subplot for top view
            side  - Subplot for side view
            front - Subplot for front view
        """
        super(ACBiWing,self).Draw(fig, top, side, front)

        self.UpperWing.Draw(fig, top, side, front)
        self.LowerWing.Draw(fig, top, side, front)
        if self.UpperWing.WingWeight is not None:
            self.UpperWing.WingWeight.Draw(fig, top, side, front)
        if self.LowerWing.WingWeight is not None:
            self.LowerWing.WingWeight.Draw(fig, top, side, front)
        if self.EndPlate is not None:
            self.EndPlate.Draw(fig, top, side, front)
        
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
        super(ACBiWing, self).Refresh()
        self.param.refreshing = True
        
        self._CheckConsistent()
        self._UpdateO_eff()
        self._UpdateBWCF()
        
        UW = self.UpperWing
        LW = self.LowerWing
        
        # Pass parameters
         

        UW.V_max_climb = LW.V_max_climb = self.V_max_climb
        UW.Alt_LO      = LW.Alt_LO      = self.Alt_LO
        UW.V_LOstall   = LW.V_LOstall   = self.V_LOstall
        UW.V_Stall     = LW.V_Stall     = self.V_Stall
        UW.FullWing    = LW.FullWing    = self.FullWing
        
        
        #
        # Determines the lift off lift of both wings assuming some lift ratio
        #
        if self.Lift_LO is not None:
            LR = self.Lift_Ratio
            UW.Lift_LO = self.Lift_LO * (1 - LR)
            LW.Lift_LO = self.Lift_LO * LR
            
        if self.TotalS is not None:
            LR = self.Lift_Ratio
            UW.S = self.S * (1 - LR)
            LW.S = self.S * LR    
            
        #
        # If an overall span is specified, the upper and lower wings are equal
        #
        if self.b is not None:
            UW.b = LW.b = self.b
            
        self._UpdateStagger()
        
        #
        # Now update the drawings
        #
        self.UpdateDrawing()
        
        #
        # Now that everything is up to date, update the Alpha2d calculation
        #
        self._UpdateAlpha2d_LO()
        
        self.param.refreshing = False

#===============================================================================
    def UpdateDrawing(self):
        """
        Updates all drawings of the wings and the endplates
        """
        self.LowerWing.UpdateDrawing()
        self.UpperWing.UpdateDrawing()

        if self.EndPlate is not None:
            self.EndPlate.Refresh() #Endplate must actually be refreshed

#===============================================================================
    def Draw2DAirfoilPolars(self, fig = 1):
        """
        Draws the lift and drag polars for the airfoil

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        UW = self.UpperWing
        LW = self.LowerWing
        
        #
        # Create the angles of attack
        #
        alpha = self.AlphaRange()

        #
        # Create the polars
        #
        Cl = self.CL(alpha)
        Cd = UW.af.Cd(alpha, UW.Re()) + LW.af.Cd(alpha, LW.Re())
        Cm = self.CM(alpha)
        
        ClSlopeAt = self.ClSlopeAt[0]
        CmSlopeAt = self.CmSlopeAt[0]

        dCl_da = self.dCL_da() / (1/ARCDEG)
        dCm_da = self.dCM_da() / (1/ARCDEG)
         
        #
        # Make the slope go through the specified ClSlopeAt. Otherwise let it go
        # though the Cl = 0 point.
        #
        ClOffset  = self.CL(ClSlopeAt) - dCl_da*ClSlopeAt / (ARCDEG)
        CmOffset  = self.CM(CmSlopeAt) - dCm_da*CmSlopeAt / (ARCDEG)
        
        #
        # These are in degrees without being a scalar
        #
        alp_da = npy.array([alpha[0] / (ARCDEG), alpha[-1] / (ARCDEG)])
        
        Cl_da  = alp_da*dCl_da + ClOffset
        Cm_da  = alp_da*dCm_da + CmOffset

        self._PlotPolars('2D Airfoil Polars', fig, alpha, Cl, Cd, Cm, Cl_da, Cm_da, alp_da)

#===============================================================================
    def Draw3DWingPolars(self, fig = 1):
        """
        Draws the lift and drag polars for the airfoil

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        #
        # Create the angles of attack
        #
        alpha = self.AlphaRange()

        #
        # Create the polars. Call ACLiftSurf other wise Tail surfaces will not be properly generated.
        #
        CL = self.CL(alpha)
        CD = self.CD(alpha)
        CM = self.CM(alpha)
        
#        ClSlopeAt = self.ClSlopeAt
#        CmSlopeAt = self.CmSlopeAt
#
#        dCl_da = self.dCL_da() / (1/ARCDEG)
#        dCm_da = self.dCM_da() / (1/ARCDEG)
#        
#        #
#        # Make the slope go through the specified ClSlopeAt and CmSlopeAt.
#        #
#        ClOffset  = self.CL(ClSlopeAt) - dCl_da*self.AlphaWing(ClSlopeAt) / (ARCDEG)
#        CmOffset  = self.CM(CmSlopeAt) - dCm_da*self.AlphaWing(CmSlopeAt) / (ARCDEG)
#        
#        #
#        # These are in degrees without being a scalar
#        #
#        alp_da = npy.array([alpha[0] / (ARCDEG), alpha[-1] / (ARCDEG)])
#        alp_da = self.AlphaWing(alp_da*ARCDEG) / (ARCDEG)
#        
#        Cl_da  = alp_da*dCl_da + ClOffset
#        Cm_da  = alp_da*dCm_da + CmOffset
        
        #
        # Convert the 2D angles if attack to 3D
        #
        alpha = self.AlphaWing(alpha)
        
        #
        # Plot the polars
        #
        self._PlotPolars('3D Wing Polars', fig, alpha, CL, CD, CM) #, Cl_da, Cm_da, alp_da)
        
#===============================================================================
    def _PlotPolars(self, title, fig, alpha, CL, CD, CM, Cl_da = None, Cm_da = None, alp_da = None):
        """
        An internal method for plotting lift and drag polars
        
        Inputs:
            title  - The title for the plot
            alpha  - Angles of attack
            CL     - CL curve
            CD     - CD curge
            CM     - CM curve
            Cl_da  - Two points definig a line with the analysis Cl slope
            alp_da - Two angles of attack for Cl_da
        """
        
        alpha = alpha / (ARCDEG)
        #
        # Draw them
        #
        pyl.subplot(221)
        pyl.plot(alpha, CL)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_l$')
        #
        # Also draw the line prepresentative of the CL slope
        #
        if Cl_da is not None:
            pyl.plot(alp_da, Cl_da)
            pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_l$')
            pyl.legend([r'$C_l$',r'$C_l$ Slope'], loc = 'best')


        pyl.subplot(222)
        pyl.plot(alpha, CD)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_d$')

        pyl.subplot(223)
        pyl.plot(alpha, CM)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_m$')
        #
        # Also draw the line prepresentative of the CL slope
        #
        if Cm_da is not None:
            pyl.plot(alp_da, Cm_da)
            pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_m$')
            pyl.legend([r'$C_m$',r'$C_m$ Slope'], loc = 'best')

        pyl.subplot(224)
        pyl.plot(CD, CL)
        pyl.xlabel(r'$C_d$'); pyl.ylabel(r'$C_l$')

#        pyl.figure(fig+42)
#        pyl.plot(alpha,CL/CD)
#        pyl.xlabel(r'$\alpha$'); pyl.ylabel('L/D')

        title += ' Re = %1.0f' % self.Re()
        pyl.annotate(title, xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)
                                
################################################################################
class ACEndPlate(ACLiftSurf):
    """
    A class for creating end plates for the box wing
    
    Note:
        AR does not exist for a winglet because the root and tip chord
          is determined from the main wings
        The span is the gap between the two wings for an equal span biwing,
          otherwise it is sqrt(gap**2 + ((b1-b2)/2)**2)
        The overall dihedral and sweep are determined by the stagger and the
          span difference between the biwings, respectively
    """
    def __init__(self, index, BiWing):
        """
        Inputs:
            index  - Index of the end plate for AVL
            BiWing - the biwing the end plat is attached to
        """
        super(ACEndPlate, self).__init__(index)
        
        self.Axis = (0, 1) #Default to a vertical end plate (A horizontal one does not make much sense...)
        
        #
        # Save of a reference to the main wings
        #
        self.param.BiWing = BiWing
        
        self.param.index = index
        
        self.TR = 1
        self.Fb = 0.5
        self.Gam = 0*ARCDEG
        self.Lam = 0*ARCDEG

#===============================================================================
    def Re(self):
        """
        Calculates the Reynolds number for the wing.

        This is a relatively cheap calculation and is therefore left
        as a calculation rather than storing a self.param.Re
        to reduce the interdependencies.
        """
        Alt = self.param.LWing.GetAlt_LO()
        V   = self.GetV_LO()
        rho = AeroUtil.DensAlt(Alt)
        mu  = AeroUtil.Mu(AeroUtil.TempAlt(Alt))
        mac = self.MAC()

        return rho * V * mac / mu

#===============================================================================
    def GetV_LO(self):
        """
        Retreives the lift of velocity for the end plate. This was taken from the 
        lower wing to begin with. The method GetV_LO allows a winglet to 'pretend'
        that it is a 'lower wing'.
        """
        return self.param.BiWing.LowerWing.GetV_LO()

#===============================================================================
    def GetAlt_LO(self):
        """
        Retreives the lift of altitude for the end plate. This was taken from the 
        lower wing to begin with. The method GetAlt_LO allows a winglet to 'pretend'
        that it is a 'lower wing'.
        """
        return self.param.BiWing.LowerWing.GetAlt_LO()
        
#===============================================================================
    def Refresh(self):
        """
        Updates all the internal quantities.
        """
        self.__dict__['dirty'] = False
        self.param.refreshing = True
        #
        # Due to the cyclic behavior of the area calculation, it must be updated before anything else.
        #
        self._CalcChord()
        self._UpdateArea()
        self._CalcSweep()
        self._CalcDihedral()
        
        #
        # Set the winglet to be at the tip of the wing.
        # This is done before calling the baseclass where it will be drawn
        #
        self.X = self.param.BiWing.LowerWing.Tip()
              
        #
        # Now update everything in the base class
        #
        super(ACEndPlate, self).Refresh()
       
        self.param.refreshing = False    

#===============================================================================
    def _GetChord(self):
        """
        Calculates the array of chords. This is simple because the root and tip chord is 
        given by the tip chord of the two wings. The returned list is not a scalar but a 
        list of numbers in units of inches.
        """
        LWing = self.param.BiWing.LowerWing
        UWing = self.param.BiWing.UpperWing
        b1,b2 = self._Wing_WingTip()
        
        TR  = self.TR
        
        try: l = len(TR)
        
        except: TR = [TR]
            
        #
        # Note the the chord locations are stored in inch units
        #
        C1 = Cr = LWing.Chord(b1) / IN   #Root chord
        Ct      = UWing.Chord(b2) / IN   #Root chord
        chord = npy.empty(len(TR)+1)
        chord[0] = Cr
        chord[len(TR)] = Ct

        # Calculate remaining chords except the last chord, which is the upper wing tip chord
        for i in range(len(TR)-1):
            C1 *= TR[i]
            chord[i+1] = C1

        return chord

#===============================================================================
    def _CalcChord(self):
        """
        Fills the internal chord array

        This is overloaded from ACLiftSurf because the root chord is given by the main wing
        """
        self.param.chord = self._GetChord()
            
#===============================================================================
    def _CalcArea(self):
        """
        Returns the calculated area
        """
        if self.dirty: self.Refresh()
        return self.param.S

#===============================================================================
    def _UpdateArea(self):
        """
        Calculates the area of the end plate
        """
        
        if not self.NoneList.has_key('S'):
            return
        #
        # The area is calculated based on the tip chord of the two wings and
        # the taper ratios. Because the area is needed for so many other calculations,
        # the chord is also calculated here so the area can be updated before the base class is Refreshed.
        #
        
        b  = self.b / IN
        TR = self.TR
        Fb = self.Fb
        
        #
        # Check that the array sizes are correct
        #
        TR, Fb  = self._CheckSameLen(TR, "Taper ratio", Fb, "span fraction")

        #
        # Get the chords
        #
        chord = self._GetChord()

        #
        # Sum up the trapezoidal areas
        #
        S = (chord[0] + chord[1])/2.*Fb[0]*b
        for i in range(1,len(Fb)):
            S += (chord[i] + chord[i+1])/2.*(Fb[i]-Fb[i-1])*b
        
        #
        # Save of the area of the end plate
        #
        self.param.S = S*IN**2

#===============================================================================
    def _CalcSpan(self):
        """
        Returns the calculated span
        """
        
        LWb   = self.param.BiWing.LowerWing.b
        LWFb  = self.param.BiWing.LowerWing.Fb
        LWGam = self.param.BiWing.LowerWing.Gam
        UWb   = self.param.BiWing.UpperWing.b
        UWFb  = self.param.BiWing.UpperWing.Fb
        UWGam = self.param.BiWing.UpperWing.Gam
        
        #
        # Check that the array sizes are correct
        #
        LWGam, LWFb  = self._CheckSameLen(LWGam, "Dihedral", LWFb, "span fraction")
        UWGam, UWFb  = self._CheckSameLen(UWGam, "Dihedral", UWFb, "span fraction")
        
        # Calculate the change in span due to the dihedral of the wings
        DelbGamLW = LWFb[0]*LWb/2.*math.tan(LWGam[0] / RAD).real
        for i in range(1, len(LWFb)):
            DelbGamLW += (LWFb[i]-LWFb[i-1])*LWb/2.*math.tan(LWGam[i] / RAD).real
        
        DelbGamUW = UWFb[0]*UWb/2.*math.tan(UWGam[0] / RAD).real
        for i in range(1, len(UWFb)):
            DelbGamUW += (UWFb[i]-UWFb[i-1])*UWb/2.*math.tan(UWGam[i] / RAD).real
        
        return self.param.BiWing.Gap*LWb + DelbGamUW - DelbGamLW

#===============================================================================
    def _CalcDihedral(self):
        """
        Returns the calculated span
        """
        
        LWb = self.param.BiWing.LowerWing.b
        UWb = self.param.BiWing.UpperWing.b
        
        #
        # Calculate the dihedral based on the given angles
        #
        b   = self.b / IN
        Gam = self.Gam
        Fb  = self.Fb
        
        #
        # Check that the array sizes are correct
        #
        Fb, Gam = self._CheckSameLen(Fb, "Span fraction", Gam, "dihedral angles")
        
        b = b/2. if self.FullWing else b
        
        dihed = self.param.dihed = npy.zeros(len(Gam)+1)
        
        dihed[0] = 0 # Never offset the root chord
        dihed[1] = b * Fb[0] * math.tan(Gam[0] / RAD).real
        for i in range(1,len(Gam)-1):
            dihed[i+1] = dihed[i] + b * (Fb[i]-Fb[i-1]) * math.tan(Gam[i] / RAD).real
            
        dihed[-1] = ((UWb - LWb)/2.) / IN

#===============================================================================
    def _CalcSweep(self):
        """
        Updates the sweep offset used to move chord locations forward/back
        """

        #
        # Change lambda if a constant edge is desired
        #
        if self._CalcConstantEdge():
            return

        BW    = self.param.BiWing
        LWswp = self.param.BiWing.LowerWing.param.sweep
        UWswp = self.param.BiWing.UpperWing.param.sweep
        
        b     = self.b / (IN)
        Lam   = self.Lam
        Fb    = self.Fb
        chord = self.param.chord
        sfc   = self.SweepFc

        b = b/2. if self.FullWing else b

        #
        # Check that the lists are the same lenght
        #
        Fb, Lam = self._CheckSameLen(Fb, "Span fraction", Lam, "sweep angles")
        
        Lam = [Lam[i] / (RAD) for i in range(len(Lam))]

        sweep = self.param.sweep = npy.zeros(len(Lam)+1)
        #
        # First sweep the quarter chord back to get a constant line
        # along the desired chord fraction to sweep about
        #
        for i in range(len(Fb)-1):
            sweep[i+1] = (0.25 - sfc)*(chord[i+1] - chord[i])

        #
        # Now apply the disreid sweep
        #
        sweep[0] = 0 # Never offset the root chord
        sweep[1] += b * Fb[0] * math.tan(Lam[0]).real
        for i in range(1,len(Lam)-1):
            sweep[i+1] += sweep[i] + b * (Fb[i]-Fb[i-1]) * math.tan(Lam[i]).real
            
        sweep[-1] += UWswp[-1]-LWswp[-1]-(BW.Stagger*BW.MAC()) / (IN)

#===============================================================================
    def _CalcConstantEdge(self):
        """
        Sets the sweep to get a constant LE or TE if desired
        """

        #
        # 'chord' must have been updated before calling this method
        #
        if self.CEdge != 'LE' and self.CEdge != 'TE':
            return False #Nothing was asked for so do nothing

        b     = self.b / IN
        Fb    = self.Fb
        chord = self.param.chord

        b = b/2. if self.FullWing else b

        # Calculate sweep angle due to main wing stagger and/or sweep
        BW    = self.param.BiWing
        LWswp = self.param.BiWing.LowerWing.param.sweep
        UWswp = self.param.BiWing.UpperWing.param.sweep
        gap     = (BW.Gap*BW.LowerWing.b) / IN

        sweep = self.param.sweep = npy.zeros(len(Fb)+1)

        sweep[0] = 0 # Never offset the root chord

        chordswp = 0
        if self.CEdge == 'LE':
            swpstag = UWswp[-1]-LWswp[-1]-(BW.Stagger*BW.MAC()) / IN - 0.25*(chord[-1]-chord[0])
            lam     = math.atan(swpstag/gap).real
            for i in range(len(Fb)):
                chordswp += 0.25*(chord[i+1] - chord[i])
                sweep[i+1] = chordswp + b*Fb[i]*math.tan(lam).real

        elif self.CEdge == 'TE':
            swpstag = UWswp[-1]-LWswp[-1]-(BW.Stagger*BW.MAC()) / IN + 0.75*(chord[-1]-chord[0])
            lam     = math.atan(swpstag/gap).real
            for i in range(len(Fb)):
                chordswp  += -1*0.75*(chord[i+1] - chord[i])
                sweep[i+1] = chordswp + b*Fb[i]*math.tan(lam).real

        return True


#===============================================================================
    def __setattr__(self, key, value):
        """
        Attempts to calculate any variables set to None
        """
        
        if key == 'TR' or key == 'Fb':
            if isinstance(value,list):
                value += [1]
            else:
                value = [value,1]
        
        if key == 'Gam' or key == 'Lam':
            if isinstance(value,list):
                value += [0*ARCDEG]
            else:
                value = [value,0*ARCDEG]
        
        super(ACEndPlate, self).__setattr__(key,value)
        
#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.
        """
        super(ACEndPlate, self)._CheckConsistent()
        
        # Aspect ratio and span should never be specified for an end plate
        self._CheckEquation(['AR'], Need = 0)
        self._CheckEquation(['b'], Need = 0)
            
#===============================================================================
    def _Wing_WingTip(self):
        """
        Returns span location of the wing-tip of the two wings
        """
        b1 = self.param.BiWing.LowerWing.b
        b2 = self.param.BiWing.UpperWing.b
        b1 = b1/2. if self.param.BiWing.FullWing else b1
        b2 = b2/2. if self.param.BiWing.FullWing else b2
        return b1,b2

################################################################################
if __name__ == '__main__':
    # Create a wing and changing some parameters
#    wing = ACMainWing(1)
#
#    wing.Lift_LO       = 27.00654 * LBF
#    wing.V_max_climb   = 65 * FT/SEC
#    wing.V_Stall       = 41.2699923511 * FT/SEC
##    wing.V_Stall_Init   = 39.2 * FT/SEC
#    wing.Alt_LO        = 600 * FT
#
##    wing.S = 900 * IN**2
#    wing.AR = 16
##    wing.b = 10. * FT
#    wing.TR = [1, 0.7, 0.7]
#    wing.Fb = [0.3, 0.8, 1]
#    wing.Gam = [0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG]
#    wing.CEdge = 'LE'
#    wing.ConstUpper = True
#    wing.Airfoil = ('../Airfoils/S1223/S1223.dat','../Airfoils/S1223')
#    wing.FullWing = True
#    wing.ClSlopeAt = 7 * ARCDEG
#    wing.CmSlopeAt = 8 * ARCDEG
##    wing.color = 'r'
#    wing.LLFile = 'LLTPolar.txt'
#    
#    wing.X = [10.0 * IN, 0.0 * IN, 10.0 * IN]
#
#    ail = wing.AddControl('ail')
#
#    wing.AddWinglet("Winglet",2)
#    Winglet = wing.Winglets.Winglet
#    Winglet.b = 10 *IN
#    Winglet.Airfoil = ('../Airfoils/NACA0012/NACA0012.dat','../Airfoils/NACA0012')
#    Winglet.Lam = [20*ARCDEG, 20*ARCDEG]
#    Winglet.Gam = [20*ARCDEG, 20*ARCDEG]
#    Winglet.Fb  = [0.6, 1]
#    Winglet.TR  = [0.9, 0.8]
#    Winglet.SweepFc = 0
#    Winglet.Symmetric = True
##    Winglet.FullWing = True
#    
#    Winglet.AddWinglet("CWinglet", 3)
#    CWinglet = Winglet.Winglets.CWinglet
#    CWinglet.name = 'CWinglet'
#    CWinglet.Axis = (-1,0)
#    CWinglet.b = 10 *IN
#    CWinglet.Airfoil = ('../Airfoils/NACA0012/NACA0012.dat','../Airfoils/NACA0012')
##    CWinglet.Airfoil = ('../Airfoils/S1223/S1223.dat','../Airfoils/S1223')
#    CWinglet.Lam = [10*ARCDEG, 20*ARCDEG]
#    CWinglet.Gam = [0*ARCDEG, 0*ARCDEG]
#    CWinglet.Fb  = [0.6, 1]
#    CWinglet.TR  = [0.8, 0.8]
#    CWinglet.SweepFc = 0
#    CWinglet.Symmetric = True
##    CWinglet.FullWing = True
    
#    CWinglet.AddWinglet("WingletDrop",4)
#    WingletDrop = CWinglet.Winglets.WingletDrop
#    WingletDrop.name = 'WingletDrop'
#    WingletDrop.Axis = (0,-1)
#    WingletDrop.b = 9.9 *IN
#    WingletDrop.Airfoil = ('../Airfoils/NACA0012/NACA0012.dat','../Airfoils/NACA0012')
#    WingletDrop.Lam = [10*ARCDEG, 20*ARCDEG]
#    WingletDrop.Gam = [0*ARCDEG, 0*ARCDEG]
#    WingletDrop.Fb  = [0.6, 1]
#    WingletDrop.TR  = [0.8, 0.8]
#    WingletDrop.SweepFc = 0
#    WingletDrop.Symmetric = True
#    CWinglet.FullWing = True

#    wing.WriteAVLWing('AVLMainWing.avl')
    
    
    biwing = ACBiWing(1, 2, 3)

    biwing.Lift_LO       = 27.00654 * LBF
    biwing.Lift_Ratio    = 0.5
    biwing.V_max_climb   = 65 * FT/SEC
    biwing.V_Stall       = 41.2699923511 * FT/SEC
    biwing.Alt_LO        = 600 * FT
    
#    biwing.b             = 6 * FT
    biwing.Gap           = 0.2
    biwing.Stagger       = 0.5
    biwing.UpperWing.b   = 6*FT
    biwing.LowerWing.b   = 4*FT
    biwing.UpperWing.Airfoil = 'S1223'
    biwing.LowerWing.Airfoil = 'S1223'
#    biwing.LowerWing.SweepFc = 1
    biwing.LowerWing.TR      = [0.9,0.8]
    biwing.LowerWing.Gam     = [0*ARCDEG, 0*ARCDEG]
    biwing.LowerWing.Lam     = [-5*ARCDEG, -5*ARCDEG]
    biwing.LowerWing.Fb      = [0.6,1]
    biwing.UpperWing.TR      = [1,1]
    biwing.UpperWing.Gam     = [0*ARCDEG, 0*ARCDEG]
    biwing.UpperWing.Lam     = [10*ARCDEG, 10*ARCDEG]
    biwing.UpperWing.Fb      = [0.6,1]
    biwing.FullWing = True
    
    biwing.CreateEndPlate()
    #used for conservative drag estimation
    biwing.EndPlate.Airfoil = 'NACA0012'
    
    #DO NOT specify an Fb of 1 for the end plate!!!
    #An Fb of 1 is at the Upper Wing and does not need to be specified
    biwing.EndPlate.Fb      = [0.2,0.5,0.8]
    biwing.EndPlate.TR      = [0.8,0.25,4]
    biwing.EndPlate.Gam     = [-30*ARCDEG,60*ARCDEG,30*ARCDEG]
    biwing.EndPlate.Lam     = [0*ARCDEG,0*ARCDEG,0*ARCDEG]
    biwing.EndPlate.SweepFc = 1
    biwing.EndPlate.CEdge   = 'LE'
    biwing.EndPlate.Symmetric = True
    
    biwing.WriteAVLWing('AVLBiWing.avl')
    

#    alpha2dw = 5 * ARCDEG

#    print 'Fuselage angle', wing.AlphaFus(alpha2dw)
#
#    print 'Wing AC', wing.Xac()
#    print 'Wing MAC', AsUnit( wing.MAC(), 'in' )
#
#    print 'Wing Cl', wing.af.Cl(alpha2dw,wing.Re())
#    print 'Wing AR', wing.AR
#    print 'Wing area: ', AsUnit( wing.S, 'in**2' )
#    print 'Wing Span', AsUnit( wing.b, 'ft' )
#    print 'Wing Lift_LO', AsUnit( wing.Lift_LO, 'lbf' )
#    print 'Wing V_Stall', AsUnit( wing.V_Stall, 'ft/s' )
#    print 'Wing MaxTE', AsUnit( wing.MaxTE(), 'in' )
#      
#    print 'Schreiner load', AsUnit( wing.LoadDistSchreiner(2* IN, 30 * FT/SEC), 'lbf/in' )
#    print 'ShearLoading!!!  ', wing.ShearLoading(5*FT, 30 * FT/SEC)
    
#    alphas = wing.AlphaRange()
#    print 'alphas', alphas
#    print 'Several Horizontal Tail angle', htail.Alpha2DHTail(5.*ARCDEG, del_c)

#    alphas = wing.AlphaRange(5)
#    alpha2dh = htail.Alpha2DHTail(alphas, 0 * ARCDEG) / (ARCDEG)
#    alphafus = wing.AlphaFus(alphas) / (ARCDEG)

#    pyl.plot(alphafus, alpha2dh)
#    pyl.show()

#    alphas = wing.AlphaRange(5)
#    hCL = htail.CL(alphas, 0 * ARCDEG)
#    alphafus = wing.AlphaFus(alphas) / (ARCDEG)

#    pyl.plot(alphafus, hCL)
#    pyl.show()

#    wing.LLTPlots(2)

#    wing.Draw2DAirfoilPolars(4)
#    wing.Draw3DWingPolars(3)
#    wing.Draw()
    biwing.Draw()

    pyl.show()

