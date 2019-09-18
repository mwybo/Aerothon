"""
A class for constructing lifting surfaces
"""
from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import OverSpecifiedError, ACComponent, ACStructural, Length, Area, Angle, Force, Unitless, g, Torque, MassTorque
from ACAirfoil import ACAirfoil
from ACMass import ACMassBox
from ACMaterial import ACMaterial
from ACAVLWing import ACAVLWing
from ACLiftingLine import ACLiftingLine
from ACXFoil import ACXFoil
import AeroUtil
from scipy.interpolate import interp1d, RectBivariateSpline, InterpolatedUnivariateSpline
from scipy import integrate
import numpy as npy
import pylab as pyl
import cmath as math
import operator as op
import os
from scalar.units import FT, IN, ARCDEG, RAD, LBM, LBF, SEC, OZF, OZM
from scalar.units import AsUnit

################################################################################
class ACLiftSurfError(Exception):
    def __init__(self, message, name):
        self.msg = message
        self.name = name
    def __str__(self):
        return self.msg + "'" + self.name + "'"

#===============================================================================
class ACControlSurfError(Exception):
    def __init__(self, message, name):
        self.msg = message
        self.name = name

    def __str__(self):
        return repr(self.msg + "'" + self.name + "'")

################################################################################
class ACControlSurf(ACStructural):
    """
    A class for constructing control surfaces

        Fc        - Fraction of chord
        Fb        - Fraction of span
        Ft        - Fraction of span to determine distance from lifting surface tip
        Ht        - Fraction of thickness for the hinge location
        Gain      - Gain of the control surface
        SgnDup    - Sign of the deflection of a duplicate surface for a symmetric lifting surface
                    i.e. for an aileron (-1) and elevator (1)
    """
    
#===============================================================================
    class CSServo(ACMassBox):
        """
        A class for managing servo's for the control surface
        
        Attributes:
            Fc  - Lifting surface chord fraction location of the servo
            Fbc - Control surface span fraction location of the servo
        """

#===============================================================================
        def __init__(self, Control):
            super(ACControlSurf.CSServo, self).__init__()
            UnitList = self.UnitList
            
            self.name               = 'Servo'
            self.__dict__['Fc']     = 0.5         ;   UnitList['Fc']     = Unitless
            self.__dict__['Fbc']    = 0.5         ;   UnitList['Fbc']    = Unitless
            self.__dict__['Torque'] = 1.0*OZF*IN  ;   UnitList['Torque'] = MassTorque
            
            self.Xat = [0, 0.5, 0.5]
            #
            # Make the mass box brown
            #
            self.Color = 'brown'
            
            self.param.Control = Control
            
#===============================================================================
        def Refresh(self):
            """
            Updates the position of the servo
            """
            self.param.refreshing = True
            
            Control = self.param.Control
            Wing    = Control.param.wing
            
            self.CopyVisibility(Wing)
            
            self.Symmetric = Control.Symmetric
                       
            SpanLoc  = Control.Root() + Control.Span()*self.Fbc
            ChordLoc = Wing.Chord(SpanLoc)*self.Fc
            
            LE = Wing.LE(SpanLoc)
            
            StreamLoc = LE + ChordLoc
            MidLoc    = Wing.Mid(SpanLoc, LE + ChordLoc)
            
            if abs(Wing.Axis[0]) == 1: #Horizontal lifting surface
                self.X = [StreamLoc, SpanLoc + Wing.X[1]*Wing.Axis[0], MidLoc]
            else:                      #Vertical lifting surface
                self.X = [StreamLoc, MidLoc, SpanLoc + Wing.X[2]*Wing.Axis[1]]
    
            self.param.refreshing = False
            #
            # Now that the position is updated, update the MassBox
            #
            super(ACControlSurf.CSServo, self).Refresh()
            
#===============================================================================
        def Mirror(self, coordinate):
            """
            Mirrors the coordinate appropriately
            """
            Control = self.param.Control
            wing    = Control.param.wing
            
            if abs(wing.Axis[0]) == 1: #Horizontal wing
                ax = 1
            else:                      #Vertical wing
                ax = 2
            
            coordinate = npy.array([coordinate[i] / IN for i in range(3)])*IN
            coordinate[ax] = -coordinate[ax] + 2*wing.X[ax]*wing.Axis[ax-1]
            return coordinate
            
#===============================================================================
    def __init__(self, name, wing):
        """
        Inputs:
             name - Name of the control surface
             wing - Wing which the surface belongs to
        """
        super(ACControlSurf, self).__init__()
        UnitList = self.UnitList
        
        self.name = name
        
        self.__dict__['Fc']        = 0.25       ;   UnitList['Fc']        = Unitless
        self.__dict__['Fb']        = 0.3        ;   UnitList['Fb']        = Unitless
        self.__dict__['Ft']        = 0.2        ;   UnitList['Ft']        = Unitless
        self.__dict__['Ht']        = 0.5        ;   UnitList['Ht']        = Unitless
        self.__dict__['Gain']      = 1.0        ;   UnitList['Gain']      = Unitless
        self.__dict__['SgnDup']    = 1.0        ;   UnitList['SgnDup']    = Unitless

        #
        # Give the control surface some material
        #
        self.__dict__['ControlMat'] = ACMaterial()
        self.ControlMat.ForceDensity = 0*LBF/IN**3
        
        #
        # The servo for the control surface
        #
        self.__dict__['Servo'] = ACControlSurf.CSServo(self)
        self.Servo.Parents.append( self )
        
        #
        # Create the drawing object for the control surface
        #
        self.Drawings.Control = ACComponent.ComponentDrawing()
        self.Drawings.Control.color = 'brown'

        self.Drawings.ControlMir = ACComponent.ComponentDrawing()
        self.Drawings.ControlMir.color = 'brown'

        #
        # Save of the wing that this control surface belongs to
        #
        self.param.wing = wing

        self._InitKdelc()

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the control surface and its parts to the Weight Table
        
        Input:
            PartName    - The name given to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WeightTable[PartName] = self
        
        Servo = self.Servo
        Control = self
        
        #
        # Because all the servo inAsUnition is actually calculated in the
        # ACLiftSurf class, they have to calculated in a temporary class here
        #
        class CServo:
            def CG(self):
                CG  = Servo.CG()              
                if Control.param.wing.FullWing and Control.RootFb() > 0.0:
                    CG   = (CG + Servo.Mirror(CG))/2              
                return CG 
            def MOI(self):
                MOI = Servo.MOI()
                if Control.param.wing.FullWing and Control.RootFb() > 0.0:
                    MOI *= 2
                return MOI          
                
        
        tServo = CServo()
        tServo.name        = self.Servo.name
        tServo.Weight      = self.Servo.Weight
        tServo.WeightGroup = self.Servo.WeightGroup
        
        if self.param.wing.FullWing and self.RootFb() > 0.0:
            tServo.Weight *= 2

        if self.param.wing.Symmetric:
            tServo.Weight *= 2
         
        WeightTable[PartName]['Servo'] = tServo
        
#===============================================================================
    def _WingChord(self,y):
        """ Assists in the drawing process"""
        return self.param.wing.Chord(y) / IN

#===============================================================================
    def _WingTE(self,y):
        """ Assists in the drawing process"""
        return self.param.wing.TE(y) / IN

#===============================================================================
    def Root(self):
        """ Gives the root location of the control surface """
        if self.dirty: self.Refresh()
        wing = self.param.wing
        b = wing.b
        b = b/2. if wing.FullWing else b
        return self.RootFb() * b

#===============================================================================
    def Tip(self):
        """ Gives the tip location of the control surface """
        if self.dirty: self.Refresh()
        wing = self.param.wing
        b = wing.b
        b = b/2. if wing.FullWing else b
        return self.TipFb() * b

#===============================================================================
    def RootChord(self):
        """ Gives the root chord of the control surface """
        return self.param.wing.Chord(self.Root())

#===============================================================================
    def TipChord(self):
        """ Gives the tip chord of the control surface """
        return self.param.wing.Chord(self.Tip())

#===============================================================================
    def RootFb(self):
        """ Gives the control surface span fraction location closes to the fuselage"""
        Fb = self.Fb
        Ft = self.Ft
        return (1. - (Fb + Ft))

#===============================================================================
    def TipFb(self):
        """ Gives the control surface span fraction location furthest away from the fuselage"""
        Ft = self.Ft
        return (1. - Ft)

#===============================================================================
    def MeanChord(self):
        """  Calculates the mean chord of the control surface """
        if self.dirty: self.Refresh()
        
        wing = self.param.wing
        
        rSpan  = self.Root()
        tSpan  = self.Tip()
        rChord = wing.Chord(rSpan)
        tChord = wing.Chord(tSpan)
        
        return (tChord + rChord)/2.

#===============================================================================
    def Span(self):
        """ Gives the span of the control surface """
        if self.dirty: self.Refresh()
        wing = self.param.wing
        b = wing.b
        b = b/2. if wing.FullWing and self.RootFb() > 0.0 else b
        return self.Fb * b


#===============================================================================
    def HingeMoment(self, Chm, ArmRatio, V = None):
        """ 
        Calculates the required hinge moment
        
        Inputs:
            Chm      - Hinge moment coefficient from xfoil
            ArmRatio - Ratio of the servo arm to control surface arm lengths (Servo Arm/Control Arm). The lengths are measured to the hinge location.
            V        - Flight velocity
        """
        if self.dirty: self.Refresh()

        c = self.param.wing.RefLen()

        if V is None:
            V = self.param.wing.GetV_LO()
        
        Alt = self.param.wing.GetAlt_LO()
        rho = AeroUtil.DensAlt(Alt)
        b   = self.Span()
        
        return Chm * 0.5 * rho * V**2 * c**2 * b * ArmRatio

#===============================================================================
    def PlotHingeMoment(self, fig, RunDir, alpha2d, Ht, des, ArmRatio, Vs, Execute, ShowLegend = True):
        """ 
        Calculates the required hinge moment
        
        Inputs:
            fig     - Figure number
            RunDir  - Directory for XFoil files
            alpha2d - 2D airfoil angle of attack
            Ht      - Fraction of thickness location of the hinge (0 - bottom, 0.5 - mid, 1 - top)            
            des     - Control surface deflection range
            ArmRatio - Ratio of the servo arm to control surface arm lengths (Servo Arm/Control Arm). The lengths are measured to the hinge location.
            Vs      - Velocity range
            Execute - True/False if XFoil should be executed
            ShowLegend - Determines if a legend should be added
        """
        if self.dirty: self.Refresh()

        wing = self.param.wing
        Airfoil = self.param.wing.af.GetRealitveFilename()
        
        xfoil = ACXFoil()
        
        xfoil.RunDir = RunDir + os.sep
        
        n = 0
        for de in des:
            n += 1
            xfoil.AddRun('Run' + str(n), Airfoil, 'Run'+str(n)+'Output.txt')
            xfoil.Runs[-1].SetFlap(self.Fc, Ht, de)
            xfoil.Runs[-1].ToggleVisc(wing.Re(V=Vs[0]))
            for V in Vs:
                xfoil.Runs[-1].AddHingeMoment(alpha2d, wing.Re(V=V))
        
        if Execute:
            xfoil.ExecuteXFoil()
        
        HingeMoments = xfoil.ReadXFoilFiles()
        
        pyl.figure(fig)
        legend = []
        for i in xrange(len(des)):
            HM = self.HingeMoment( HingeMoments[i], ArmRatio, Vs )
            pyl.plot( Vs/(FT/SEC), HM/g/(IN*OZM) )
            legend.append( AsUnit(des[i], 'deg') )
            
        pyl.xlabel('V (ft/s)')
        pyl.ylabel('Torque Required (in*oz)')
        pyl.gca().yaxis.set_label_coords(-0.1,0.5)
        
        pyl.axhline(y=self.Servo.Torque / (IN*OZM), linewidth=3,  color='r', linestyle = '-.')
        ylim = pyl.ylim()
        pyl.ylim(ymin = 0, ymax = max(ylim[1]*1.1, self.Servo.Torque*1.1 / (IN*OZM)))
        
        if ShowLegend:
            pyl.legend(legend, loc='best')
        pyl.grid()

#===============================================================================
    def _CalcWeight(self):
        """ 
        Calculates a weight for the control surface
        """
        if self.dirty:
            self.Refresh()
                    
        Weight  = self.param.MassBox.Volume() * self.ControlMat.ForceDensity
        
        if self.param.wing.FullWing:
            Weight *= 2

        if self.param.wing.Symmetric:
            Weight *= 2
            
        return Weight

#===============================================================================
    def CG(self):
        """ 
        Returns the CG for the control surface
        """
        if self.dirty:
            self.Refresh()
            
        CG = self.param.MassBox.CG()
        if self.param.wing.FullWing:
            CG += self.param.MassBox2.CG()
            CG /= 2
            
        return CG
        
#===============================================================================
    def MOI(self):
        """ 
        Returns the moments of inertia for the control surface
        """
        if self.dirty:
            self.Refresh()
    
        MOI = self.param.MassBox.MOI()
        if self.param.wing.FullWing:
            MOI += self.param.MassBox2.MOI()
            
        return MOI
                        
#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """Must refresh if the wing is dirty as well"""
        if self.param.wing.dirty:
            self.Refresh()

        wing  = self.param.wing
        Servo = self.Servo
        
        Servo.Draw(fig, top, side, front)
        
#        self.param.MassBox.CopyVisibility(self)
#        self.param.MassBox.Draw(fig, top, side, front)
        
        #
        # If the control surfaces is not at the root and the wing is a full wing
        # generate a mirror copy of the servo
        #
        if wing.FullWing and self.RootFb() > 0.0:
            
            #
            # Mirror the servo appropriately
            #
            if abs(wing.Axis[0]) == 1: #Horizontal wing
                ax = 1
            else:                      #Vertical wing
                ax = 2

            Servo.X[ax] = -Servo.X[ax] + 2*wing.X[ax]*wing.Axis[ax-1]
            super(ACControlSurf.CSServo, Servo).Refresh()
            Servo.Draw(fig, top, side, front)
            
            #
            # This will put it back in the right position
            #
            Servo.Refresh() 
        
        super(ACControlSurf, self).Draw(fig, top, side, front)

#===============================================================================
    def Refresh(self):
        """
        Updates all the drawing arrays and other internal quantities.
        """
 
        super(ACControlSurf, self).Refresh()
        self.param.refreshing = True
        self.Servo.param.refreshing = True
        self._dCl_d()

        wing = self.param.wing
        
        self.CopyVisibility(wing)

        Control = self.Drawings.Control
        ControlMir = self.Drawings.ControlMir
        #
        # Update the slope of the airfoil curve
        #
        Re = self.param.wing.Re()
        self.param.dCl_dalpha = wing.af.dCl_da(Re, wing.ClSlopeAt)

        #
        # Refresh the servo
        #
        self.Servo.Symmetric = wing.Symmetric
        #self.Servo.Refresh()
               
        #
        # Fill the draw array
        #
        b = wing.b / IN
        Fc = self.Fc
        Fb = self.Fb
        Ft = self.Ft
        X = wing.X
        Chord = self._WingChord
        TE = self._WingTE
        self.Symmetric = wing.Symmetric

        b = b/2. if wing.FullWing else b

        #
        # Compute the spanwise locations
        #
        span = npy.array(
               [b*(1. - (Fb + Ft)),
                b*(1. - (Fb + Ft)),
                b*(1. - Ft),
                b*(1. - Ft),
                b*(1. - (Fb + Ft))])

        edges = npy.array(
               [TE(span[0]*IN),
                TE(span[1]*IN) - Fc * Chord(span[1]*IN),
                TE(span[2]*IN) - Fc * Chord(span[2]*IN),
                TE(span[3]*IN),
                TE(span[4]*IN)])

        if abs(wing.Axis[0]) == 1:
            # Top view
            self.ShowSide = False
            Control.yt = span*wing.Axis[0]  + X[1] / IN
            Control.xt = edges
            if wing.FullWing:
                ControlMir.yt = -Control.yt + 2*X[1]/ IN*wing.Axis[0]
                ControlMir.xt = Control.xt
        else:
            # Side view
            self.ShowTop = False
            Control.zs = span*wing.Axis[1]  + X[2]/ IN
            Control.xs = edges
            if wing.FullWing:
                ControlMir.zs = -Control.zs + 2*X[2]/ IN*wing.Axis[1]
                ControlMir.xs = Control.xs

        # No front view for the controls


        #
        # Setup a 'representative' mass box of the control surface for MOI and CG
        #
        MassBox            = ACMassBox()
        MassBox2           = ACMassBox()
        MassBox.Xat        = (0, 0.5, 0.5)
        MassBox2.Xat       = (0, 0.5, 0.5)
        MassBox.Axis       = [0,  wing.Axis[0],  wing.Axis[1]]
        MassBox2.Axis      = [0, -wing.Axis[0], -wing.Axis[1]]
        MassBox.Symmetric  = wing.Symmetric
        MassBox2.Symmetric = wing.Symmetric
        
        
        rSpan  = span[1]*IN
        tSpan  = span[2]*IN
        rChord = wing.Chord(rSpan)
        tChord = wing.Chord(tSpan)
        xr = wing.TE(rSpan) - Fc*rChord/2
        xt = wing.TE(tSpan) - Fc*tChord/2
        b  = b*IN

        L = Fb*b
        W = Fc * (rChord + tChord)/2
        H = (wing.Thickness(rSpan, xr) + wing.Thickness(tSpan, xt))/2
        
        MassBox.LWH  = [L, W, H]
        MassBox2.LWH = [L, W, H]
        
        RootFb = self.RootFb()
        mid = (wing.Mid(rSpan, xr) + wing.Mid(tSpan, xt))/2
        
        #
        # Compute the X location of the box
        #
        if abs(wing.Axis[0]) == 1:
            BoxX = [(xr + xt)/2, wing.X[1] + RootFb*b*wing.Axis[0], mid]
        else:
            BoxX = [(xr + xt)/2, mid, wing.X[2] + RootFb*b*wing.Axis[1]]

        MassBox.X = BoxX

        #
        # Mirror the second box if it is a full wing
        #
        if wing.FullWing:
            if abs(wing.Axis[0]) == 1: # Horizontal wing
                MassBox2.X = [BoxX[0], -BoxX[1] + 2*X[1]*wing.Axis[0], BoxX[2]]
            else: # Vectical wing
                MassBox2.X = [BoxX[0], BoxX[1], -BoxX[2] + 2*X[2]*wing.Axis[1]]


        self.param.MassBox = MassBox
        
        if wing.FullWing:
            self.param.MassBox2 = MassBox2
        else:
            self.param.MassBox2 = None
        #
        # The weight must be assigned now because the MassBox may be used to calculate it! Can we say catch 22?
        #
        if wing.FullWing:
            MassBox.Weight  = self.Weight/2
            MassBox2.Weight = self.Weight/2
        else:
            MassBox.Weight  = self.Weight

        self.param.refreshing = False
        self.Servo.param.refreshing = False

#===============================================================================
    def _InitKdelc(self):
        """
        Data file for non-linear correction factor (K') data.
        Taken from Fig. 9.9 of Nicolai.
        """

        # Column values
        # Control surface chord fractions
        Fc    = npy.array([    0.1,  0.15,   0.2,     0.25,    0.30,    0.4,   0.5  ])

        # Row values
        # Control deflection in degrees
        del_c = npy.array([   10.,   20.,    30.,     40.,     50.,     60.    ])

        # Cl correction factor
        kdel_c = npy.array([[   1,    0.9,    0.715,   0.65,    0.6,     0.575  ],
                            [   1,    0.89,   0.7,     0.62,    0.575,   0.54   ],
                            [   1,    0.875,  0.675,   0.6,     0.55,    0.51   ],
                            [   1,    0.85,   0.65,    0.575,   0.525,   0.4875 ],
                            [   1,    0.8,    0.6125,  0.55,    0.5,     0.4625 ],
                            [   1,    0.75,   0.575,   0.515,   0.475,   0.4375 ],
                            [   1,    0.7,    0.55,    0.49,    0.45,    0.425  ]])


        self.param.kdel_c = RectBivariateSpline(Fc, del_c, kdel_c, kx = 2, ky = 2)

#===============================================================================
    def _dCl_d(self):
        """
        Sets up the interpolation of delta Cl values associated with a control deflection
        Taken from page 9.9 of Nicolai.

        This is technically for a NACA0014 airfoil... However, the curve varies
        little between airfoils.
        """

        #if self.Fc < 0.05 or self.Fc > 0.5:
        #    raise ACControlSurfError("Control surface fraction must be between 0.05 and 0.5 for ",self.name)

        # Chord fractions of the control surfaces
        Fc     = [0.05,    0.1,    0.15,     0.2,    0.25,     0.3,     0.35,    0.4,    0.45,    0.5]
        # Delta Cl ceofficients
        dCl_da = [1,       2.6,    3.2,      3.7,    4.2,      4.65,    4.9,     5.3,    5.65,    5.9]

        # Create the interpolating function
        #dCl_d = interp1d(Fc, dCl_da, bounds_error=False, fill_value=np.nan)
        dCl_d = InterpolatedUnivariateSpline(Fc, dCl_da, k =1)

        # Calculate the delta Cl ceofficients associated with the chord fraction
        self.param.dCl_d = dCl_d(self.Fc) * 1/RAD

#===============================================================================
    def _DAlphaZero(self, del_c):
        """
        Calculates alpha zero lift offset
        Taken from page 9.9 of Nicolai.

        Inputs:
            del_c - The control surface deflectin
        """

        dCl_dalpha = self.param.dCl_dalpha

        dCl_d = self.param.dCl_d
        kdel_c = self.param.kdel_c

        Fc = self.Fc
        
        del_c = self.MakeList(del_c)
        
        da0 = []
        for dc in del_c:
            kd = kdel_c(Fc, abs(dc / ARCDEG))[0,0]
            da0.append( (-dCl_d/dCl_dalpha * dc * kd) / ARCDEG )

        if len(da0) == 1:
            return da0[0]*ARCDEG
        else:
            return npy.array(da0)*ARCDEG

#===============================================================================
    def DCl(self, del_c):
        """
        Calculates the Cl offset due to the control deflection
        Taken from page 9.9 of Nicolai.

        Inputs:
            del_c - The control surface deflection
        """
        if self.dirty or self.param.wing.dirty:
            self.Refresh()

        #Multiply by the span fraction to reduce the effictiveness of the control surface for partial span
        return -(self.param.dCl_dalpha * self._DAlphaZero(del_c)) * self.Span()/self.param.wing.b

#===============================================================================
    def DCd(self, del_c):
        """
        Calculates the Cd offset due to the control deflection

        TODO: This needs to be implemented. Hopefully this can be found in Nocolai. Near page 9.9 maybe?

        Inputs:
            del_c - The control surface deflection
        """
        if self.dirty or self.param.wing.dirty:
            self.Refresh()

        return 0

#===============================================================================
    def DCm(self, del_c):
        """
        Calculates the Cm offset due to the control deflection

        TODO: This needs to be implemented. Hopefully this can be found in Nocolai. Near page 9.9 maybe?

        Inputs:
            del_c - The control surface deflection
        """
        if self.dirty or self.param.wing.dirty:
            self.Refresh()

        return 0

################################################################################
class ACLiftSurf(ACComponent, ACAVLWing):
    """
    
    Attributes:
        Axis        - (y, z) Defining axis in the y, z plane
                      (1, 0) gives a horizontal lifting surface
                      (0, 1) gives a vertical lifting surface
        AR          - Aspect ratio
        b           - Span of the wing ( independent of FullWing)
        S           - Planform area
        i           - Incidence angle
        TR          - Array of taper ratios
        Lam         - Array of sweeps about the quarter chord for each section
        SweepFc     - Fraction of the chord to sweep about
        Gam         - Array of dihedral angles
        Fb          - Array of fraction of semi span for each section
        CEdge       - ''  , LE and TE determined by Lam
                      'LE', constant leading edge
                      'TE', constant trailing edge
        ConstUpper  - True/False constant upper surface
        FullWing    - True/False If true a full wing will be generated, otherwise only a half wing
        Airfoil     - Tuple of XFOIL coordinate file and polar files directory
        Inverted    - True/False Determines if the lifting surface is inverted
        o_eff       - Oswald efficiency
        FWCF        - Finite wing correction factor
        LLFile      - A Lifting Line Theory polar generated by XFLR5
        LLRe        - Reynolds number of the lifting line theory polar

        ClSlopeAt   - Angle of attack to calculate the Cl slope of the lifting surface
        CmSlopeAt   - Angle of attack to calculate the Cm slope of the lifting surface

        Weight      - Actual weight of the wing
        WingWeight  - An ACWingWeight class used to calculate the weight of the wing

        PointLoads  - A dictionary of point loads to modify the load distribution
    """

#===============================================================================
    def __init__(self, index):
        """
        Inputs:
            index - Index of the lifting surface for AVL
        """
        ACComponent.__init__(self)
        UnitList = self.UnitList

        # Initialize the dictionary with available properties
        self.__dict__['X']          = [0 * IN, 0* IN, 0 * IN] ;   UnitList['X']           = Length
        self.__dict__['Axis']       = (1, 0)                  ;   UnitList['Axis']        = Unitless
        self.NoneList['AR']         = None                    ;   UnitList['AR']          = Unitless
        self.NoneList['b']          = None                    ;   UnitList['b']           = Length
        self.NoneList['S']          = None                    ;   UnitList['S']           = Area
        self.__dict__['i']          = 0 * ARCDEG              ;   UnitList['i']           = Angle
        self.__dict__['CEdge']      = ''
        self.__dict__['ConstUpper'] = False
        self.__dict__['FullWing']   = False
        self.__dict__['TR']         = [1]                     ;   UnitList['TR']          = Unitless
        self.__dict__['Lam']        = [0 * ARCDEG]            ;   UnitList['Lam']         = Angle
        self.__dict__['SweepFc']    = 0.25                    ;   UnitList['SweepFc']     = Unitless
        self.__dict__['Gam']        = [0 * ARCDEG]            ;   UnitList['Gam']         = Angle
        self.__dict__['Fb']         = [1]                     ;   UnitList['Fb']          = Unitless
        self.__dict__['Airfoil']    = ''
        self.__dict__['Inverted']   = False
        self.NoneList['o_eff']      = None                    ;   UnitList['o_eff']       = Unitless
        self.__dict__['FWCF']       = 0.98                    ;   UnitList['FWCF']        = Unitless
        self.NoneList['LLFile']     = None
        self.NoneList['LLRe']       = None                    ;   UnitList['LLRe']        = Unitless
        self.NoneList['ClSlopeAt']  = None                    ;   UnitList['ClSlopeAt']   = Angle
        self.NoneList['CmSlopeAt']  = None                    ;   UnitList['CmSlopeAt']   = Angle
        self.NoneList['Weight']     = None                    ;   UnitList['Weight']      = Force
        self.NoneList['WingWeight'] = None

        self.__dict__['PointLoads'] = {} # A dictionary for storing point loads

        self.__dict__['WeightGroup'] = 'None'

        self.__dict__['LLT'] = ACLiftingLine(self)
        self.LLT.Alpha3d = 5*ARCDEG
        
        #
        # Objects available for drawing
        #
        self.Drawings.LiftSurf = ACComponent.ComponentDrawing()
        self.Drawings.LiftSurf.color = 'r'

        self.Drawings.LiftSurfMir = ACComponent.ComponentDrawing()
        self.Drawings.LiftSurfMir.color = 'r'

        #
        # Root and tip airfoils
        #
        self.Drawings.AirfoilR = ACComponent.ComponentDrawing()
        self.Drawings.AirfoilR.color = 'r'

        self.Drawings.AirfoilT = ACComponent.ComponentDrawing()
        self.Drawings.AirfoilT.color = 'r'

        #
        # Internal lists with inAsUnition stored in units of inches
        #
        self.param.chord      = [] #Chord lengths
        self.param.span       = [] #Span locations of chords
        self.param.sweep      = [] #Offset added to each span location for sweep
        self.param.dihed      = [] #Dihedral offset
 

        #
        # Internal variables to reduce calculations. These are never set
        # but it is faster to compute and save them rather than recompute
        # every time they are needed
        #
        self.param.MAC   = None  #Mean aerodynamic chord
        self.param.index = index #AVL index of this lifting surface

        #
        # An class for storing control surfaces
        #
        class Controls:
            pass
        self.__dict__['Controls'] = Controls()

        #
        # An class for storing winglets
        #
        class Winglets:
             pass
        self.__dict__['Winglets'] = Winglets()

        #
        # The actual airfoil object of the wing
        #
        self.__dict__['af']       = None

#===============================================================================
    def AddControl(self, ControlName):
        """
        Adds a control surface to the lifting surface

        Inputs:
            ControlName - The name of the control surface
        """
        self.Controls.__dict__[ControlName] = ACControlSurf(ControlName, self)
        self.__dict__[ControlName] = self.Controls.__dict__[ControlName]
        self.__dict__[ControlName].Parents.append( self )
        return self.__dict__[ControlName]

#===============================================================================
    def AddPointLoad(self, y, Load):
        """
        Adds a point load to modify the load distribution on the wing

        Inputs:
            y    - Spanwise location of the load
            Load - The pointwise load
        """
        self.PointLoads[y] = Load

#===============================================================================
    def SetWeightCalc(self, WingWeight):
        """
        Adds a weight calculator to the lifting surface

        Inputs:
            WingWeight - A wing weight object
        """
        self.WingWeight = WingWeight(self)
        #self.WingWeight.Parents.append( self )
        self.Parents.append( self.WingWeight )
        return self.WingWeight

#===============================================================================
    def GetWeightCalc(self):
        """
        Retreives the weight calculator to the lifting surface
        """
        return self.WingWeight

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the lift surface and its parts to the Weight Table
        
        Input:
            PartName    - The name give to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WeightTable[PartName] = self
        
        if self.WingWeight is not None:
            self.WingWeight.AddToWeightTable(self.name, WeightTable[PartName])
        
        self._AddSubPartsWeightTable(self.Winglets, WeightTable[PartName])
        self._AddSubPartsWeightTable(self.Controls, WeightTable[PartName])
        

#===============================================================================
    def _AddSubPartsWeightTable(self, SubPartDict, WeightTable):
        """
        Adds sup-parts to the weight table
        """  
        for key, SubPart in SubPartDict.__dict__.iteritems():
            SubPart.AddToWeightTable(key, WeightTable)


#===============================================================================
    def AddWinglet(self, Wingletname, index):
        """
        Adds a winglet to the lifting surface

        Inputs:
            Wingletname - The name of the winglet
            index       - The avl index for the winglet lifting surface
        """
        self.Winglets.__dict__[Wingletname] = ACWinglet(Wingletname, index, self)
        self.Winglets.__dict__[Wingletname].Parents.append( self )
        return self.Winglets.__dict__[Wingletname]

#===============================================================================
#
# Interface methods to retrieve inAsUnition about the lifting surface
#
#===============================================================================
    def Chord(self, y):
        """
        Calculates the chord of the lifting surface given
        a spanwise location. Negative values are acceptable
        only if Sym is true.

        Inputs:
            y - The wing span location
        """
        #
        # Update chord and span locations if needed
        #
        if self.dirty:
            self.Refresh()

        #
        # Make sure it is a valid spanwise location
        #
        y = self._CheckSpanLoc(y)

        #
        # Interpolate the chord length for the desired position
        #
        ay = npy.abs(y / IN)
        return npy.interp(ay, self.param.span, self.param.chord) * IN

#===============================================================================
    def LE(self, y):
        """
        Returns the leading edge given a spanwise location

        Inputs:
            y - span location
        """
        #
        # 'span' and 'sweep' must exist and have units of inches
        #
        if self.dirty:
            self.Refresh()

        #
        # Get the chord
        #
        C = self.Chord(y) / IN
        yi = y / IN
        
        #
        # The Chord method validated the spanwsie position. But it must still be flipped if it is negative
        #
        if isinstance( yi, npy.ndarray):
            if npy.all(yi < 0):
                y = -y       
        else:
            if yi < 0:
                y = -y       
        #
        # Interpolate the sweep offset
        #
        sweep = npy.interp(y / IN, self.param.span, self.param.sweep);

        # Sweep about the desired chord, but make sure that the root chuarter chord
        # is still the reference location
        return (-0.25*C + sweep) * IN + self.X[0]

#===============================================================================
    def TE(self, y):
        """
        Returns the trailing edge given a spanwise location

        Inputs:
            y - span location
        """
        # This one is thankfully very simple
        return self.LE(y) + self.Chord(y)

#===============================================================================
    def MaxTE(self):
        """
        Returns maximum trailing edge location
        """
        if self.dirty: self.Refresh()
        
        Spans = self.param.span*IN
        MaxTE = 0*IN
        for span in Spans:
            MaxTE = max(MaxTE, self.TE(span))
        
        return MaxTE
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

        #
        # 'span' and 'dihed' must exist and have units of inches
        #
        if self.dirty:
            self.Refresh()
            
        C  = self.Chord(y)     # Get the chord
        
        if x is None:
            af_y = self.af.Lower()    # Min value of the airfoil thickness
        else:
            ax = (x - self.LE(y)) / C
            if ax < 0.0 or ax > 1.0:
                message = "Chord location out of bounds\n" + \
                          "LE: " + str(self.LE(y)) + " x: " + str(x) + " TE: " + str(self.TE(y)) + "\n in "
                raise ACLiftSurfError(message, self.name)
            af_y = self.af.Lower(ax) #Get the airfoil lower value at this location

        #
        # Interpolate the dihedral offset
        #
        ay  = y / IN
        mid = npy.interp(ay, self.param.span, self.param.dihed);
        
        if abs(self.Axis[0]) == 1: 
            yz = self.X[2] #Horizontal wing
        else:
            yz = self.X[1] #Vertical wing

        # Calculate the bottom accounting for the thickness of the airfoil
        return mid * IN + af_y * C + yz

#===============================================================================
    def Upper(self, y, x = None):
        """
        Returns the upper side given a spanwise location

        Inputs:
            y - span location
            x - chord location (Max y value if None)
                Note that the airfoil max and min y values may
                not correspond to the same x location
        """
        #
        # 'span' and 'dihed' must exist and have units of inches
        #
        if self.dirty:
            self.Refresh()

        C = self.Chord(y)        # Get the chord

        if x is None:
            af_y = self.af.Upper()  # Maximum height of the airfoil
        else:
            ax = (x - self.LE(y)) / C 
            if ax < 0.0 or ax > 1.0:
                message = "Chord location out of bounds\n" + \
                          "LE: " + str(self.LE(y)) + " x: " + str(x) + " TE: " + str(self.TE(y)) + "\n in "
                raise ACLiftSurfError(message, self.name)
            af_y = self.af.Upper(ax) #Get the airfoil upper value at this location

        #
        # Interpolate the dihedral offset
        #
        ay = y / IN
        mid = npy.interp(ay, self.param.span, self.param.dihed);

        if abs(self.Axis[0]) == 1: 
            yz = self.X[2] #Horizontal wing
        else:
            yz = self.X[1] #Vertical wing

        # Calculate the top accounting for the thickness of the airfoil
        return mid * IN + af_y * C + yz

#===============================================================================
    def Thickness(self, y, x = None):
        """
        Returns the thickness given a spanwise location

        Inputs:
            y - span location
            x - chord location (Default = Airfoil TC)
        """
        if self.dirty:
            self.Refresh()
        
        if x is None:
            return self.af.TC * self.Chord(y)
        else:
            return self.Upper(y, x) - self.Lower(y, x)

#===============================================================================
    def Mid(self, y, x = None):
        """
        Returns the mid point given a spanwise location

        Inputs:
            y - span location
            x - chord location
        """
        if x is not None:
            return (self.Upper(y, x) + self.Lower(y, x))/2
        
        if self.dirty: self.Refresh()
        #
        # Interpolate the dihedral offset to get the mid point
        #
        ay  = y / IN
        mid = npy.interp(ay, self.param.span, self.param.dihed);

        if abs(self.Axis[0]) == 1: 
            yz = self.X[2] #Horizontal wing
        else:
            yz = self.X[1] #Vertical wing

        return mid * IN + yz

#===============================================================================
    def Tip(self):
        """
        Returns the quarter chord of the tip of the lifting surface
        """
        b = self.b
        b = b/2. if self.FullWing else b
        Ct = self.Chord(b)
        
        Tip = [0,0,0]
        if abs(self.Axis[0]) == 1: # Horizontal Wing
            Tip[0] = self.LE(b) + 0.25*Ct
            Tip[1] = self.X[1] + b * self.Axis[0]
            Tip[2] = self.Mid(b)
        else:                      # Vertical Wing
            Tip[0] = self.LE(b) + 0.25*Ct
            Tip[1] = self.Mid(b)
            Tip[2] = self.X[2] + b*self.Axis[1]
            
        return Tip
    
#===============================================================================
    def Volume(self):
        """
        Calculates the volume of the wing
        """
        chord  = self.param.chord * IN
        span   = self.param.span * IN
        AFArea = self.af.Area

        Vol = 0*IN**3
        for i in range(len(chord)-1):
            Vol += (chord[i]**2 + chord[i+1]**2)*AFArea/2 * (span[i+1]-span[i])

        return 2*Vol if self.FullWing else Vol

#===============================================================================
    def SWet(self):
        """
        Calculates an appriximate wetted area
        """
        #
        # 'chord' must exist and have units of inches
        #
        if self.dirty:
            self.Refresh()

        chord = self.param.chord
        t = 0 * IN #Average thickness of the wing
        for i in range(len(chord)):
            t += self.af.TC * chord[i] * IN
        t *= 1./float(len(chord))

        S = self.S
        b = self.b

        return 2. * S + t * b #Calculate the approximated wetted area

#===============================================================================
    def MAC(self):
        """
        Retrieves the mean aerodynamic chord
        """
        if self.dirty:
            self.Refresh()

        if self.param.MAC is None:
            self._MAC()

        return self.param.MAC

#===============================================================================
    def _MAC(self):
        """
        Calculates the mean aerodynamic chord

        The equation is given by
        mac=\frac{2}{S}\int_{0}^{b/2}c\left(y\right)^{2}dy
        =\frac{2}{S}\left[\frac{a^{2}y^{3}}{3}+aby^{2}+b^{2}y\right]_{0}^{b/2}
        """
        #
        # 'chord' and 'span' must exist and have units of inches
        #

        S = self.S
        chord = self.param.chord
        span = self.param.span
        mac = 0

        for i in range(len(chord)-1):
            a = (chord[i+1] - chord[i])/(span[i+1]-span[i])
            b = chord[i]
            y = span[i+1]-span[i]
            mac += a**2 * y**3 / 3. + a * b * y**2 + b**2 * y

        fac = 2. if self.FullWing else 1.

        #
        # IN**3 becase of integral of C**2 (IN**2)
        #
        self.param.MAC = mac * IN**3 * fac/S

#===============================================================================
    def DownWash(self, alpha2d, Re = None):
        """
        Calculates the shed vortex induced downwas angle assuming an elliptocal loading.

        Ref: Abbott and Doenhoff pg 7, Prandtl wing theory

        Inputs:
            alhpah2d - Angle of attack of the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        return self.af.Cl(alpha2d, Re)/(math.pi * self.o_eff * self.AR) * RAD

#===============================================================================
    def CLmax(self, Re = None):
        """
        Calculates max CL with the finite wing correction factor

        Inputs:
            Re - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        return self.af.Clmax(Re) * self.FWCF

#===============================================================================
    def CLmin(self, Re = None):
        """
        Calculates min CL with the finite wing correction factor

        Inputs:
            Re - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        return self.af.Clmin(Re) * self.FWCF

#===============================================================================
    def Cl(self, alpha2d, y, Re):
        """
        Calculates local airfoil Cl

        Inputs:
            alhpah2d - Angle of attack the airfoil
            y        - Spanwise location
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        #Right now y is not used as a single airfoil is assumed....
        return self.af.Cl(alpha2d, Re)

#===============================================================================
    def CL(self, alpha2d, Re = None):
        """
        Calculates CL with the finite wing correction factor

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        S = self.S
        CL = self.af.Cl(alpha2d, Re) * self.FWCF
        CLWinglet = self._WingletCL(self.Winglets, alpha2d)

        return CL + (CLWinglet/S)

#===============================================================================
    def _WingletCL(self, Winglets, alpha2d):
        """
        Calculates area weighted CL of the winglets
        
        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        try:
            CL = npy.zeros(len(alpha2d))*IN**2
        except:
            CL = 0*IN**2
            
        for Winglet in Winglets.__dict__.itervalues():
            S    = Winglet.S if not Winglet.Symmetric else 2*Winglet.S 
            Axis = Winglet.Axis
            i    = Winglet.i
            CL  += Winglet.CL(alpha2d + i)*S*abs(Axis[0])

        return CL
#===============================================================================
    def CDi(self, alpha2d, Re = None):
        """
        Calculates vortex induced drag

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()
        
        #
        # Must compute the CL of the wing like this. Otherwise 
        # the induced drag will use the lift from all winglets, which is not correct.
        #
        CL = self.af.Cl(alpha2d, Re) * self.FWCF
        
        return CL**2/(math.pi * self.o_eff * self.AR)

#===============================================================================
    def Cd(self, alpha2d, y, Re):
        """
        Calculates local airfoil Cd

        Inputs:
            alhpah2d - Angle of attack the airfoil
            y        - Spanwise location
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        #Right now y is not used as a single airfoil is assumed....
        return self.af.Cd(alpha2d, Re)

#===============================================================================
    def CD(self, alpha2d, Re = None):
        """
        Calculates CD including the induced drag
        
        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        #
        # Compute an area wegithed drag coefficient between the wing
        # and winglets
        # 
        S = self.S
        CD = self.af.Cd(alpha2d, Re) + ACLiftSurf.CDi(self, alpha2d, Re)
        CDWinglet = self._WingletCD(self.Winglets, alpha2d)
        
        return CD + (CDWinglet/S)

#===============================================================================
    def _WingletCD(self, Winglets, alpha2d):
        """
        Calculates area weighted CD of the winglets
        
        Inputs:
            alhpah2d - Angle of attack the wing
        """
        try:
            CD = npy.zeros(len(alpha2d))*IN**2
        except:
            CD = 0*IN**2
            
        for Winglet in Winglets.__dict__.itervalues():
            S    = Winglet.S if not Winglet.Symmetric else 2*Winglet.S
            Axis = Winglet.Axis
            i    = Winglet.i
            CD  += Winglet.CD((alpha2d + i)*abs(Axis[0]))*S

        return CD
#===============================================================================
    def dCL_da(self, alpha2d = None, Re = None):
        """
        Calculates the derivative dCL/da with the finite wing correction factor

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        if alpha2d is None:
            alpha2d = self.ClSlopeAt

        return self.af.dCl_da(Re, alpha2d) * self.FWCF

#===============================================================================
    def dCDi_da(self, alpha2d = None, Re = None):
        """
        Calculates induced drag

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        return 2.*self.dCl_da(Re, alpha2d)/(math.pi * self.o_eff * self.AR)

#===============================================================================
    def dCD_da(self, alpha2d = None, Re = None):
        """
        Calculates the derivative dCD/da

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()

        return self.af.dCd_da(Re, alpha2d) + self.dCDi_da(alpha2d, Re)

#===============================================================================
    def CM(self, alpha2d, Re = None):
        """
        Calculates CM

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()
            
        Sw   = self.S
        macw = self.MAC()
        Xacw = self.Xac()
        
        CM  = self.af.Cm(alpha2d, Re)
        CM += self._WingletCM(self.Winglets, alpha2d, Sw, macw, Xacw)
        
        return CM

#===============================================================================
    def _WingletCM(self, Winglets, alpha2d, Sw, macw, Xacw):
        """
        Calculates area weighted CM of the winglets
        
        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        try:
            CM = npy.zeros(len(alpha2d))
        except:
            CM = 0
            
        for Winglet in Winglets.__dict__.itervalues():
            S    = Winglet.S if not Winglet.Symmetric else 2*Winglet.S
            Axis = Winglet.Axis
            i    = Winglet.i
            Re   = Winglet.Re()
            CM  += Winglet.af.Cm(alpha2d + i, Re)*abs(Axis[0])*S/Sw
            CM  += ((Xacw - Winglet.Xac())*S/(Sw*macw))*abs(Axis[0])
            
            CM  += Winglet._WingletCM(Winglet.Winglets, alpha2d, Sw, macw, Xacw)

        return CM
#===============================================================================
    def dCM_da(self, alpha2d = None, Re = None):
        """
        Calculates the derivative CM w.r.t angle of attack

        Inputs:
            alhpah2d - Angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        if Re is None:
            Re = self.Re()
        
        if alpha2d is None:
            alpha2d = self.CmSlopeAt

        return self.af.dCm_da(Re, alpha2d)

#===============================================================================
    def AlphaRange(self, nalpha = 30):
        """
        Returns an appropriate range of angles

        Inputs:
            nalpha - Number of angles of attack
        """
        if self.dirty: self.Refresh()

        return self.af.AlphaRange(nalpha)

#===============================================================================
    def AlphaClmax(self, Re = None):
        """
        Returns the angle of attack with the maximum Cl

        Inputs:
            Re - Reynolds number
        """
        if self.dirty: self.Refresh()
        
        if Re is None:
            Re = self.Re()

        return self.af.AlphaClmax(Re)

#===============================================================================
    def AlphaClmin(self, Re = None):
        """
        Returns the angle of attack with the minimum Cl

        Inputs:
            Re - Reynolds number
        """
        if self.dirty: self.Refresh()
        
        if Re is None:
            Re = self.Re()

        return self.af.AlphaClmin(Re)

#===============================================================================
    def AlphaWing(self, alpha2d, Re = None):
        """
        Computes the three dimensional wing angle of attack for plotting purposes
        
        Inputs:
            alhpah2d - 2D angle of attack the wing
            Re       - Reynolds number
        """
        if self.dirty: self.Refresh()

        return alpha2d + self.DownWash(alpha2d, Re)

#===============================================================================
    def CG(self):
        """
        Calculates center of gravity of the lifting surface

        Defined to be about the quarter chord of the MAC
        """
        if self.dirty: self.Refresh()
        
        if self.WingWeight is not None:
            Weight = self.WingWeight.Weight
            CG     = self.WingWeight.CG()*Weight
            for Winglet in self.Winglets.__dict__.itervalues():
                CG     += Winglet.CG()*Winglet.Weight
                Weight += Winglet.Weight

            for Control in self.Controls.__dict__.itervalues():
                CG     += Control.CG()*Control.Weight
                Weight += Control.Weight
                
                Servo   = Control.Servo
                tmpCG   = Servo.CG()
                CG     += tmpCG*Servo.Weight
                Weight += Servo.Weight
                if self.FullWing and Control.RootFb > 0.0:
                    tmpCG   = Servo.Mirror(tmpCG)
                    CG     += tmpCG*Servo.Weight
                    Weight += Servo.Weight
            
            return CG/Weight
        
        else:
            if self.FullWing:
                CG = self.GetX()
            else:
                x = (self.MAC()*0.25) / IN
                b  = self.b / IN
                if abs(self.Axis[0]) == 1: #Horizontal surface
                    b *= self.Axis[0]
                    CG = npy.array([x, b/2, 0])*IN
                else:                      #Vertical surface
                    b *= self.Axis[1]
                    CG = npy.array([x, 0, b/2])*IN
                
                CG += self.GetX()
    
            return CG

 #===============================================================================
    def MOI(self):
        """
        Calculates the moments of inertia for the lifting surface
        """

        if self.WingWeight is not None:
            MOI = self.WingWeight.MOI()
        else:
            MOI = npy.array([0,0,0]) * LBM*IN**2
        
        CG = self.CG()
        
        def AddMOI(PCG, PM, PMOI):
            dCG = PCG - CG
            MOI[0] += PMOI[0] + PM*(dCG[1]**2 + dCG[2]**2) #Ixx
            MOI[1] += PMOI[1] + PM*(dCG[0]**2 + dCG[2]**2) #Iyy
            MOI[2] += PMOI[2] + PM*(dCG[1]**2 + dCG[0]**2) #Izz

        #
        # Add moments of inertia for the control surfaces and their servos
        #
        for Control in self.Controls.__dict__.itervalues():
            
            AddMOI(Control.CG(), Control.Weight / g, Control.MOI())
            if self.FullWing:
                AddMOI(Control.CG(), Control.Weight / g, Control.MOI())
                
            Servo = Control.Servo
            AddMOI(Servo.CG(), Servo.Weight / g, Servo.MOI())
            if self.FullWing and Control.RootFb > 0.0:
                AddMOI(Servo.CG(), Servo.Weight / g, Servo.MOI())

        #
        # Add moments of inertia for the winglets
        # This will recursevely add up the MOI for winglets on winglets
        #
        for Winglet in self.Winglets.__dict__.itervalues():
            AddMOI(Winglet.CG(), Winglet.Weight / g, Winglet.MOI()) 

        return MOI

#===============================================================================
    def Xac(self, Re = None):
        """
        Calculates the location of the aerodynamic center

        This assumes that Cm was calculated about the quarter
        chord of the airfoil.

        Inputs:
            alpha2d - Angle of attack used to calculated alpha derivatives
            Re      - Reynolds number
        """
        if self.dirty:
            self.Refresh()

        if Re is None:
            Re = self.Re()

        alphaCl = self.ClSlopeAt
        alphaCm = self.CmSlopeAt

        #
        # Aerodynamic center is relative the Quarter chord because thats where XFOIL 
        # copmutes the moments about
        #
        return self.X[0] - self.MAC()*self.af.dCm_da(Re, alphaCm) / self.af.dCl_da(Re, alphaCl)

#===============================================================================
    def EllipChord(self, y):
        " Calculats the chord to generate an elliptic chord distribution "
        #
        # Update chord if needed
        #
        if self.dirty:  self.Refresh()

        #
        # Make sure it is a valid spanwise location
        #
        self._CheckSpanLoc(y)

        Cr = self.param.chord[0] * IN
        b = self.b/2 if self.FullWing else self.b
        
        return Cr * npy.sqrt(1 - (y / b)**2).real


#===============================================================================
    def ChordSchreiner(self, y):
        """
        Calculates the Schreiner chord distribution used to approximate a 
        non-elliptic lift distribution
        
        Inputs:
            y - Span location            
        """
        Cw = self.Chord(y)
        Cellip = self.EllipChord(y)
        return (Cw + Cellip) / 2.

#===============================================================================
    def LoadDistSchreiner(self, y, V):
        """
        Calculates the Schreiner approximated load at a span location at CLmax

        Inputs:
            V - Velocity
            y - Span location            
        """   
        Alt_LO         = self.GetAlt_LO()
        q              = AeroUtil.q(Alt_LO, V)
        CLmax          = self.CLmax()
        ChordSchreiner = []
        
        for ay in y:
            ChordSchreiner.append(self.ChordSchreiner(ay) / IN)
        
        ChordSchreiner = ChordSchreiner[0]*IN if len(ChordSchreiner) == 1 else npy.array(ChordSchreiner)*IN
        
        return q * CLmax * ChordSchreiner
         
 #===============================================================================
    def ShearLoading(self, y, V, LoadDist = None):
        """ 
        Calculats the shear loading distribution as function of y (span) 
        
        Inputs:
            y - Span location(s)
            V - Veclocity
        """
        #
        # Update chord if needed
        #
        if self.dirty:  self.Refresh()
        
        b = self.b/2 if self.FullWing else self.b
        b = b / IN

        #
        # This assumes that ShearLoad and y are lists
        #
        if LoadDist is None:
            LoadDist = self.LoadDistSchreiner(y, V)
        
        LoadDistInt = interp1d(y / IN, LoadDist / (LBF/IN), kind='linear')

        def Load(yloc):
            return LoadDistInt(yloc)
        
        ShearLoads = []
        for ay in y:
            ShearLoads.append(integrate.romberg(Load, ay / IN, b, tol=1e-2))
        
        y = y / IN
        
        yl     = self.ToUnumList(self.PointLoads.keys(), IN)
        pLoads = self.ToUnumList(self.PointLoads.values(), LBF)
        
        yi     = yl / IN
        pLoads = pLoads / LBF
        
        for i in xrange(len(yi)):
            pL               = pLoads[i]
            pt               = y.searchsorted(yi, side='right')
            ShearLoads[:pt] += pL 
         
        return ShearLoads[0]*LBF if len(ShearLoads)==1 else npy.array(ShearLoads)*LBF
               
#===============================================================================
    def MomentLoading(self, y = None, V = None, ShearLoad = None):
        """
        Calculats the moment loading distribution as function of y (span) 
        
        Inputs:
            y - Span location(s)
            V - Veclocity
        """
        #
        # Update chord if needed
        #
        if self.dirty:  self.Refresh()

        b = self.b/2 if self.FullWing else self.b
        b = b / IN

        if y is None:
            y = npy.linspace(0.0, b, 40)*IN
        
        if V is None:
            raise ACLiftSurfError("Must specify a velocity for MomentLoading")
        #
        # This assumes that ShearLoad and y are lists
        #
        if ShearLoad is None:
            ShearLoad = self.ShearLoading(y, V)
        
        ShearInt = interp1d(y / IN, ShearLoad / LBF, kind='linear')
        
        def Shear(yloc):
            return ShearInt(yloc)
        
        Moment = []
        for ay in y:
            Moment.append(integrate.romberg(Shear, ay / IN, b, tol=1e-2))
        
        return Moment[0]*LBF*IN if len(Moment)==1 else npy.array(Moment)*LBF*IN
       
#===============================================================================
    def PlotSpanLoading(self, V, LdMult = 1, fig = 1):
        """
        Plots distributions of the spanwise loading, shear, and moment

        Inputs:
            V   - Velocity to compute loads
            fig - Figure Number
            LdMult - Generic Load multiplier for either safety factor or GLoads
        """
        if self.dirty: self.Refresh()
        
        b = self.b/2. if self.FullWing else self.b
        
        y = npy.linspace(0.0, b / IN, 40)*IN
        Load   = self.LoadDistSchreiner(y, V) * LdMult
        Shear  = self.ShearLoading(y, V, Load)
        Moment = self.MomentLoading(y, V, ShearLoad = Shear) / (LBF*IN)
        
        Load  = Load / (LBF/IN)
        Shear = Shear / (LBF)
        y = y / IN
        
        gridarg = {'color':'k', 'linestyle':'--', 'linewidth':1}
        
        pyl.figure(fig)
        ax1 = pyl.subplot(311); pyl.grid(gridarg)
        pyl.setp( ax1.get_xticklabels(), visible=False)
        pyl.title("Spanwise Loading for " + self.name)
        pyl.plot(y, Load);   pyl.ylabel("Load (lbf/in)"); #pyl.xlabel("Semi Span (in)");
        ax2 = pyl.subplot(312, sharex=ax1); pyl.grid(gridarg)
        pyl.setp( ax2.get_xticklabels(), visible=False)
        pyl.plot(y, Shear);  pyl.ylabel("Shear Load (lbf)"); #pyl.xlabel("Semi Span (in)"); 
        pyl.subplot(313, sharex=ax1); pyl.grid(gridarg)
        pyl.plot(y, Moment); pyl.ylabel("Moment (lbf in)"); pyl.xlabel("Semi Span (in)"); 
    
#===============================================================================
    def PlotNormalizedSpanLoading(self, V, LdMult = 1, fig = 1):
        """
        Plots normalized distributions of the spanwise loading, shear, and moment

        Inputs:
            V   - Velocity to compute loads
            fig - Figure Number
            LdMult - Generic Load multiplier for either safety factor or GLoads
        """
        if self.dirty: self.Refresh()
        
        b = self.b/2. if self.FullWing else self.b
        
        y = npy.linspace(0.0, b / IN, 40)*IN
        Load   = self.LoadDistSchreiner(y, V) * LdMult
        Shear  = self.ShearLoading(y, V, Load)
        Moment = self.MomentLoading(y, V, ShearLoad = Shear) / (LBF*IN)
        
        Load  = Load / (LBF/IN)
        Shear = Shear / LBF
        
        MaxLoad   = max(Load)
        MaxShear  = max(Shear)
        MaxMoment = max(Moment)
        
        y = y / IN
        
        pyl.figure(fig)
        pyl.title("Normalized Spanwise Loading for " + self.name)
        pyl.plot(y, Load/MaxLoad);   
        pyl.plot(y, Shear/MaxShear);  
        pyl.plot(y, Moment/MaxMoment); 
        pyl.xlabel("Semi Span (in)"); pyl.ylabel("Loads")
        
        AsUnit = '%1.2f'
        
        pyl.legend([('Load '   + AsUnit + ' (lbf/in)') % MaxLoad, 
                    ('Shear '  + AsUnit + ' (lbf)   ') % MaxShear, 
                    ('Moment ' + AsUnit + ' (lbf in)') % MaxMoment], loc = 'best')
                            
#===============================================================================
    def LLTPlots(self, fig = 1):
        """
        Plots the finite wing coefficients along with a Lifting Line Calculation
        """
        if self.dirty: self.Refresh()

        LLFile = self.LLFile

        f = open(LLFile,'r')
        #
        # Remove the eader from the file
        #
        line = f.readline()
        i = 0
        #
        # Find the header line before the data
        #
        while line.find('_________') == -1:
            line = f.readline()
            if len(line) == 0:
                raise ACLiftSurfError(self.LLFile + " does not appear to be a valid lifting line polar file.\n Error occured in", self.name)

        #
        # Read in the polar
        #
        LLalpha = []
        LLCl = []
        LLCd = []
        LLCm = []
        for line in f:
            vals = line.lstrip()
            vals = vals.split('  ')
            if len(vals) == 1:
                break #Just in case there is a newline at the end of the file

            LLalpha.append(float(vals[0]))
            LLCl.append(float(vals[1]))
            LLCd.append(float(vals[4]))
            LLCm.append(float(vals[5]))

        f.close()

        #
        # It is faster to compute Re once up front rather than 3 times
        # in the respective coefficient functions
        #
        Re = self.LLRe if self.LLRe is not None else self.Re()

        alphas = self.AlphaRange()
        CL = self.CL(alphas, Re)
        CD = self.CD(alphas, Re)
        CM = self.CM(alphas, Re)

        alphas = (alphas + self.DownWash(alphas)) / ARCDEG

        pyl.figure(fig)

        #
        # Plot all the polars to compare with the LLT data
        #
        pyl.subplot(221)
        pyl.plot(alphas, CL)
        pyl.plot(LLalpha, LLCl, 'r+')
        pyl.xlabel(r'$\alpha$ $[^o]$'); pyl.ylabel('$C_L$')
        pyl.legend(('Computed','LLT'),loc='upper left')

        pyl.subplot(222)
        pyl.plot(alphas, CD)
        pyl.plot(LLalpha, LLCd, 'r+')
        pyl.xlabel(r'$\alpha$ $[^o]$'); pyl.ylabel('$C_D$')
        pyl.legend(('Computed','LLT'),loc='upper center')

        pyl.subplot(223)
        pyl.plot(alphas, CM)
        pyl.plot(LLalpha, LLCm, 'r+')
        pyl.xlabel(r'$\alpha$ $[^o]$'); pyl.ylabel('$C_M$')
        pyl.legend(('Computed','LLT'),loc='upper right')

        pyl.subplot(224)
        pyl.plot(CD, CL)
        pyl.plot(LLCd, LLCl, 'r+')
        pyl.xlabel('Cd'); pyl.ylabel('Cl')
        pyl.legend(('Computed','LLT'),loc='center right')

#===============================================================================
    def Draw2DAirfoilPolars(self, fig = 1):
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
        alpha = self.af.AlphaRange()

        #
        # Create the polars
        #
        Re = self.Re()
        Cl = self.af.Cl(alpha, Re)
        Cd = self.af.Cd(alpha, Re)
        Cm = self.af.Cm(alpha, Re)
        
        ClSlopeAt = npy.array(self.ClSlopeAt)    
        CmSlopeAt = npy.array(self.CmSlopeAt)

        dCl_da = self.af.dCl_da(Re, ClSlopeAt)
        dCm_da = self.af.dCm_da(Re, CmSlopeAt)
        a0     = self.af.AlphaZeroLift(Re)
        
        #
        # Make the slope go through the specified ClSlopeAt. Otherwise let it go
        # though the Cl = 0 point.
        #
        ClOffset  = self.af.Cl(ClSlopeAt, Re) - dCl_da*ClSlopeAt \
                                         if self.ClSlopeAt is not None else -dCl_da*a0
 
        CmOffset  = self.af.Cm(CmSlopeAt, Re) - dCm_da*CmSlopeAt \
                                         if self.CmSlopeAt is not None else -dCm_da*a0
        #
        # These are in degrees without being a scalar
        #
        alp_da = npy.array([self.af.AlphaMin*ARCDEG, self.af.AlphaMax*ARCDEG])
        
        Cl_da  = alp_da*dCl_da + ClOffset
        Cm_da  = alp_da*dCm_da + CmOffset
        
        alp_da /= ARCDEG

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
        alpha = self.af.AlphaRange()

        #
        # Create the polars. Call ACLiftSurf other wise Tail surfaces will not be properly generated.
        #
        Re = self.Re()
        CL = ACLiftSurf.CL(self, alpha, Re = Re)
        CD = ACLiftSurf.CD(self, alpha, Re = Re)
        CM = ACLiftSurf.CM(self, alpha, Re = Re)
                
        alpha  = self.AlphaWing(alpha)
        
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
        
        alpha = alpha / ARCDEG
        
        #
        # Draw them
        #
        pyl.subplot(321)
        pyl.plot(alpha, CL)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_l$')
        #
        # Also draw the line prepresentative of the CL slope
        #
        if Cl_da is not None:
            pyl.plot(alp_da, Cl_da)
            pyl.legend([r'$C_l$',r'$C_l$ Slope'], loc = 'best')


        pyl.subplot(322)
        pyl.plot(alpha, CD)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_d$')

        pyl.subplot(323)
        pyl.plot(alpha, CM)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_m$')
        #
        # Also draw the line prepresentative of the CM slope
        #
        if Cm_da is not None:
            pyl.plot(alp_da, Cm_da)
            pyl.legend([r'$C_m$',r'$C_m$ Slope'], loc = 'best')

        pyl.subplot(324)
        pyl.plot(CD, CL)
        pyl.xlabel(r'$C_d$'); pyl.ylabel(r'$C_l$')

        pyl.subplot(325)
        pyl.plot(alpha,CL/CD)
        pyl.xlabel(r'$\alpha$'); pyl.ylabel('L/D')

        title += ' Re = %1.0f' % self.Re()
        pyl.annotate(title, xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)                

#===============================================================================
    def PlotDragBuildup(self, fig = 1):
        """
        Draws how the drag buildup is calculated for the lifting surface

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        #
        # Create the angles of attack
        #
        alpha = self.af.AlphaRange()

        #
        # Create the polars. Call ACLiftSurf other wise Tail surfaces will not be properly generated.
        #
        Re = self.Re()
        CL = ACLiftSurf.CL(self, alpha, Re = Re)
        CDi = ACLiftSurf.CDi(self, alpha, Re = Re)
        CDp = self.af.Cd(alpha, Re)
                
        alpha  = self.AlphaWing(alpha)
        
        gridarg = {'color':'k', 'linestyle':'--', 'linewidth':1}
        
        legend = []
        
        pyl.plot(CDi, CL); legend.append("Induced")
        pyl.plot(CDp, CL); legend.append("Profile")
        pyl.plot(CDp+CDi, CL); legend.append("Total")
        pyl.xlabel(r'$C_D$'); pyl.ylabel(r'$C_L$')
        pyl.legend(legend,loc = 'best')
        pyl.gca().yaxis.set_label_coords(-0.05,0.5)
        pyl.grid(gridarg)
                        
#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACLiftSurf, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'S':
            return self._CalcArea()
        elif key == 'AR':
            return self._CalcAR()
        elif key == 'b':
            return self._CalcSpan()
        elif key == 'Weight':
            return self._CalcWeight()
        elif key == 'o_eff':
            return self._CalcOEff()
        #
        # If no calculation exist return None
        #
        return None


#===============================================================================
    def Refresh(self):
        """
        Updates all internal calculation arrays.

        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACLiftSurf, self).Refresh()
        self.param.refreshing = True

        #
        # This ensures that the variable gets updated
        #
        self.param.MAC = None

        #
        # Update all internal arrays and variables
        #
        self._CheckConsistent()
        self._CalcAirfoils()
        self._CalcChord()
        self._CalcChordSpanLoc()
        self._CalcDihedral()
        self._CalcSweep()
        self._MAC()

        #
        # Now update the drawing of the lifting surface
        #
        self.UpdateDrawing()
        
        #self.WingWeight.Refresh()

        self.param.refreshing = False

#===============================================================================
    def UpdateDrawing(self):
        """
        Updates all drawings of the lifting surface
        """
        #
        # The wing is special and may need to be refereshed at this point
        #
        if not self.param.refreshing and self.dirty:
            self.Refresh()
            
        #
        # Update the control surfaces 
        #
        for Control in self.Controls.__dict__.itervalues():
            Control.Refresh()

        #
        # This will indirectly recursively refresh all the winglets
        #
        for Winglet in self.Winglets.__dict__.itervalues():
            Winglet.Refresh()

        #
        # Update the wing weight calculation
        #
        if self.WingWeight is not None and not self.WingWeight.param.refreshing:
            self.WingWeight.Refresh()

        LiftSurf    = self.Drawings.LiftSurf
        LiftSurfMir = self.Drawings.LiftSurfMir
        AirfoilR    = self.Drawings.AirfoilR
        AirfoilT    = self.Drawings.AirfoilT

        X = self.GetX() / IN
        span = self.param.span
        b = self.b
        b = b/2. if self.FullWing else b
        
        #
        # Get root and tip airfoil properties
        #
        Cr   = self.Chord(0 * IN) / IN #Root chord
        LEr  = self.LE(0 * IN) / IN    #Root leading edge
        Midr = self.Mid(0 * IN) / IN   #Root lower surface edge

        Ct   = self.Chord(b) / IN      #Tip chord
        LEt  = self.LE(b) / IN         #Tip leading edge
        Midt = self.Mid(b) / IN        #Tip lower surface edge
        
        #
        # Get the airfoil coordiates at the incidence angle
        #
        xaf, yaf = self.af.XYAtAngle(self.i)
       
        #
        # Fill the plooting arrays
        # This could be slicker by doing some rotations and projections...
        #
        npt = 2*len(span)

        if abs(self.Axis[0]) == 1: #Horizontal surface
            # Top view
            LiftSurf.xt = npy.empty(npt)
            for i in range(len(span)):
                LiftSurf.xt[i]      = self.LE(span[i] * IN) / IN
                LiftSurf.xt[-(i+1)] = self.TE(span[i] * IN) / IN

            LiftSurf.yt = npy.ones(npt) * X[1]
            for i in range(len(span)):
                LiftSurf.yt[i]      += span[i]*self.Axis[0]
                LiftSurf.yt[-(i+1)] += span[i]*self.Axis[0]

            # Front view
            LiftSurf.yf = LiftSurf.yt

            LiftSurf.zf = npy.empty(npt)
            for i in range(len(span)):
                LiftSurf.zf[i]      = self.Upper(span[i] * IN) / IN
                LiftSurf.zf[-(i+1)] = self.Lower(span[i] * IN) / IN

            #
            # Generate a mirror of the half wing to draw a full wing
            #
            if self.FullWing:
                LiftSurfMir.xt =  LiftSurf.xt
                LiftSurfMir.yt = -LiftSurf.yt + 2*X[1]

                LiftSurfMir.yf = -LiftSurf.yf + 2*X[1]
                LiftSurfMir.zf =  LiftSurf.zf

            #
            # Draw the airfoil for a Side view
            #
            AirfoilR.xs = xaf * Cr + LEr
            AirfoilR.zs = yaf * Cr + Midr

            AirfoilT.xs = xaf * Ct + LEt
            AirfoilT.zs = yaf * Ct + Midt

        else: #Vertical surface

            # Side view
            LiftSurf.xs = npy.empty(npt)
            for i in range(len(span)):
                LiftSurf.xs[i]      = self.LE(span[i] * IN) / IN
                LiftSurf.xs[-(i+1)] = self.TE(span[i] * IN) / IN

            LiftSurf.zs = npy.ones(npt) * X[2]
            for i in range(len(span)):
                LiftSurf.zs[i]      += span[i]*self.Axis[1]
                LiftSurf.zs[-(i+1)] += span[i]*self.Axis[1]

            # Front view
            LiftSurf.zf = LiftSurf.zs

            LiftSurf.yf = npy.empty(npt)
            for i in range(len(span)):
                LiftSurf.yf[i]      = self.Upper(span[i] * IN) / IN
                LiftSurf.yf[-(i+1)] = self.Lower(span[i] * IN) / IN
            
            #
            # Generate a mirror of the half wing to draw a full wing
            #
            if self.FullWing:
                LiftSurfMir.xs =  LiftSurf.xs
                LiftSurfMir.zs = -LiftSurf.zs + 2*X[2]

                LiftSurfMir.yf =  LiftSurf.yf
                LiftSurfMir.zf = -LiftSurf.zf + 2*X[2]

            #
            # Draw the airfoil for a Top view
            #
            AirfoilR.xt = xaf * Cr + LEr 
            AirfoilR.yt = yaf * Cr + Midr

            AirfoilT.xt = xaf * Ct + LEt
            AirfoilT.yt = yaf * Ct + Midt
 
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
        super(ACLiftSurf,self).Draw(fig, top, side, front)
          
        for Control in self.Controls.__dict__.itervalues():
            Control.Draw(fig, top, side, front)
             
        self._DrawWinglets(self.Winglets, fig, top, side, front)

#===============================================================================
    def _DrawWinglets(self, Winglets, fig, top, side, front):
        """
        Draws winglets of a lifting surface. This is needed due to to recusive nature of winglets

        Inputs:
            Winglets - A class of winglets
            fig      - Integer number of figure to draw into
            top      - Subplot for top view
            side     - Subplot for side view
            front    - Subplot for front view
        """        
        for Winglet in Winglets.__dict__.itervalues():
            Winglet.Draw(fig, top, side, front)
            self._DrawWinglets(Winglet.Winglets, fig, top, side, front)
       

#===============================================================================
    def _CheckConsistent(self):
        """
        For example, only two of S, b, and AR should be specified, not all three.
        """
        NoneList = self.NoneList
        Axis = self.Axis
        
        if not ((abs(Axis[0]) == 1 and Axis[1] == 0) or \
                (abs(Axis[1]) == 1 and Axis[0] == 0)):
            raise ACLiftSurfError("Axis must be set to (+-1,0) for a Horizontal or (0,+-1) for Vertical lifting surface in ", self.name)
        

#===============================================================================
#
# Internal functions to the class for Aerodynamic calculations
#
#===============================================================================
    def Re(self):
        """
        Virtual function to calculate Reynolds number. This function must be implimented
        in a derived class.
        """
#        return 5e5 #A hack for debugging
        raise ACLiftSurfError("Re not implimented in ", self.__class__.__name__)
        pass

#===============================================================================
    def _CalcOEff(self):
        """
        Calculates the oswald efficiency using lifint line theory
        """
        if not self.FullWing:
            return 0.97
        else:
            return self.LLT.OswaldEff()

#===============================================================================
    def _CalcArea(self):
        """
        Calculates the area based on span and aspect ration.
        
        This is more of a backup function. Typically area is calculated based on 
        a desired lift of a volume coefficient
        """
        return self.b**2/self.AR

#===============================================================================
    def _CalcAR(self):
        "Calculates aspect ratio."
        return self.b**2 / self.S

#===============================================================================
    def _CalcSpan(self):
        "Calculates wing span."
        b2 = self.AR * self.S / (IN**2)
        return math.sqrt(b2).real * IN

#===============================================================================
    def _CalcWeight(self):
        """
        Calculates the weight of the wing if the Wing Weight class exists
        """
        Weight = 0*LBF
        if self.WingWeight is not None:
            Weight = self.WingWeight.Weight
            
        Weight += self._WingletWeight()
        Weight += self._ControlWeight()
        
        return Weight

#===============================================================================
    def _WingletWeight(self):
        """
        Adds the weight of the winglets
        """
        Weight = 0*LBF
        for Winglet in self.Winglets.__dict__.itervalues():
            Weight += Winglet.Weight
        
        return Weight

#===============================================================================
    def _ControlWeight(self):
        """
        Adds the weight of the control surfaces
        """
        Weight = 0*LBF
        for Control in self.Controls.__dict__.itervalues():
            Weight += Control.Weight
            ServoWeight = Control.Servo.Weight
            if self.FullWing and Control.RootFb > 0.0:
                ServoWeight *= 2
            if self.Symmetric:
                ServoWeight *= 2
            Weight += ServoWeight
        
        return Weight

#===============================================================================
    def PrintCompufoilGenRibs(self):
        """
        Dumps information for generating Compufoil GenRibs
        """
        
        y_old = 0*IN
        
        def WriteSection(y):
            print
            print 'Section    : ', AsUnit( y - y_old, 'in' )
            print 'Chord Root : ', AsUnit( self.Chord(y_old), 'in' )
            print 'Chord Tip  : ', AsUnit( self.Chord(y), 'in' )
            print "Sweep      : ", AsUnit( self.LE(y) - self.LE(y_old), 'in' )
            print '25% C Root : ', AsUnit( self.Chord(y_old) * 0.25, 'in' )
            print '25% C Tip  : ', AsUnit( self.Chord(y) * 0.25, 'in' )
            
        b = self.b
        b = b/2. if self.FullWing else b

        print
        print "Airfoil : ", self.Airfoil
        print "Span    : ", AsUnit(self.b, 'in')
        
        #
        # Write out all the sections of the wing
        #
        for i in range(1, len(self.param.span)):
            span = self.param.span[i]*IN

            #
            # Defining spanwise locations of the wing are simple
            #
            WriteSection(span)
            y_old = span
            

            #
            # Control surfaces require a little more work...
            #
            for control in self.Controls.__dict__.itervalues():
                rootFb = control.RootFb()
                tipFb  = control.TipFb()
                
                rootChord = control.RootChord()
                tipChord  = control.TipChord()
                
                print
                print "Control Surface : ", control.name
                print "Root Loc and Chord : ", AsUnit( rootChord*(1 - control.Fc), 'in' ), AsUnit( rootChord*control.Fc, 'in' )
                print "Tip  Loc and Chord : ", AsUnit( tipChord*(1 - control.Fc), 'in' ), AsUnit( tipChord*control.Fc, 'in' )
                                    
                # Check to see if the root is in between this and the next span fraction
                if self.Fb[i-1] < rootFb and self.Fb[i] > rootFb:
                    WriteSection(rootFb*b)
                    y_old = rootFb*b
                                                    
                # Check to see if the tip is in between this and the next span fraction
                if self.Fb[i-1] < tipFb and self.Fb[i] > tipFb:
                    WriteSection(tipFb*b)
                    y_old = tipFb*b
                


#===============================================================================
    def WriteAVLWing(self, filename):
        """
        Writes out an AVL input file for the lifting surface

        Inputs:
            filname - The file name of the AVL file to be written
        """
        fAVL = open(filename, 'w')
        
        mach = 0. #Mach number is just hard coded to zero for now...

        iYsym = 0  # This is an integer
        iZsym = 0  # This is an integer
        Zsym  = 0. # This is a floating point

        Sref = self.S / (IN**2)
        Cref = self.MAC() / IN
        Bref = self.b / IN

        Xref = self.X[0] / IN
        Yref = 0.0 if self.Symmetric else self.X[1] / IN
        Zref = self.X[2] / IN

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
    def _WriteAVLWinglets(self, Winglets, fAVL):
        """
        Writes winglets of a lifting surface. This is needed due to to recusive nature of winglets

        Inputs:
            Winglets - A class of winglets
            fAVL     - The open file to write to
        """        
        for Winglet in Winglets.__dict__.itervalues():
            Winglet._WriteAVLSurface(fAVL)
            self._WriteAVLWinglets(Winglet.Winglets, fAVL)


#===============================================================================
    def _WriteAVLSurface(self, fAVL, Sspace = None):
        """
        Writes out the lifting surface to an AVL file

        Inputs:
            fAVL - The open file to write to
            Sspace - Spanwise distribution if latice points (-2.0 to cluster at tip, -1.0 to cluster at both ends)
        """
        if self.dirty:
            self.Refresh()

        i = self.i / ARCDEG
        spanloc = self.param.span
        
        FullWing = self.FullWing and not (self.X[1] / IN == 0 and abs(self.Axis[0]) == 1)

        fAVL.write('#==============================================================\n')
        fAVL.write('#\n')
        fAVL.write('SURFACE\n')
        fAVL.write(self.name + '\n')
        
        #
        # If it is a full wing, cluster the vortecies at both ends by setting Sspace = -1, otherwise Sspan = -2
        #
        if Sspace is None:
            if FullWing:
                Sspace = -1.0
            else:
                Sspace = -2.0
            
        fAVL.write('10  1.0  20  '+str(Sspace)+'  |  Nchord   Cspace   Nspan  Sspace\n')
        #
        # Determine of symmetry is appropriate
        #
        if self.Symmetric or (self.FullWing and self.X[1] / IN == 0 and abs(self.Axis[0]) == 1):
            fAVL.write('\nYDUPLICATE\n' + str(0.0) + '\n')
                
        fAVL.write('\nANGLE\n' + str(i) + '\n')
        fAVL.write('\nINDEX\n' + str(self.param.index) + '\n')


        #
        # If this is a full wing, frist write the negative half
        #
        if FullWing:
            self._WriteAVLSpanSections(-spanloc[:0:-1], fAVL)
 
        #
        # Write out the positive spanwise locations
        #
        self._WriteAVLSpanSections(spanloc, fAVL)
        
        #
        # Write the winglets for this lifting surface
        #
        self._WriteAVLWinglets(self.Winglets, fAVL)


#===============================================================================
    def _WriteAVLSpanSections(self, spanloc, fAVL):
        """
        Writes out spanwise locations to an AVL file

        Inputs:
            spanloc - The spanwise locations to write
            fAVL    - The open file to write to
        """
        
        b = self.b
        b = b/2. if self.FullWing else b
        #
        # Write out all the sections of the wing
        #
        for i in range(len(spanloc)):
            span = spanloc[i]

            #
            # Defining spanwise locations of the wing are simple
            #
            self._WriteAVLSection(fAVL, span*IN)

            #
            # Control surfaces require a little more work...
            #
            for control in self.Controls.__dict__.itervalues():
                rootFb = control.RootFb()
                tipFb  = control.TipFb()

                #
                # Chech if the control surface lines up with the root
                #
                if i == 0:  
                    if rootFb < 0.00001: #The section has already been written, just write the control surface
                        self._WriteAVLControl(fAVL, control)
                    elif rootFb < self.Fb[0]-0.00001: #The control surface and section must be written
                        self._WriteAVLSection(fAVL, rootFb*b)
                        self._WriteAVLControl(fAVL, control)

                    if tipFb < self.Fb[0]-0.00001: #The control surface and section must be written
                        self._WriteAVLSection(fAVL, tipFb*b)
                        self._WriteAVLControl(fAVL, control)
                    continue

                #Check to see if the current section is between the end points of the control surface
                if rootFb < self.Fb[i-1] and self.Fb[i-1] < tipFb:
                    self._WriteAVLControl(fAVL, control)

                # Check to see if the span location is the tip or root of the control surface
                elif abs(self.Fb[i-1] - rootFb) < 0.001 or abs(self.Fb[i-1] - tipFb) < 0.001:
                    self._WriteAVLControl(fAVL, control)
                                    
                # Check to see if the root is in between this and the next span fraction
                if self.Fb[i-1] < rootFb and self.Fb[i] > rootFb:
                    self._WriteAVLSection(fAVL, rootFb*b)
                    self._WriteAVLControl(fAVL, control)
                                                    
                # Check to see if the tip is in between this and the next span fraction
                if self.Fb[i-1] < tipFb and self.Fb[i] > tipFb:
                    self._WriteAVLSection(fAVL, tipFb*b)
                    self._WriteAVLControl(fAVL, control)
                    continue

                #Finish of the section if it is at the tip
#                if i == len(spanloc) - 1 and abs(tipFb - 1.) < 0.001:
#                    self._WriteAVLControl(fAVL, control)
                        

#===============================================================================
    def _WriteAVLSection(self, fAVL, span):
        """
        Adds the properties to the AVL section

        Inputs:
            fAVL - The open file to write to
            span - scalar of the span location
        """

        Xle   = self.LE(span) / IN
        if abs(self.Axis[0]) == 1:
            Yle   = (self.X[1] + span * self.Axis[0]) / IN
            Zle   =  self.Mid(span) / IN
        else:
            Yle   =  self.Mid(span) / IN
            Zle   = (self.X[2] + span * self.Axis[1]) / IN

        Chord = self.Chord(span) / IN
        AFIL  = self.af.GetFilename()

        fAVL.write('#-----------------------\n')
        fAVL.write('SECTION\n')
        fAVL.write(str(Xle) + ' ' + str(Yle) + ' ' + str(Zle) + '   ' + str(Chord) +
                   '  0 0 0 | Xle Yle Zle   Chord   Ainc   [ Nspan Sspace ]\n')
        fAVL.write('AFIL\n')
        fAVL.write('"' + AFIL + '"\n')

#===============================================================================
    def _WriteAVLControl(self, fAVL, control):
        """
        Adds the control properties to a AVL section

        Inputs:
            fAVL    - The open file to write to
            control - The control surface
        """
        Gain    = control.Gain
        Fc      = control.Fc
        SgnDup  = control.SgnDup
        
        fAVL.write('CONTROL\n')
        fAVL.write(control.name + '  ' + str(Gain) + '  ' + str(1.-Fc) + '  0. 0. 0.  ' + str(SgnDup) +
                   ' | Name  gain  Xhinge  XYZhvec  SgnDup\n')
    
#===============================================================================
    def _CheckSpanLoc(self, y):
        """
        Checks if a spanwise location is valid
        """
        #
        # Get the wing span
        #
        b = self.b

        br = b       # 'Right' wing end
        bl = 0 * IN  #  'Left' wing end
        if self.FullWing:
            br /= 2.
            bl = -br

        yi = y / IN
        
        message = "Invalid wing span location '" + str(y) + \
                                  "' in range " + str(br) + " to " + str(bl) + \
                                  " for "
        
        if isinstance( yi, npy.ndarray ):
            if (npy.all( npy.round(yi, 4) > round(br/IN, 4) )  or 
                npy.all( npy.round(yi, 4) < round(bl/IN, 4) ) ):
                raise ACLiftSurfError(message, self.name)
        else:
            if ( round(yi, 4) > round(br/IN, 4) ) or ( round(yi, 4) < round(bl/IN, 4) ):
                raise ACLiftSurfError(message, self.name)
        #
        # All span positions are stored in a positive coordinate.
        #
        if isinstance( yi, npy.ndarray ):
            return -y if npy.all(yi < 0) else y
        
        return -y if yi < 0 else y
    
#===============================================================================
    def _CalcChordSpanLoc(self):
        """
        Calculates the spanwise location of each chord
        """
        b = self.b / IN
        Fb = self.Fb

        b = b/2. if self.FullWing else b

        #
        # The last chord location must always be 1
        #
        if Fb[-1] != 1:
            raise ACLiftSurfError('Fb[-1] must be set to 1.', self.name)

        self.param.span = npy.empty(len(Fb)+1)
        self.param.span[0] = 0.
        for i in range(len(Fb)):
            self.param.span[i+1] = Fb[i] * b


#===============================================================================
    def _CalcChord(self):
        """
        Fills the internal chord array

        The root chord is calculated based on the recursion

        S=fac\sum_{i}S_{i}
        S_{i}=Fb_{i}\frac{b}{2}\frac{\left(1+TR_{i}\right)}{2}C_{i}
        C_{i+1}=C_{i}TR_{i}
        """

        b   = self.b / IN
        TR  = self.TR
        Fb  = self.Fb
        Gam = self.Gam
        S   = self.S / (IN**2)


        #
        # If a constant upper surface is desired just set the dihedral to zero (should be close enough...)
        #
        if self.ConstUpper:
            Gam = [0. * ARCDEG for i in range(len(Fb))]

        #
        # Check that the array sizes are correct
        #
        TR, Fb  = self._CheckSameLen(TR, "Taper ratio", Fb, "span fraction")
        TR, Gam = self._CheckSameLen(TR, "Taper ratio", Gam, "dihedral angles")

        #
        # The last chord location must always be 1
        #
        if Fb[-1] != 1:
            raise ACLiftSurfError('Fb[-1] must be set to 1.', self.__name__)

        #
        # First calculate the root chord
        #
        # Compute the denomenator
        Si = b  * Fb[0] * (1 + TR[0])/2 * math.cos(Gam[0] / RAD).real
        for i in range(1, len(TR)):
            TR1 = 1
            for j in range(i):
                TR1 *= TR[j]

            Si += b * (Fb[i]-Fb[i-1]) * TR1 * (1 + TR[i])/2 / math.cos(Gam[i] / RAD).real

        #
        # Note the the chord locations are stored in inch units
        #
        C = Cr = S/Si   #Root chord
        self.param.chord = npy.empty(len(TR)+1)
        self.param.chord[0] = Cr

        # Calculate remaining chords
        for i in range(len(TR)):
            C *= TR[i]
            self.param.chord[i+1] = C


#===============================================================================
    def _CalcConstantUpper(self):
        """
        Calculates dihedral angles to get a constant a upper surface
        A constant lower surface is given by setting all the dihedral angles to zero
        """

        if not self.ConstUpper:
            return False

        Fb  = self.Fb

        af_ymax = npy.max(self.af.y) #Maximum height of the airfoil

        #
        # 'chord' must be calculated and be in inches
        #
        dihed = self.param.dihed = npy.zeros(len(Fb)+1)
        chord = self.param.chord
        dihed[0] = 0
        for i in range(len(Fb)):
            dihed[i+1] = dihed[i] + af_ymax*(chord[i]-chord[i+1])

        return True

#===============================================================================
    def _CalcDihedral(self):
        """
        Updates the dihedral offset
        """
        #
        # Calculate the dihedrals to get a constant upper surfase if desired
        #
        if(self._CalcConstantUpper()):
            return

        #
        # Calculate the dihedral based on the given angles
        #
        b   = self.b / IN
        Gam = self.Gam
        Fb  = self.Fb

        #
        # Check that the lists are the same length
        #
        Fb, Gam = self._CheckSameLen(Fb, "Span fraction", Gam, "dihedral angles")

        b = b/2. if self.FullWing else b

        dihed = self.param.dihed = npy.zeros(len(Gam)+1)

        dihed[0] = 0 # Never offset the root chord
        dihed[1] = b * Fb[0] * math.tan(Gam[0] / RAD).real
        for i in range(1,len(Gam)):
            dihed[i+1] = dihed[i] + b * (Fb[i]-Fb[i-1]) * math.tan(Gam[i] / RAD).real


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

        sweep = self.param.sweep = npy.zeros(len(Fb)+1)

        sweep[0] = 0 # Never offset the root chord

        if self.CEdge == 'LE':
            for i in range(len(Fb)):
                sweep[i+1] = sweep[i] + 0.25*(chord[i+1] - chord[i])

        elif self.CEdge == 'TE':
            for i in range(len(Fb)):
                sweep[i+1] = sweep[i] - 0.75*(chord[i+1] - chord[i])

        return True


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

        b     = self.b / IN
        Lam   = self.Lam
        Fb    = self.Fb
        chord = self.param.chord
        sfc   = self.SweepFc

        b = b/2. if self.FullWing else b

        #
        # Check that the lists are the same lenght
        #
        Fb, Lam = self._CheckSameLen(Fb, "Span fraction", Lam, "sweep angles")
        
        Lam = [Lam[i] / RAD for i in range(len(Lam))]

        sweep = self.param.sweep = npy.zeros(len(Lam)+1)
        #
        # First sweep the quarter chord back to get a constant line
        # along the desired chord fraction to sweep about
        #
        for i in range(len(Fb)):
            sweep[i+1] = (0.25 - sfc)*(chord[i+1] - chord[i])

        #
        # Now apply the desired sweep
        #
        sweep[0] = 0 # Never offset the root chord
        sweep[1] += b * Fb[0] * math.tan(Lam[0]).real
        for i in range(1,len(Lam)):
            sweep[i+1] += sweep[i] + b * (Fb[i]-Fb[i-1]) * math.tan(Lam[i]).real


#===============================================================================
    def _CalcAirfoils(self):
        """
        Updates the airfoil sections if they have changed
        """
        airfoil = self.Airfoil

        if self.af is None or self.af.name != airfoil:
            self.af = ACAirfoil(airfoil)
            
        self.af.Inverted = self.Inverted

################################################################################
class ACWinglet(ACLiftSurf):
    """
    A class for creating winglets
    
    Note:
        AR does not exist for a winglet because the root chord is determined from the main wing
    
    """
    def __init__(self, name, index, Wing):
        """
        Inputs:
            name  - The name of the winglet
            index - Index of the winglet for AVL
            Wing  - The wing this winglet belongs to
        """
        super(ACWinglet, self).__init__(index)
        
        #
        # Save of the name of the winglet
        #
        self.name = name
        
        self.Axis = (0, 1) #Default to a vertical winglet (A horizontal one does not make much sense...)
        
        #
        # Save of a reference to the main wing
        #
        self.param.Wing = Wing

#===============================================================================
    def Re(self):
        """
        Calculates the Reynolds number for the wing.

        This is a relatively cheap calculation and is therefore left
        as a calculation rather than storing a self.param.Re
        to reduce the interdependencies.
        """
        Alt = self.param.Wing.GetAlt_LO()
        V   = self.GetV_LO()
        rho = AeroUtil.DensAlt(Alt)
        mu  = AeroUtil.Mu(AeroUtil.TempAlt(Alt))
        mac = self.MAC()

        return rho * V * mac / mu

#===============================================================================
    def GetV_LO(self):
        """
        Retreives the lift of velocity for the winglet. This was taken from the 
        main wing to begin with. The method GetV_LO allows a winglet to 'pretend'
        that it is a 'main wing'.
        """
        return self.param.Wing.GetV_LO()

#===============================================================================
    def GetAlt_LO(self):
        """
        Retreives the lift of altitude for the winglet. This was taken from the 
        main wing to begin with. The method GetV_LO allows a winglet to 'pretend'
        that it is a 'main wing'.
        """
        return self.param.Wing.GetAlt_LO()
        
#===============================================================================
    def Refresh(self):
        """
        Updates all the internal quantities.
        """
        self.__dict__['dirty'] = False
        self.param.refreshing = True
        
        self.CopyVisibility(self.param.Wing)
        
        #
        # Due to the cyclic behavior of the area calculation, it must be updated before anything else.
        #
        self._UpdateArea()
        self._UpdateSpan()
        
        #
        # Set the winglet to be at the tip of the wing.
        # This is done before calling the baseclass where it will be drawn
        #
        self.X = self.param.Wing.Tip()
              
        #
        # Now update everything in the base class
        #
        super(ACWinglet, self).Refresh()
       
        self.param.refreshing = False    

#===============================================================================
    def _GetChord(self):
        """
        Calculates the array of chords. This is simple because the root chord is 
        given by tip chord of the main wing. The returned list is not a scalar but a 
        list of numbers in units of inches.
        """
        Wing = self.param.Wing
        b   = self._Wing_WingTip()
        
        TR  = self.TR
        
        if op.isNumberType(TR):
            TR = [TR]
            
        #
        # Note the the chord locations are stored in inch units
        #
        C = Cr = Wing.Chord(b) / IN   #Root chord
        chord = npy.empty(len(TR)+1)
        chord[0] = Cr

        # Calculate remaining chords
        for i in range(len(TR)):
            C *= TR[i]
            chord[i+1] = C

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
        Calculates the area of the winglet
        """
        
        if not self.NoneList.has_key('S'):
            return
        #
        # The area is calculated based on the tip chord of the MainWing and
        # the taper ratios. Because the area is needed for so many other calculations,
        # the chord is also calculated here so the area can be updated before the base class is Refreshed.
        #
        
        b   = self.b / IN
        TR  = self.TR
        Fb  = self.Fb
        
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
        # Save of the area of the winglet
        #
        self.param.S = S*IN**2

#===============================================================================
    def _CalcSpan(self):
        """
        Returns the calculated span
        """
        if self.dirty: self.Refersh()
        return self.param.b

#===============================================================================
    def _UpdateSpan(self):
        """
        Calculates the span of the winglet
        """
        
        if not self.NoneList.has_key('b'):
            return
       
        #
        # The span is calculated based on the Area calculation.
        #
        
        TR  = self.TR
        Fb  = self.Fb
        S   = self.S / IN**2
        
        #
        # Check that the array sizes are correct
        #
        TR, Fb  = self._CheckSameLen(TR, "Taper ratio", Fb, "span fraction")

        #
        # Get the chords
        #
        chord = self._GetChord()

        #
        # This formula is based on how the area is calculated (Look at _UpdateArea)
        #
        b = S
        for i in range(len(TR)):
            b /= (chord[i] + chord[i+1])/2.*Fb[i]
       
        #
        # Save of the span of the winglet
        #
        self.param.b = b*IN

#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.
        """
        super(ACWinglet, self)._CheckConsistent()
                
        self._CheckEquation(['S','TR','b'])
        self._CheckEquation(['AR'], Need = 0)
            
#===============================================================================
    def _Wing_WingTip(self):
        """
        Returns span location of the wing-tip of the main wing
        """
        b = self.param.Wing.b
        return b/2. if self.param.Wing.FullWing else b
    
################################################################################
if __name__ == '__main__':
    # Create a lifting surface and try changing some parameters
    wing = ACLiftSurf(1)

    ail = wing.AddControl('ail')

    #
    # Because ACLiftSurf cannot calculate it's area (by design)
    # it is the first quantity set
    #
    wing.S = 900 * IN**2

    wing.AR = 16
    wing.b  = 10 * FT
    wing.TR = (1, 0.7, 0.6)
    wing.Fb = (0.4, 0.8, 1)
    wing.Airfoil = 'S1223'
#    wing.Set('Lam',[0., 0., 0.])
#    wing.Set('Gam',[0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG])
    wing.CEdge = 'LE'
    wing.ConstUpper = True
    wing.LLFile = 'LLTPolar.txt'

#    wing.LLTPlots(2)
    
    print 'Moment Loading', wing.MomentLoading(0*IN, 40*FT/SEC)

    Re = 5e5

#    print 'Wing area',wing.S
#    print 'Wetted area',wing.SWet()
#    print 'MAC',wing.MAC()
#    print 'CL',wing.CL(5*ARCDEG, Re)
#    print 'CDi',wing.CDi(5*ARCDEG, Re)
#    print 'CD',wing.CD(5*ARCDEG, Re)
#    print 'Down Wash', wing.DownWash(5*ARCDEG, Re)
#    print 'CL max',wing.CLmax(Re)
#    print 'CL min',wing.CLmin(Re)
#
#    tail = ACLiftSurf('VerticalTail')
#    tail.X = (50 * IN, 0* IN, 10 *IN)
#    tail.Axis = (0 , 1)
#    tail.S = 100 * IN**2
#    tail.b = 10 * IN
#    tail.TR = (0.7)
#    tail.Airfoil = ('../Airfoils/NACA1410.dat','../Airfoils/NACA1410')

    wing.FullWing = True
      
    wing.Draw()
    tail.Draw()
    pyl.show()
