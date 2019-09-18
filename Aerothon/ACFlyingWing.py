"""
A class for flying wings
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import M, FT, IN, ARCDEG, RAD, LBF, SEC, LBM, MIN
from scalar.units import AsUnit
from ACBase import ACComponent, Angle, Unitless, Length, Force, Time, Velocity, g
from ACAircraft import ACAircraft
from ACTableBase import ACTableBase
from ACWing import ACMainWing
from ACTail import  ACTailSurf, ACHorizontalTailSurf
from ACLandingGear import ACLandingGear
from ACFuselage import ACFuselage
from ACWeightTable import ACWeightTable
import pylab as pyl
import numpy as npy
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import AeroUtil
import math



################################################################################
class ACFlyingWing(ACAircraft):
    """
    A class for working with a flying wing aircraft

    Attributes:
        CLHTSlopeAt      - Fuselage angle of attack used to evaluate slope of the horizontal tail CL
        RotationAngle    - The angle the aircraft can rotate before the tail hits the ground
        NoseGearOffset   - Offset of the nose gear from the front bulkhead (default 0*in)
        VTailPos         - Horizontal location of the vertical tail in percentage of the horizonalt tail semi-span
        EngineAlign      - Allignment of the engine (1 - top, 0 - middle, -1 - bottom)
        WingXOmega       - Relaxation factor for positioning the wing. Higher number relaxes the solution process. (default 1)
        WingXMaxIt       - Maximum number of iterations used to position the wing (default 10)
    """
#===============================================================================
    def __init__(self):
        super(ACFlyingWing, self).__init__()
        UnitList = self.UnitList

        #
        # Aircraft attributes
        #
        self.NoneList['DWSlopeAt']        = None        ; UnitList['DWSlopeAt']      = Angle
        self.__dict__['RotationAngle']    = 0 * ARCDEG  ; UnitList['RotationAngle']  = Angle
        self.__dict__['NoseGearOffset']   = 0*IN        ; UnitList['NoseGearOffset'] = Length
        self.__dict__['VTailPos']         = 0           ; UnitList['VTailPos']       = Unitless
        self.__dict__['EngineAlign']      = 0           ; UnitList['EngineAlign']    = Unitless
        self.__dict__['WingXOmega']       = 1           ; UnitList['WingXOmega']     = Unitless
        self.__dict__['WingXMaxIt']       = 10          ; UnitList['WingXMaxIt']     = Unitless
        
        self.__dict__['dCL_de'] = 1/ARCDEG              ; UnitList['dCL_de']     = Unitless
        self.__dict__['dCM_de'] = 1/ARCDEG              ; UnitList['dCM_de']     = Unitless
        
        #
        # Aircraft components
        #
        self.__dict__['Wing']       = None #ACMainWing(1)
        self.__dict__['VTail']      = ACTailSurf(self.Wing, 3)
        self.__dict__['NoseGear']   = ACLandingGear()
        self.__dict__['MainGear']   = ACLandingGear()
        self.__dict__['Fuselage']   = None #ACFuselage()
        self.__dict__['Propulsion'] = None
        #
        # Make sure the aircraft gets dirtied if any of the following are dirtied
        #
        self.VTail.Parents.append( self )
        self.NoseGear.Parents.append( self )
        self.MainGear.Parents.append( self )
        
        #
        # Give some approprate names
        #
        self.VTail.name = "Vertical Tail"
        self.NoseGear.name   = "Nose Gear"
        self.MainGear.name   = "Main Gear"
       
        #
        # Main Gear is always symmetric
        #
        self.MainGear.Symmetric = True 
        
        rud = self.VTail.AddControl('Rudder')

        #
        # Rudder typically extends the entire vertical tail
        #
        rud.Fb = 1.
        rud.Fc = 0.25
        rud.Ft = 0.
        
        #
        # Add a drawing for the rotation angle
        #
        self.Drawings.Rotation = ACComponent.ComponentDrawing()
        self.Drawings.Rotation.color = 'r'

        self.Drawings.Tipping = ACComponent.ComponentDrawing()
        self.Drawings.Tipping.color = 'k'


#===============================================================================
    def CM(self, alpha2dw, del_e = 0*ARCDEG, Xcg = None):
        """
        Calculates the total moment coefficient of the aircraft

        CM = CLw*(Xcg-Xac)/macw + CMw + eta*Sht/Sw*CLht*(Xcg-Xht)/macw + CMht*Sht/Sw

        Inputs:
            alhpah2dw - Angle of attack the Wing
            del_e     - Elevator deflection
            Xcg       - Center of gravity
        """
        if self.dirty:  self.Refresh()

        Wing  = self.Wing

        if Xcg is None:
            Xcg = self.Xcg()

        Sw   = Wing.S
        CLw  = Wing.CL(alpha2dw) + del_e*self.dCL_de
        CMw  = Wing.CM(alpha2dw) + del_e*self.dCM_de
        Xac  = Wing.Xac()
        macw = Wing.MAC()

        return CLw*(Xcg-Xac)/macw + CMw

#===============================================================================
    def dCM_da(self, del_e = 0*ARCDEG, Xcg = None):
        """
        Calculates the derivative of total moment coefficient w.r.t angle of attack of the aircraft

        Inputs:
            del_e     - Elevator deflection
            Xcg       - Center of gravity of the aircraft
        """
        if self.dirty:  self.Refresh()

        Wing  = self.Wing
        
        a2dw = self.param.CmSlopeAt
        a3dw = self.CMSlopeAt

        if Xcg is None:
            Xcg = self.Xcg()

        return ((self.CM(a2dw[1], del_e, Xcg = Xcg) - \
                 self.CM(a2dw[0], del_e, Xcg = Xcg))/(a3dw[1] - a3dw[0]))

#===============================================================================
    def Cmadot(self, Xcg = None):
        """
        TODO: Document what Cmadot is.
        
        It is for controls
        """
        if self.dirty:  self.Refresh()
       
        return 0
        #
        # Calculate the slope of the wing downwash
        #
        deda = self.dDW_da()
        
        #
        # Calculate the CL slope for the horizontal tail
        #
        dCL_daht = self.dCLht_da()

        
        HTail = self.HTail
        Wing  = self.Wing
        if Xcg is None:
            Xcg   = self.Xcg()

        eta  = HTail.eta
        Xht  = HTail.Xac()
        VCht = HTail.VC

        Xac  = Wing.Xac()
        macw = Wing.MAC()

        return -(2.*deda*eta*dCL_daht*VCht*(Xht - Xcg)/macw) / (1/RAD)

#===============================================================================
    def Cmq(self, Xcg = None):
        """
        TODO: Document what Cmq is.
        
        It is for controls
        """
        if self.dirty:  self.Refresh()
        
        return 0
    
        #
        # Calculate the CL slope for the horizontal tail
        #
        dCL_daht = self.dCLht_da()

        
        HTail = self.HTail
        Wing  = self.Wing
        if Xcg is None:
            Xcg   = self.Xcg()

        eta  = HTail.eta
        Xht  = HTail.Xac()
        VCht = HTail.VC

        Xac  = Wing.Xac()
        macw = Wing.MAC()

        return -(2.*eta*dCL_daht*VCht*(Xht - Xcg)/macw) / (1/RAD)
    
#===============================================================================
    def Czq(self):
        """
        TODO: Document what Czq is.
        
        It is for controls
        """
        if self.dirty:  self.Refresh()

        return 0
    
        #
        # Calculate the CL slope for the horizontal tail
        #
        dCLht_da = self.dCLht_da()

        HTail = self.HTail
        eta  = HTail.eta
        VCht = HTail.VC

        return -(2.*eta*dCLht_da*VCht) / (1/RAD)

#===============================================================================
    def dDW_da(self):
        """
        Calculates the wing downwash slope
        """
        return self._da(self.WingDownWash, self.param.dwSlopeAt, self.DWSlopeAt) 

#===============================================================================
    def dCLht_da(self):
        """
        Calculates the CL slope of the horizontal tail
        """
        return self._da(self.HTail.CL, self.param.ClhtSlopeAt, self.CLHTSlopeAt)

#===============================================================================
    def dCL_da(self):
        """
        Calculates the CL slope of the aircraft
        """
        return self._da(self.CL, self.param.ClSlopeAt, self.CLSlopeAt)
        
#===============================================================================
    def _PositionWingX(self):
        """
        Positions the wing such that the neutral point has the desired static margin.

        The neutral point is defined as the location where dCM_da == 0.
        """
        
        #
        # The incidence angle on the tail must be zero for this
        #
        Xnp = self.Xnp()
        
        self.Wing.X[0] += Xnp - self.Wing.Xac()
        
        self.Wing.SetDirty()

#===============================================================================
    def _Get_del_e_trim(self, Alpha2dw, Xcg = None):
        """
        Calculates the elevator deflection to trim the aircraft at a given angle of attack.

        Inputs:
            Alpha2dw - 2D angle of attack of the wing
            Xcg      - Center of gravity for the calculation
        """
        
        if Xcg is None:
            Xcg = self.Xcg()
            
        #
        # Setup the equation
        #
        def y(del_e):
            return self.CM(a, Xcg = Xcg, del_e = del_e*ARCDEG)

        #
        # Give the initial guesses and solve the equation
        #
        del_e_trim = []
        for a in Alpha2dw:
            try:
                
                if len(del_e_trim) == 0:
                    x0 = 0.0
                else:
                    x0 = del_e_trim[-1]
                    
                x1 = x0 + 1.0
                
                #del_e = SecantSolver(x0 = x0, x1 = x1, y = y, tol = 1.e-3, itmax = 500)
                del_e = brentq(f = y, a = -20, b = 20, rtol = 1.e-4, maxiter = 500)
                
                del_e_trim.append(del_e)

            except:
                del_e_trim.append(0)
                print "WARNING: Failed to compute trimmed elevator at 2D angle of attack ", AsUnit( a, "deg", '%3.3f' )
#            raise ACAircraftError("\nFailed to computed a trimmed elevator deflection." +
#                                  "\nThis typically means that the horizontal tail and/or elevator is inadequate.")
        #
        # Return the trimmed elevator deflection
        #
        return del_e_trim[0]*ARCDEG if len(del_e_trim) == 1 else npy.array(del_e_trim)*ARCDEG

#===============================================================================
    def _Update_del_e_trim(self):
        """
        Calculates a series of elevator deflections and saves them to an
        interpolating array to improve efficiency
        """
        AlphaRange = self.Wing.AlphaRange()

        del_e_trim = self._Get_del_e_trim(AlphaRange) / (ARCDEG)

        self.param.del_e_trim = interp1d(AlphaRange / (ARCDEG), del_e_trim, kind='linear')

#===============================================================================
    def del_e_trim(self, a2dw, Xcg = None):
        """
        Returns the trimmed elevator angle for a given wing angle of attack

        Inputs:
            a2dw - 2D angle of attack of the main wing
            Xcg  - A center of gravity
        """
        if self.dirty: self.Refresh()

        #
        # An interpolation is saved off for the CG of the aircraft.
        # If a different CG is required, the point will be calculated
        #
        if Xcg is None or Xcg == self.Xcg():
            return self.param.del_e_trim(a2dw / (ARCDEG)) * ARCDEG
        else:
            return self._Get_del_e_trim(a2dw, Xcg)

#===============================================================================
    def _PositionParts(self):
        """
        Positions all the parts on the aircraft
        """
        self._PositionLandingGear()
        self._PositionWing()
        self._PositionWingX()
        #
        # Wing needs to be positioned before the vertical tail
        #
        self._PositionTail()
        self._PositionPropulsion()

#===============================================================================
    def _PositionWing(self):
        """
        Positions the wing on the fuselage
        """
        Fuselage = self.Fuselage
        
        #
        # Position the wing on the fuselage
        #
        self.Wing.X[1] = 0*IN
        self.Wing.X[2] = Fuselage.Bottom(self.Xcg()) + Fuselage.Height(self.Xcg()) * self.WingFuseFrac

#===============================================================================
    def _PositionFuselage(self):
        """
        Positions the fuselage
        """
        pass

#===============================================================================
    def _PositionTail(self):
        """
        Positions the fuselage on the horizontal tail
        """
        Fuselage = self.Fuselage
        MainGear = self.MainGear
        Wing     = self.Wing
        VTail    = self.VTail
        A_Rot    = self.RotationAngle / RAD
        
        # Update CG
        self.VTail.SetAircraftXcg(self.Xcg())
                        
        TE = VTail.MaxTE()
        
        if Fuselage.Tail.Align is None:
            if VTail.FullWing:
                TE = VTail.TE(VTail.b/2)
            else:
                TE = VTail.MaxTE()

            Z_vtail = (TE - MainGear.X[0])*math.tan(A_Rot)
    
            if VTail.FullWing:
                Z_vtail += VTail.b/2
                
            Z_vtail = max(Z_vtail, Fuselage.Sections[-1].FrontBulk.Top())
            Z_vtail = max(Z_vtail, Fuselage.Sections[-1].BackBulk.Top() - Fuselage.Sections[-1].BackBulk.Height/2)
        else:
            Z_vtail = Fuselage.Sections[-1].BackBulk.Top() - Fuselage.Sections[-1].BackBulk.Height/2
                                

        VTail.X[1] = Wing.Tip()[1]*self.VTailPos
        VTail.X[2] = Z_vtail
        
        itmax = 50
        it = 0
        while( abs( Wing.TE(VTail.X[1]) - VTail.MaxTE() ) > 0.01*IN and it < itmax):
            VTail.L = Wing.TE(VTail.X[1]) - self.Xcg() - VTail.Chord(0*IN)*3./4
            it += 1
        
        Fuselage.TailBulk.X[0] = Wing.TE(0*IN)
        Fuselage.TailBulk.X[2] = Wing.X[2] - Fuselage.Sections[-1].BackBulk.Height/2
        

#===============================================================================
    def _PositionLandingGear(self):
        """
        Positions the nose and main landing gear
        """
        Fuselage   = self.Fuselage
        NoseGear   = self.NoseGear
        MainGear   = self.MainGear
        TAngle     = self.TippingAngle / RAD
        A_GR       = self.Alpha_Groundroll / RAD
        
        Fuselage.XOffset = self.Propulsion.Length()

        CG         = self.DesignCG()
        HCG        = CG[2] - Fuselage.Bottom(CG[0]) + MainGear.GearHeight
        X_MainGear = CG[0] + HCG*math.tan(TAngle)
                
        #
        # The fuselage is generated relative to the foremost bulkhead
        #
        NoseBulk    = Fuselage.GetNoseBulk()
        NoseSection = Fuselage.GetNoseSection()
        
        X_NoseGear = NoseBulk.X[0] + self.NoseGearOffset
        
        dZ_MG = NoseSection.Bottom(X_NoseGear) - Fuselage.Bottom(X_MainGear)
        dZ_NG = (X_MainGear - X_NoseGear)*math.tan(A_GR)
        
        NoseGear.StrutL = None
        NoseGear.GearHeight = MainGear.GearHeight + dZ_NG + dZ_MG
        
        #
        # Position the height of the nose bulkhead of the fuselage based on 
        # which gear makes the fuselage the highest
        #
        NoseBulk.X[2] = max(NoseGear.GearHeight, dZ_MG + MainGear.GearHeight) + NoseBulk.Height/2. + ( NoseBulk.Bottom() - NoseSection.Bottom(X_NoseGear) )
    
        Fuselage.Refresh()
        #
        # Finally position the landing gear
        #
        NoseGear.X[0] = X_NoseGear
        NoseGear.X[2] = NoseSection.Bottom(X_NoseGear)
    
        MainGear.X[0] = X_MainGear
        MainGear.X[2] = Fuselage.Bottom(X_MainGear)
        
        NoseGear.SetDirty()
        MainGear.SetDirty()

#===============================================================================
    def OverturnAngle(self):
        """
        Analyze Overturn Angle and throw warning if criteria not met
        Reymer p265 Fig 11.5, p266 paragraph 2
        """ 
        if self.dirty: self.Refresh()
        
        Fuselage   = self.Fuselage
        NoseGear   = self.NoseGear
        MainGear   = self.MainGear
        
        CG         = self.TotalCG()
        HCG        = (CG[2] - Fuselage.Bottom(CG[0]) + MainGear.GearHeight) / IN
        
        MainGearX = MainGear.WheelBase()
        NoseGearX = NoseGear.WheelBase()
        
        #Find Area of triangle formed by maingear, nosegear, and CG on ground plane
        x1 = MainGearX[0] / IN
        y1 = MainGearX[1] / IN
        x2 = NoseGearX[0] / IN
        y2 = NoseGearX[1] / IN
        x3 = CG[0] / IN
        y3 = CG[1] / IN
        LG_Area = 0.5 * ((x1 - x3)*(y2 - y1) - (x1 - x2)*(y3 - y1))
        LG_Area = math.fabs(LG_Area)
        WhlDist = ((x2-x1)**2 + (y2-y1)**2)**0.5
        WlCGLen = LG_Area / (0.5 * WhlDist)
        OTAngle = math.atan(HCG / WlCGLen) * RAD
        
        if OTAngle / (ARCDEG) > 63:
            print "    *** Overturn Angle is Greater than the suggested value of 63 degrees ***"
        
        return OTAngle

#===============================================================================
    def _PositionPropulsion(self):
        """
        Positions the propulsion system on the nose bulkhead
        """
        Fuselage   = self.Fuselage
        
        NoseBulk = Fuselage.GetNoseBulk()
        
        dirty = self.Propulsion.dirty
        
        EngineAlign = self.EngineAlign
        
        self.Propulsion.Engine.X[0] = Fuselage.GetNoseBulk().X[0]
        self.Propulsion.Engine.X[2] = Fuselage.GetNoseBulk().Top()*(1. + EngineAlign)/2.  + Fuselage.GetNoseBulk().Bottom()*(1. - EngineAlign)/2.
        
        self.Propulsion.Engine.Axis = [1, 0., 0.]
        self.Propulsion.Engine.Xat  = [1, 0.5, 0.5]

        #
        # Preserve the 'dirty' state of the propulsion
        #
        self.Propulsion.dirty = dirty
        self.Propulsion.Engine.dirty = dirty
        
        #
        # Also update the Vmax and altitude for the propulsion system
        #
        if self.Propulsion.Vmax != self.VmaxPlt:
            self.Propulsion.Vmax = self.VmaxPlt

        if self.Propulsion.Alt != self.GetAlt_LO():
            self.Propulsion.Alt = self.GetAlt_LO()
                  
#===============================================================================
    def Xnp(self):
        """
        Retrieves the neutral point of the aircraft.
        """
        if self.dirty:  self.Refresh()
        
        mac          = self.Wing.MAC()
        StaticMargin = self.StaticMargin
        Xnp          = self.Xcg() + mac*StaticMargin

        return Xnp

#===============================================================================
    def Xcg(self):
        """
        Retrieves the center of gravity of the aircraft.
        """
        if self.dirty:  self.Refresh()

        return self.Fuselage.AircraftCG()[0]

#===============================================================================
    def DesignCG(self):
        """
        Retrieves the desired design center of gravity of the aircraft.
        """
        if self.dirty:  self.Refresh()

        return self.Fuselage.AircraftCG()
    
#===============================================================================
    def TotalCG(self):
        """
        Calculates the empty CG of the aircraft
        """
        if self.dirty:  self.Refresh()
         
        EmptyWeight = self.EmptyWeight
        EmptyCG = self.EmptyCG()
        
        Payload = self.Fuselage.Payload
        
        TotalWeight = self.TotalWeight
        TotalCG = (Payload.Weight*Payload.CG() + EmptyWeight*EmptyCG)/TotalWeight

        return TotalCG

#===============================================================================
    def PayloadSize(self):
        """
        Calculates the payload height for the aircraft
        """
        if self.dirty: self.Refresh()
        
        return self.Fuselage.Payload.LWH
    
#===============================================================================
    def SetWing(self, Wing):
        """
        Assigns a new wing to the aircraft
        
        Inputs:
            Wing - The new wing to assign to the aircraft
        """
        #
        # Set the wing
        #
        self.Wing = Wing
        #
        # Update the references in the horizontal and vertical tail surfaces
        #
        self.VTail.SetWing( Wing )
        
        #
        # Make sure the aircraft is dirtied if the wing is dirtied
        #
        self.Wing.Parents.append( self )
         
        #
        # Give a better guess for the initial X value of the wing
        #
        if self.Fuselage is not None:
            self.Wing.X[0] = self.Fuselage.AircraftCG()[0]
            self.Wing.SetDirty()
               
#===============================================================================
    def SetFuselage(self, Fuselage):
        """
        Assigns a new fuselage to the aircraft
        
        Inputs:
            Fuselage - The new fuselage to assign to the aircraft
        """
        self.Fuselage = Fuselage
        #
        # Make sure the aircraft is dirtied if the fuselage is dirtied
        #
        self.Fuselage.Parents.append( self )

        #
        # Give a better guess for the initial X value of the wing
        #
        if self.Wing is not None:
            self.Wing.X[0] = self.Fuselage.AircraftCG()[0]
            self.Wing.SetDirty()

#===============================================================================
    def SetPropulsion(self, Propulsion):
        """
        Assigns a new propulsion system to the aircraft
        
        Inputs:
            Propulsion - The new propulsion to assign to the aircraft
        """
        self.Propulsion = Propulsion
        #
        # Make sure the aircraft is dirtied if the propulsion is dirtied
        #
        self.Propulsion.Parents.append( self )

#===============================================================================
    def WriteAVLMainWing(self, filename):
        """
        Writes out an AVL input file for just the main wing

        Inputs:
            filename - The name of the AVL file to be written
        """
        fAVL = open(filename, 'w')

        self._WriteAVLHeader(fAVL)
        self.Wing._WriteAVLSurface(fAVL)

        fAVL.close()

#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that the aircraft is not over or under specified
        """
        super(ACFlyingWing, self)._CheckConsistent()

        self._CheckEquation(['CLHTSlopeAt'], Need = 1)
        self._CheckEquation(['DWSlopeAt'], Need = 1)

        if self.DWSlopeAt[1] <= self.DWSlopeAt[0]:
            raise ACAircraftError(self.name + ": DWSlopeAt[1] must be grater than DWSlopeAt[0]")
               
#===============================================================================
    def Refresh(self):
        """
        Updates all the internal quantities.
        """
        super(ACFlyingWing, self).Refresh()
        self.param.refreshing = True
#        self.HTail.param.refreshing = True
#        self.VTail.param.refreshing = True
#        self.Wing.param.refreshing = True

        self._CheckConsistent()
        
        # 
        # Compute the 2D angle of attack given the fuselage angle of attack
        #
        self.param.A_GR        = self.Wing.Alpha2D(self.Alpha_Groundroll)
        self.param.a_Zero_CM   = self.Wing.Alpha2D(self.Alpha_Zero_CM)
        self.param.dwSlopeAt   = self.Wing.Alpha2D(self.DWSlopeAt)
        
        #
        # Position the part of the aircraft
        #
        self._PositionParts()
        
        
        #
        # Position the compute the incidence angle of the tail
        #
        #self._Updateiht() #This is now done as part of positioning the wing
        
        #
        # Need trimmed conditins to get lift of velocity
        #
        self._Update_del_e_trim()

        #
        # The lift of velocity is needed by many calculations
        #
        self._UpdateV_LO() 
        
        #
        # Update now the trimmed angle of attack can be found
        #
        self._UpdateAlphaTrim()
        self._UpdateVmax()
        self._UpdateV_max_climb()

        #
        # Now that everything is positioned and Refreshed, update the moments of inertia
        # 
        self._UpdateEmptyWeight()
        self._UpdateMOI()
        
        #
        # Now update the payload
        #
        self.Fuselage.Payload.Weight = self.PayloadWeight()
        self.Fuselage.Payload.Refresh()
                
        #
        # Update the drawings for the aircraft
        #
        self.UpdateDrawing()

        self.param.refreshing = False
#        self.HTail.param.refreshing = False
#        self.VTail.param.refreshing = False
#        self.Wing.param.refreshing = False

#===============================================================================
    def UpdateDrawing(self):
        """
        Updates all drawings of the tail aircraft
        """
        
        #
        # This is only applicable to the side view
        #
        Rotation  = self.Drawings.Rotation
        Tipping  = self.Drawings.Tipping
        WheelBase = self.MainGear.WheelBase()
        
        dx = (self.VTail.MaxTE() - WheelBase[0]) / (IN)
        
        #
        # Draw the rotation angle
        #
        RotAngle  = self.RotationAngle / (RAD)

        Rotation.xs = npy.empty(2)
        Rotation.zs = npy.empty(2)
                
        Rotation.xs[0] = WheelBase[0] / (IN)
        Rotation.zs[0] = 0
        
        Rotation.xs[1] = self.VTail.MaxTE() / (IN)
        Rotation.zs[1] = math.tan(RotAngle)*dx

        #
        # Draw the tipping angle
        #
        TipAngle  = self.TippingAngle / (RAD)

        Tipping.xs = npy.empty(2)
        Tipping.zs = npy.empty(2)

        Tipping.xs[0] = WheelBase[0] / (IN)
        Tipping.zs[0] = 0
        
        Tipping.xs[1] = self.VTail.MaxTE() / (IN)
        Tipping.zs[1] = math.tan(TipAngle)*dx
        
        #
        # Make sure the drawings are up to date
        #
        self.Wing.UpdateDrawing()
        self.VTail.UpdateDrawing()
        self.Fuselage.UpdateDrawing()
        self.NoseGear.UpdateDrawing()
        self.MainGear.UpdateDrawing()
        self.Propulsion.UpdateDrawing()
        

#===============================================================================
    def _CLComponents(self, a2dw, del_e):
        """
        Calculates all the components to build up the CL of the aircraft
        
        Inputs:
            a2dw  - 2D angle of attack of the main wing
            del_e - Elevator deflection
            V     - Velocity 
        """
        Wing = self.Wing

        #
        # Create a list of all lift components
        #
        CL = []
        CLLegend = []

        Wing = self.Wing

        Sw = Wing.S
        CLw = Wing.CL(a2dw)
        CL.append(CLw)
        CLLegend.append('Wing')
               
        return npy.array(CL), CLLegend
        
#===============================================================================
    def _CL(self, a2dw, del_e):
        """
        Calculates the total area weighted CL 
        of the aircraft with the given elevator deflection

        Inputs:
            a2dw - 2D angle of attack of the main wing
            del_e - Elevator deflection
        """
        CL, CLLegen = self._CLComponents(a2dw, del_e)
        
        CL = sum(CL)
        
        try:
            return CL[0] if len(CL) == 1 else CL
        except:
            return CL

#===============================================================================
    def _CDComponents(self, a2dw, del_e, V = None):
        """
        Calculates the total area weighted CD 
        of the aircraft with the given elevator deflection

        TODO: Check if the there should be sin(DHW) here or not
        
        Inputs:
            a2dw  - 2D angle of attack of the main wing
            del_e - Elevator deflection
            V     - Velocity 
        """
        Wing = self.Wing
        VTail = self.VTail
        NoseGear = self.NoseGear
        MainGear = self.MainGear
        Fuselage = self.Fuselage

        if V is None:
            V = self.GetV_LO()
               
        if len( self.MakeList(V) ) > 1 and len( self.MakeList(a2dw) ) > 1:
            raise ACAircraftError("Cannot compute drag from a list of angles of attack and a list of velocities")
        #
        # Create a list of all drag components
        #
        CD = []
        CDLegend = []

        Sw = Wing.S
        CDw = Wing.CD(a2dw)
        CD.append(CDw)
        CDLegend.append('Wing')
                        
        Sv  = VTail.S
        CDv = VTail.CD(a2dw*0, 0*ARCDEG) #No rudder deflection for now
        CD.append((CDv*Sv/Sw) )
        CDLegend.append('Vertical Tail')
        
        Dng = NoseGear.CD(Wing,a2dw)
        Dmg = MainGear.CD(Wing,a2dw)
        Dfus = Fuselage.CD(a2dw, V, self.GetAlt_LO(), self.Wing)
        CD.append( (Dng+Dmg+Dfus)*self.DragKFactor/Sw )
        CDLegend.append('Fuselage')

        return CD, CDLegend        

#===============================================================================
    def _CD(self, a2dw, del_e, V = None):
        """
        Calculates the total area weighted CD 
        of the aircraft with the given elevator deflection
        
        Inputs:
            a2dw  - 2D angle of attack of the main wing
            del_e - Elevator deflection
            V     - Velocity
        """
        
        CD, CDLegen = self._CDComponents(a2dw, del_e, V)
        
        CD = sum(CD)
        
        #return CD[0] if len(CD) == 1 else CD
        return CD.real

#===============================================================================
    def _CMComponents(self, alpha2dw, del_e = 0*ARCDEG, Xcg = None, Reht = None, i = None):
        """
        CM = CLw*(Xcg-Xac)/macw + CMw + eta*Sht/Sw*CLht*(Xcg-Xht)/macw + CMht*Sht/Sw
        
        Inputs:
            alhpah2dw - Angle of attack the Wing
            del_e     - Elevator deflection
            Xcg       - Center of gravity
            Reht      - Reynolds number of the horizontal tail
            i         - Incidence angle of the horizontal tail
        """
        Wing  = self.Wing

        if Xcg is None:
            Xcg = self.Xcg()

        #
        # Create a list of all moment components
        #
        CM = []
        CMLegend = []
        
        Sw   = Wing.S
        CLw  = Wing.CL(alpha2dw)
        CMw  = Wing.CM(alpha2dw)
        Xac  = Wing.Xac()
        macw = Wing.MAC()

        CM.append( CLw*(Xcg-Xac)/macw + CMw )
        CMLegend.append('Wing')
                           
        return npy.array(CM), CMLegend
            
#===============================================================================
    def CL(self, a2dw):
        """
        Returns the total CL of the aircraft with a 0 deg elevator deflection

        Inputs"
            a2dw - 2D angle of attack of the main wing
        """
        return self._CL(a2dw, 0*ARCDEG)

#===============================================================================
    def CLTrim(self, a2dw, Xcg = None):
        """
        Calculates the total trimmed  CL of the aircraft with the elevator deflected
        to produce a zero moment on the aircraft

        Inputs"
            a2dw - 2D angle of attack of the main wing
            Xcg  - A center of gravity position
        """
        del_e = self.del_e_trim(a2dw, Xcg)
        return self._CL(a2dw, del_e)

#===============================================================================
    def CD(self, a2dw, V = None):
        """
        Returns the total CD of the aircraft with a 0 deg elevator deflection

        Inputs"
            a2dw - 2D angle of attack of the main wing
            V    - Velocity
       """
        return self._CD(a2dw, 0*ARCDEG, V)

#===============================================================================
    def CDTrim(self, a2dw, Xcg = None, V = None):
        """
        Calculates the total trimmed  CD of the aircraft with the elevator deflected
        to produce a zero moment on the aircraft

        Inputs"
            a2dw - 2D angle of attack of the main wing
            Xcg  - A center of gravity position
            V    - Velocity
        """
        del_e = self.del_e_trim(a2dw, Xcg)
        return self._CD(a2dw, del_e, V)

#===============================================================================
    def WingDownWash(self, Alpha2dw):
        """
        Returns the downwas produced by the wing including the downwash factor on the tail
        
        Inputs:
            Alpha2dw - Two dimensional main wing angle of attack
        """
        if self.dirty: self.Refresh()
        
        return self.Wing.DownWash(Alpha2dw)

#===============================================================================
    def Thrust(self, V):
        """
        Returns the thrust produced by the propulsino system
        
        Inputs:
            V - Velocity of he aircraft
        """
        if self.dirty: self.Refresh()
        
        return self.Propulsion.T(V)

#===============================================================================
    def GetS(self):
        """
        Returns the area of the wing.
        """
        if self.dirty: self.Refresh()
        
        return self.Wing.S

#===============================================================================
    def Getb(self):
        """
        Returns the span of the wing.
        """
        if self.dirty: self.Refresh()
        
        return self.Wing.b
    
#===============================================================================
    def Getq_LO(self):
        """
        Returns the dynamic pressure.
        """
        if self.dirty: self.Refresh()
        
        h = self.Wing.GetAlt_LO()
        v = self.Wing.GetV_LO()
        return AeroUtil.q(h,v)
    
#===============================================================================
    def GetMAC(self):
        """
        Returns the Mean Aerodynamic Chord of the wing.
        """
        if self.dirty: self.Refresh()
        
        return self.Wing.MAC()

#===============================================================================
    def GetV_LO(self):
        """
        Returns the lift of velocity of the aircraft at take off conditions.
        """
        if self.dirty: self.Refresh()
        
        return self.param.V_LO

#===============================================================================
    def GetV_Stall(self):
        """
        Returns the lift of velocity of the aircraft at take off conditions.
        """
        if self.dirty: self.Refresh()
        
        return self.Wing.V_Stall
    
#===============================================================================
    def GetAlt_LO(self):
        """
        Returns the lift of altitude of the aircraft at take off conditions.
        """
        if self.dirty: self.Refresh()
        
        return self.Wing.GetAlt_LO()

#===============================================================================
    def GetAlpha2d_LO(self):
        """
        Returns the 2d angle of attack of the aircraft at take of conditions.
        """
        if self.dirty: self.Refresh()
        return self.Wing.GetAlpha2d_LO()

#===============================================================================
    def GetAlphaFus_LO(self):
        """
        Returns the 3d angle of attack of the aircraft at take of conditions.
        """
        if self.dirty: self.Refresh()
        return self.Wing.GetAlpha3d_LO()
        
#===============================================================================
    def GetCD_LO(self):
        """
        Returns the drag coefficient of the aircraft at take of conditions.
        """
        if self.dirty: self.Refresh()
        
        Alpha2d_LO = self.GetAlpha2d_LO()
        return self.CDTrim(Alpha2d_LO)[0]

#===============================================================================
    def GetCL_LO(self):
        """
        Returns the lift coefficient of the aircraft at take of conditions.
        """
        if self.dirty: self.Refresh()
        
        Alpha2d_LO = self.GetAlpha2d_LO()
        return self.CLTrim(Alpha2d_LO)
        
#===============================================================================
    def GetHT_Design_CL(self):
        """
        Returns the lift coefficient of the horizonatl tail at zero CM. CL is relative to the tail area.
        """
        if self.dirty: self.Refresh()
        
        return 0

 #===============================================================================
    def GetHT_Design_CD(self):
        """
        Returns the drag coefficient of the horizonatl tail at zero CM. CD is relative to the tail area.
        """
        if self.dirty: self.Refresh()
        
        return 0
    
#===============================================================================
    def AlphaRange(self, nalpha = 30):
        """
        Returns a range of angles of 2D angles of attack.
        """
        return self.Wing.AlphaRange(nalpha)
        
#===============================================================================
    def GetControlSurfaces(self):
        """
        Returns the names of all the control surfaces of this aircraft.
        """
        ControlNames = []        
        ControlNames += self._GetControlNames(self.Wing)
        ControlNames += self._GetControlNames(self.VTail)

        return ControlNames    

#===============================================================================
    def PlotCMPolars(self, fig = 1, del_e = [0*ARCDEG], XcgOffsets = [-0.5*IN, 0.5*IN]):
        """
        Plots CM with the CG at the neutral point, the nominal point, and some offset points

        Inputs:
            fig        - Figure number to draw the plots on
            del_e      - A list of elevator deflection
            XcgOffsets - A list of CG offsets to display the range of stability
        """
        alphas = self.Wing.AlphaRange()

        Xnp = self.Xnp()
        Xcg = self.Xcg()

        CMnp  = self.CM(alphas, del_e = 0 * ARCDEG, Xcg = Xnp)
        CMnom = self.CM(alphas, del_e = 0 * ARCDEG, Xcg = Xcg)
        del_e_trim = self.del_e_trim(alphas)

        CMs = []
        legend = ['$X_{np}$','Nominal']
        for offset in XcgOffsets:
            CMs.append(self.CM(alphas, del_e = 0 * ARCDEG, Xcg = Xcg+offset))
            sgn = '+' if offset > 0*IN else ' '
            legend.append('$X_{cg}$ ' + sgn + AsUnit(offset, 'in'))

        #
        # Translate the angles of attack to be relative to the fuselage
        #
        alphafus = self.Wing.AlphaFus(alphas) / (ARCDEG)

        #
        # Arguments for the grid of the plots
        #
        gridarg = {'color':'k', 'linestyle':'--', 'linewidth':1}

        pyl.figure(fig)
        pyl.subplot(131)
        pyl.grid(gridarg)

        pyl.plot(alphafus,CMnp)
        pyl.plot(alphafus,CMnom)
        for i in range(len(XcgOffsets)):
            pyl.plot(alphafus,CMs[i])

        ylabel_offset = -0.1

        pyl.title(r'$C_M$ Polar with CG Shift')
        pyl.legend(legend, loc = 'best', frameon = False)
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel('$C_M$'); pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)

        #
        # Plot the effect of elevator deflections
        #
        CMs = []
        legend = []
        for de in del_e:
            CMs.append(self.CM(alphas, del_e = de, Xcg = Xcg))
            sgn = '+' if de > 0*ARCDEG else ' '
            legend.append('$\delta_e' + sgn + str(de / (ARCDEG)) + '^o$')

        pyl.subplot(132)
        pyl.grid(gridarg)

        for i in range(len(del_e)):
            pyl.plot(alphafus,CMs[i])

        pyl.title(r'$C_M$ Polar with Elevator Deflection')
        pyl.legend(legend, loc = 'best', frameon = False)
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel('$C_M$'); pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
 
 
        #
        # Plot the trimmed elevator deflections
        #
        pyl.subplot(133)
        pyl.plot(alphafus, del_e_trim / (ARCDEG))
        
        pyl.title('Trimmed Elevator Deflection')
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel(r'Trimmed Elevator Deflection $\delta_{e}$'); pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
        pyl.grid(gridarg)
           
#===============================================================================
    def PlotPolarsSlopes(self, fig = 1):
        """
        Plots polars along with the approximated slopes

        Inputs:
            fig - The figure to plot to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        #
        # Create the angles of attack
        #
        alpha = self.Wing.AlphaRange()

        #
        # Create the polars.
        #
        CL   = self.CL(alpha)
        CM   = self.CM(alpha)
        DW   = self.WingDownWash(alpha) / (ARCDEG)
        
        dCL_da   = self.dCL_da() / (1/ARCDEG)
        dCM_da   = self.dCM_da() / (1/ARCDEG)
        dDW_da   = self.dDW_da() #This is already a number
        
        #
        # Get the 2D wing angles of attack
        #
#        ClSlopeAt   = (self.param.ClSlopeAt[1] + self.param.ClSlopeAt[0])/2
#        CmSlopeAt   = (self.param.CmSlopeAt[1] + self.param.CmSlopeAt[0])/2
#        dwSlopeAt   = (self.param.dwSlopeAt[1] + self.param.dwSlopeAt[0])/2
#        ClhtSlopeAt = (self.param.ClhtSlopeAt[1] + self.param.ClhtSlopeAt[0])/2
#
#        CLSlopeAt   = (self.CLSlopeAt[1] + self.CLSlopeAt[0])/2
#        CMSlopeAt   = (self.CMSlopeAt[1] + self.CMSlopeAt[0])/2
#        CLHTSlopeAt = (self.CLHTSlopeAt[1] + self.CLHTSlopeAt[0])/2
#        DWSlopeAt   = (self.DWSlopeAt[1] + self.DWSlopeAt[0])/2
        
        ClSlopeAt   = self.param.ClSlopeAt[0]
        CmSlopeAt   = self.param.CmSlopeAt[0]
        dwSlopeAt   = self.param.dwSlopeAt[0]

        CLSlopeAt   = self.CLSlopeAt[0]
        CMSlopeAt   = self.CMSlopeAt[0]
        DWSlopeAt   = self.DWSlopeAt[0]
        
        #
        # Make the slope go through the specified CLSlopeAt.
        #
        CLOffset   = self.CLTrim(ClSlopeAt)        - dCL_da  *CLSlopeAt / (ARCDEG)
        CMOffset   = self.CM(CmSlopeAt)            - dCM_da  *CMSlopeAt / (ARCDEG)
        DWOffset   = self.WingDownWash(dwSlopeAt) / (ARCDEG) - dDW_da*DWSlopeAt / (ARCDEG)
        
        alpha = self.Wing.AlphaFus(alpha) / (ARCDEG)
        #
        # These are in degrees without being a _scalar
        #
        alp_da = npy.array([min(alpha), max(alpha)])
        
        CL_da    = alp_da*dCL_da   + CLOffset
        CM_da    = alp_da*dCM_da   + CMOffset
        DW_da    = alp_da*dDW_da   + DWOffset

        #
        # Plot them
        #
        pyl.subplot(311)
        pyl.plot(alpha, CL)
        pyl.plot(alp_da, CL_da)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_l$')
        pyl.legend([r'$C_l$',r'$C_l$ Slope'], loc = 'best')


        pyl.subplot(312)
        pyl.plot(alpha, DW)
        pyl.plot(alp_da, DW_da)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'DownWash')
        pyl.legend([r'Down Wash',r'Down Wash Slope'], loc = 'best')

        pyl.subplot(313)
        pyl.plot(alpha, CM)
        pyl.plot(alp_da, CM_da)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_m$')
        pyl.legend([r'$C_m$',r'$C_m$ Slope'], loc = 'best')


        pyl.annotate("Aircraft Slope Approximations", xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)   

#===============================================================================
    def PlotTrimmedPolars(self, fig = 1):
        """
        Plots the lift and drag polars for the aircraft

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        #
        # Create the angles of attack
        #
        alpha = self.Wing.AlphaRange()

        #
        # Create the polars.
        #
        CL = self.CL(alpha)
        CD = self.CD(alpha)
        LoD = CL/CD

        CLWing = self.Wing.CL(alpha)
        CDWing = self.Wing.CD(alpha)
        CMWing = self.Wing.CM(alpha)
        LoDWing = CLWing/CDWing

        
        CLTrim = self.CLTrim(alpha)
        CDTrim = self.CDTrim(alpha)
        CM     = self.CM(alpha)
        
        LoDTrim = CLTrim/CDTrim
        
        
        alpha_LO = self.GetAlpha2d_LO()
        CL_LO = self.CLTrim(alpha_LO)
        CD_LO = self.CDTrim(alpha_LO)
        CM_LO = self.CM(alpha_LO)

        alpha_LoD_max = alpha[LoDTrim.argmax()]
        CL_LoD_max = self.CLTrim(alpha_LoD_max)
        CD_LoD_max = self.CDTrim(alpha_LoD_max)
        CM_LoD_max = self.CM(alpha_LoD_max)

        alpha_stall = alpha[CLTrim.argmax()]
        CL_stall = CLTrim.max()
        CD_stall = self.CDTrim(alpha_stall)
        CM_stall = self.CM(alpha_stall)
        
        CL_max_climb = CL_stall/self.Wing.V_LOstall**2
        
        alpha_max_climb = alpha[CLTrim.searchsorted(CL_max_climb, side='left')]
        CL_max_climb = self.CLTrim(alpha_max_climb)
        CD_max_climb = self.CDTrim(alpha_max_climb)
        CM_max_climb = self.CM(alpha_max_climb)
        
        alpha           = self.Wing.AlphaFus(alpha) / (ARCDEG)
        alpha_LoD_max   = self.Wing.AlphaFus(alpha_LoD_max) / (ARCDEG)
        alpha_stall     = self.Wing.AlphaFus(alpha_stall) / (ARCDEG)
        alpha_LO        = self.Wing.AlphaFus(alpha_LO) / (ARCDEG)
        alpha_max_climb = self.Wing.AlphaFus(alpha_max_climb) / (ARCDEG)
        
        format_Total = 'b'
        format_Trim = 'r'
        format_Wing = 'k--'
        
        pointformat_LO = 'bo'
        pointformat_stall = 'ro'
        pointformat_LoD = 'go'
        pointformat_climb = 'ko'
        
        gridarg = {'color':'k', 'linestyle':'--', 'linewidth':1}
        #
        # Draw them
        #
        xlabel_offset = -0.1
        ylabel_offset = -0.11
        
        pyl.subplot(321)
        line1, line2, line3 = pyl.plot( alpha, CL, format_Total
                                      , alpha, CLTrim, format_Trim
                                      , alpha, CLWing, format_Wing)

        pyl.figlegend((line1, line2, line3), ('Total','Total Trim','Wing'), 'upper center', ncol=3, frameon = False)
                
        point_LO, point_stall, point_LoD, point_climb = pyl.plot( alpha_LO, CL_LO, pointformat_LO
                                                                , alpha_stall, CL_stall, pointformat_stall
                                                                , alpha_LoD_max, CL_LoD_max, pointformat_LoD
                                                                , alpha_max_climb, CL_max_climb, pointformat_climb )

        pyl.figlegend( (point_LO, point_stall, point_LoD, point_climb)
                     , ['LO','Stall','LoD','Max Climb'], loc = 'lower center', ncol=4, numpoints = 1, frameon = False)
        
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_L$')
        pyl.grid(gridarg)
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)

        pyl.subplot(322)
        pyl.plot(alpha, CD, format_Total)
        pyl.plot(alpha, CDTrim, format_Trim)
        pyl.plot(alpha, CDWing, format_Wing)
        pyl.plot(alpha_LoD_max, CD_LoD_max, pointformat_LoD )
        pyl.plot(alpha_stall, CD_stall, pointformat_stall )
        pyl.plot(alpha_LO, CD_LO, pointformat_LO )
        pyl.plot(alpha_max_climb, CD_max_climb, pointformat_climb )
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_D$')
        pyl.grid(gridarg)
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
 
        pyl.subplot(323)
        pyl.plot(CDTrim, CL, format_Total)
        pyl.plot(CDTrim, CLTrim, format_Trim)
        pyl.plot(CDTrim, CLWing, format_Wing)
        pyl.plot(CD_LoD_max, CL_LoD_max, pointformat_LoD )
        pyl.plot(CD_stall, CL_stall, pointformat_stall )
        pyl.plot(CD_LO, CL_LO, pointformat_LO )
        pyl.plot(CD_max_climb, CL_max_climb, pointformat_climb )
        pyl.xlabel(r'$C_D$'); pyl.ylabel(r'$C_L$')
        pyl.grid(gridarg)
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)

        pyl.subplot(324)
        pyl.plot(alpha, LoD, format_Total)
        pyl.plot(alpha, LoDTrim, format_Trim)
        pyl.plot(alpha, LoDWing, format_Wing)
        pyl.plot(alpha_LoD_max, CL_LoD_max/CD_LoD_max, pointformat_LoD )
        pyl.plot(alpha_stall, CL_stall/CD_stall, pointformat_stall )
        pyl.plot(alpha_LO, CL_LO/CD_LO, pointformat_LO )
        pyl.plot(alpha_max_climb, CL_max_climb/CD_max_climb, pointformat_climb )
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel('L/D')
        pyl.grid(gridarg)
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)

        pyl.subplot(325)
        pyl.plot(alpha, CM, format_Total)
        #pyl.plot(alpha, CMWing, format_Wing)
        pyl.plot(alpha_LoD_max, CM_LoD_max, pointformat_LoD )
        pyl.plot(alpha_stall, CM_stall, pointformat_stall )
        pyl.plot(alpha_LO, CM_LO, pointformat_LO )
        pyl.plot(alpha_max_climb, CM_max_climb, pointformat_climb )
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.grid(gridarg)
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)

#        pyl.annotate("Aircraft Trimmed Polars", xy=(.025, .975),
#                xycoords='figure fraction',
#                horizontalalignment='left', verticalalignment='top',
#                fontsize=20)           

        #
        # Add a table to the sixth subplot
        #
        Form = '%1.1f'
        Form2 = '%1.2f'
        colLabels = [r'$\alpha^o$',r'$C_L$',r'$C_D$',r'$C_M$',r'L/D']
        rowLabels = [r'LO', 
                     r'Stall', 
                     r'L/D max', 
                     r'Max Climb']
        
        cellText  = [[ Form % alpha_LO       , Form % CL_LO       , Form2 % CD_LO       , Form2 % CM_LO       , Form % (CL_LO/CD_LO) ],
                     [ Form % alpha_stall    , Form % CL_stall    , Form2 % CD_stall    , Form2 % CM_stall    , Form % (CL_stall/CD_stall) ],
                     [ Form % alpha_LoD_max  , Form % CL_LoD_max  , Form2 % CD_LoD_max  , Form2 % CM_LoD_max  , Form % (CL_LoD_max/CD_LoD_max) ],
                     [ Form % alpha_max_climb, Form % CL_max_climb, Form2 % CD_max_climb, Form2 % CM_max_climb, Form % (CL_max_climb/CD_max_climb) ]]
                
        colWidths = self._GetTableColWidths(rowLabels, cellText)
        
        pyl.subplot(326)
#        pyl.title("Aircraft Properties")
        pyl.axis('off')
        table = pyl.table(cellText  = cellText , cellLoc = 'center',
                          rowLabels = rowLabels, rowLoc  = 'left',
                          colWidths = colWidths,
                          colLabels = colLabels, colLoc  = 'center',
                          loc = 'upper right')

        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1., 1.7)

           
#===============================================================================
    def PlotDragBuildup(self, fig = 1, nalpha = 30):
        """
        Draws the drag components for the aircraft

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        #
        # Create the angles of attack
        #
        alpha = self.AlphaRange(nalpha)

        #
        # Create the polars.
        #
        del_e = self.del_e_trim(alpha)
        CD, CDLegend = self._CDComponents(alpha, del_e)
        
        CDLegend.append('Total Drag')
        
        alpha = self.Wing.AlphaFus(alpha) / (ARCDEG)
        #
        # Draw them
        #
        #pyl.subplot(111)
        for i in range(len(CD)):
            pyl.plot(alpha, CD[i])
        pyl.plot(alpha, sum(CD))

        gridarg = {'color':'k', 'linestyle':'--', 'linewidth':1}
        
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_D$')
        pyl.gca().xaxis.set_label_coords(0.5,-0.05)
        pyl.gca().yaxis.set_label_coords(-0.075,0.5)
        pyl.legend(CDLegend, loc = 'best', labelspacing = 0.0)
        pyl.grid(gridarg)

        pyl.annotate("Aircraft Drag Buildup", xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)

#===============================================================================
    def PlotCLCMComponents(self, fig = 1, del_es = [0*ARCDEG]):
        """
        Draws the components for the lift and moment polars of the aircraft

        Inputs:
            fig    - The figure to draw to
            del_es - A list of elevator deflections
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        Wing  = self.Wing
        Xcg = self.Xcg()
        
        #
        # Create the angles of attack
        #
        alpha2dw = self.AlphaRange()
        del_e_trim = self.del_e_trim(alpha2dw)

        #
        # Create a list of all moment components
        #
        CMs = []
        legend = []
        CLs = []
        
        Sw   = Wing.S
        CLw  = Wing.CL(alpha2dw)
        CMw  = Wing.CM(alpha2dw)
        Xac  = Wing.Xac()
        macw = Wing.MAC()

        CMs.append( CLw*(Xcg-Xac)/macw + CMw )
        CLs.append(CLw)
        legend.append('Wing')
        
        #
        # Moment for the tail
        #
        eta  = HTail.eta
        Sht  = HTail.S
        CMht = HTail.CM(alpha2dw)
        Xht  = HTail.X[0]

        def TailCM(CLht):
            return CMht*Sht/Sw + eta*Sht/Sw*CLht*(Xcg-Xht)/macw
        
        CLht = HTail.CL(alpha2dw, del_e_trim)
        CMs.append(TailCM(CLht))
        CLs.append(CLht)
        legend.append('HTail del_e = trimmed')

        for del_e in del_es:
            CLht = HTail.CL(alpha2dw, del_e)
            CLs.append(CLht)
            CMs.append(TailCM(CLht))
            legend.append('HTail del_e = ' + str( AsUnit(del_e, 'deg') ) )
        
        #
        # Draw them
        #
        alphafus = self.Wing.AlphaFus(alpha2dw) / (ARCDEG)
        
        pyl.subplot(131)
        for i in range(len(CLs)):
            pyl.plot(alphafus, CLs[i])

        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_L$')
        pyl.legend(legend, loc = 'best', labelspacing = 0.0)

        pyl.subplot(132)
        for i in range(len(CMs)):
            pyl.plot(alphafus, CMs[i])

        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.legend(legend, loc = 'best', labelspacing = 0.0)

        pyl.subplot(133)
        for i in range(len(del_es)):
#            if del_es[i] == 0*ARCDEG
            pyl.plot(alphafus, CMs[i])

        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.legend(legend, loc = 'best', labelspacing = 0.0)

        pyl.annotate("Aircraft CL and CM Buildup", xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)
        
#===============================================================================
    def PlotTailLoad(self, fig = 1, del_e_min = -10 * ARCDEG, del_e_max = 10 * ARCDEG, Vmax = None ):
        """
        Draws the lift polars for the horizontal tail

        Inputs:
            fig - The figure to draw to
            del_e_min - Minimum elevator deflection
            del_e_min - Maximum elevator deflection
            Vmax      - Maximum velocity to calculate the tail load
        """
        
        if self.dirty: self.Refresh()
        
        if Vmax is None: Vmax = self.Vmax()
        
        pyl.figure(fig)

        h = self.GetAlt_LO()
        q = AeroUtil.q(h,Vmax)
        S = self.GetS()
        qS = q*S / LBF
        
        #
        # Get the trimmed load on the horizontal tail
        #
        V = self.VelocityRange(0*FT/SEC, self.VmaxPlt)
        LiftComp, leg = self._LiftComponents(V)
        Lifth = LiftComp[1] / LBF
        
        DragComp, leg = self._DragComponents(V)
        Dragh = DragComp[1] / LBF
 
        #
        # Create the angles of attack
        #
        alpha = self.Wing.AlphaRange()

        #
        # Create the polars.
        #
        del_e_trim = self.del_e_trim(alpha)
        CLh_min  = self.HTail.CL(alpha, del_e = del_e_min )
        CLh_trim = self.HTail.CL(alpha, del_e = del_e_trim)
        CLh_max  = self.HTail.CL(alpha, del_e = del_e_max )
        
        CDh_trim = self.HTail.CD(alpha, del_e = del_e_trim)
                
        #
        # CL/CD of the tail at zero CM
        #
        CLhi = self.GetHT_Design_CL()
        CDhi = self.GetHT_Design_CD()
        Alpha_Zero_CM = self.Alpha_Zero_CM / ARCDEG
                
        alpha = self.Wing.AlphaFus(alpha) / ARCDEG
        V = V / (FT/SEC)
        #
        # Plot them
        #
        pyl.subplot(321)
#        pyl.plot(alpha, CLh_min )
        pyl.plot(alpha, CLh_trim)
#        pyl.plot(alpha, CLh_max )
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel(r'Horizontal Tail $C_L$')
        pyl.plot([Alpha_Zero_CM], [CLhi], 'o')
        pyl.legend(['CL', 'Zero CM' ], loc = 'best', numpoints=1)

        pyl.subplot(322)
        pyl.plot(alpha, CDh_trim)
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel(r'Horizontal Tail $C_D$')
        pyl.plot([Alpha_Zero_CM], [CDhi], 'o')
        pyl.legend(['CD', 'Zero CM' ], loc = 'best', numpoints=1)

        
        pyl.subplot(323)
        pyl.plot(V, Lifth)
        pyl.xlabel(r'Velocity (ft/sec)'); pyl.ylabel(r'Trimmed Horizontal Tail Load (lbf)')

        pyl.subplot(324)
        pyl.plot(V, Dragh)
        pyl.xlabel(r'Velocity (ft/sec)'); pyl.ylabel(r'Trimmed Horizontal Tail Drag (lbf)')

        pyl.subplot(325)
        pyl.plot(alpha, del_e_trim / ARCDEG)
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel(r'Trimmed Elevator Deflection $\delta_{e}$')

        pyl.subplot(326)
        pyl.plot(alpha, CLh_min  * qS)
        pyl.plot(alpha, CLh_trim * qS)
        pyl.plot(alpha, CLh_max  * qS)
        pyl.legend(['$\delta_{e}$ = ' + AsUnit(del_e_min, 'deg', '%2.0f'), '$\delta_{e}$ Trim', '$\delta_{e}$ = ' + AsUnit(del_e_max, 'deg', '%2.0f') ], loc = 'best')
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel(r'Horizontal Tail Load (lbf)  (V = ' + AsUnit(Vmax, 'ft/s', '%3.0f') + ')' )
 
        pyl.annotate("Horizontal Tail Charachteristics", xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)   
   
#===============================================================================
    def PlotVTailLoad(self, fig = 1):
        """
        Draws the force on the vertical tail for varying delta rudder

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        #
        # Find the CL's
        #
        del_r = npy.linspace(5,30,6)*ARCDEG
        V     = self.VelocityRange(0*FT/SEC, self.VmaxPlt)
        Vplt  = V / (FT/SEC)
        
        #
        # Calculate lift for the varying rudder deflections
        #
        leg = []
        for r in del_r:
            CLv   = self.VTail.CL(0*ARCDEG, r)
            h     = self.GetAlt_LO()
            q     = AeroUtil.q(h,V)
            S     = self.VTail.S 
            Lv    = (CLv*q*S) / (LBF)
            
            leg.append(r'$\delta_r$ = %2.0f' % (r / ARCDEG))
            pyl.plot(Vplt, Lv)
            
        #
        # Plot them
        #
        pyl.title('Vertical Tail Loads')
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Vertical Tail Force (lbf)')
        pyl.legend(leg, loc = 'best')   

#===============================================================================
    def PlotHTailLoad(self, fig = 1):
        """
        Draws the force on the vertical tail for varying delta rudder

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)

        #
        # Find the CL's
        #
        del_e = npy.linspace(-25,25,6)*ARCDEG
        V     = self.VelocityRange(0*FT/SEC, self.VmaxPlt)
        Vplt  = V / (FT/SEC)
        
        #
        # Calculate lift for the varying rudder deflections
        #
        leg = []
        for e in del_e:
            CL = self.HTail.CL(0*ARCDEG, e)
            h  = self.GetAlt_LO()
            q  = AeroUtil.q(h,V)
            S  = self.HTail.S 
            L  = (CL*q*S) / (LBF)
            
            leg.append(r'$\delta_e$ = %2.0f' % (e / ARCDEG))
            pyl.plot(Vplt, L)
        
        #
        # Get lift off tail load 
        #
        alpha      = self.GetAlphaFus_LO()
        del_e_trim = self.del_e_trim(alpha)
        
        CL = self.HTail.CL(alpha, del_e_trim)
        h  = self.GetAlt_LO()
        q  = AeroUtil.q(h,self.GetV_LO())
        S  = self.HTail.S 
        L  = (CL*q*S) / (LBF)
        
        CL = self.HTail.CL(alpha, -25*ARCDEG)
        L2  = (CL*q*S) / (LBF)
        
        CL = self.HTail.CL(alpha, 25*ARCDEG)
        L3  = (CL*q*S) / (LBF)
        
        print
        print 'Trim lift off horizontal tail force: ', L
        print
        print 'Lift off horizontal tail force at -25 degree deflection: ', L2
        print
        print 'Lift off horizontal tail force at 25 degree deflection: ', L3
        print
        
        #
        # Plot them
        #
        pyl.title('Horizontal Tail Loads')
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Horizontal Tail Force (lbf)')
        pyl.legend(leg, loc = 'best')  
