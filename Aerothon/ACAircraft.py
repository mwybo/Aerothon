"""
A collection of classes for aircrafts
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import M, FT, IN, ARCDEG, RAD, LBF, SEC, LBM, MIN
from scalar.units import AsUnit
from ACBase import ACComponent, Angle, Unitless, Length, Force, Time, Velocity, g
from ACTableBase import ACTableBase
from ACWing import ACMainWing
from ACTail import  ACTailSurf, ACHorizontalTailSurf
from ACLandingGear import ACLandingGear
from ACFuselage import ACFuselage
from ACMotor import ACMotor
#from ACPropulsion import ACPropulsion
from ACWeightTable import ACWeightTable
#import pylab as pyl
import matplotlib.pyplot as pyl
import numpy as npy
#from RootSolvers import SecantSolver, SecantError, BisectionSolver
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import AeroUtil
import math


################################################################################
class ACAircraftError(Exception):
    """
    Error handling to the ACAircraft
    """
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg
    def __repr__(self):
        return self.msg

################################################################################
class ACAircraft(ACComponent, ACTableBase):
    """
    A class base class for constructing aircraft assembling classes
    
    Attributes:
        Xnp_Init         - An initial guess for the neutral point
        EmptyWeight      - The empty weight of the aircraft
        TotalWeight      - The total weight of the aircraft
        CLSlopeAt        - Fuselage angle of attack used to evaluate slope of CL
        CMSlopeAt        - Fuselage angle of attack used to evaluate slope of CL
        StaticMargin     - Percent of MAC of the wing which the CG is placed in front of the neutral point
        Alpha_Zero_CM    - Fuselage angle of attack that will have a zero CM (used to set the horizontal tail incidence angle)
        Mu_r             - Rolling friction coefficient
        WingFuseFrac     - Vertical fuselage fraction used to position the wing (0 - Bottom, 1 - Top)
        Alpha_Groundroll - The ground roll angle of the fuselage
        TippingAngle     - The angle after which the CG of the aircraft moves behind the main landing gear
        RotationTime     - Estimated time for the aircraft to rotate on the runway during takeoff (default 0.5 sec)

        DragKFactor      - Multiplier to increase the drag of the fuselage and landing gear (default 1.15)

        VmaxPlt          - Maximum velocity for plotting purposes (default 100 FT/SEC)

        iZsym            - (True/False) Turns ground effect in AVL on or off (default False)
        Zsym             - Height of the AVL aircraft above the ground (default 0.0)

    """

#===============================================================================
    def __init__(self):
        super(ACAircraft, self).__init__()
        
        UnitList = self.UnitList
        
        self.NoneList['Xnp_Init']    = None   ; UnitList['Xnp_Init']    = Length
        self.NoneList['EmptyWeight'] = None   ; UnitList['EmptyWeight'] = Force
        self.NoneList['TotalWeight'] = None   ; UnitList['TotalWeight'] = Force
        self.NoneList['CLSlopeAt']   = None   ; UnitList['CLSlopeAt']   = Angle
        self.NoneList['CMSlopeAt']   = None   ; UnitList['CMSlopeAt']   = Angle
 
        self.__dict__['StaticMargin']     = 0.1          ; UnitList['StaticMargin']     = Unitless
        self.__dict__['Alpha_Zero_CM']    = 5 * ARCDEG   ; UnitList['Alpha_Zero_CM']    = Angle
        self.__dict__['Mu_r']             = 0.04         ; UnitList['Mu_r']             = Unitless
        self.__dict__['WingFuseFrac']     = 1            ; UnitList['WingFuseFrac']     = Unitless
        self.__dict__['Alpha_Groundroll'] = 0 * ARCDEG   ; UnitList['Alpha_Groundroll'] = Angle
        self.__dict__['TippingAngle']     = 0 * ARCDEG   ; UnitList['TippingAngle']     = Angle
        self.__dict__['RotationTime']     = 0.5 * SEC    ; UnitList['RotationTime']     = Time
        
        self.__dict__['DragKFactor']      = 1.15         ; UnitList['DragKFactor']      = Unitless

        self.__dict__['VmaxPlt']          = 100 * FT/SEC ; UnitList['VmaxPlt']          = Velocity
        
        self.__dict__['iZsym']            = False        ; UnitList['iZsym']            = Unitless
        self.__dict__['Zsym']             = 0.0*IN       ; UnitList['Zsym']             = Length
        #
        # Used for lift of velocity calculation
        #
        self.param.V_LO = None
        
#===============================================================================
    def _WriteAVLSurfaces(self, fAVL):
        """
        Writes out the AVL surface of all lifting surfaces of the aircraft

        Inputs:
            fAVL - The open file to write to
        """
        for part in self.__dict__.itervalues():
            if hasattr(part, '_WriteAVLSurface'):
                part._WriteAVLSurface(fAVL)

#===============================================================================
    def _WriteAVLHeader(self, fAVL, Xcg = None):
        """
        Writes out an AVL input file

        Inputs:
            fAVL - The AVL file to be written
        """
        
        CG = self.DesignCG()
            
        mach = 0. #Mach number is just hard coded to zero for now...

        iYsym = 0                       # This is an integer
        iZsym = 1 if self.iZsym else 0  # This is an integer
        Zsym  = float(self.Zsym/IN)     # This is a floating point

        Sref = self.Wing.S / (IN**2)
        Cref = self.Wing.MAC() / (IN)
        Bref = self.Wing.b / (IN)

        Xref = Xcg / IN if Xcg is not None else self.Xcg() / (IN)
        Yref = CG[1] / (IN)
        Zref = CG[2] / (IN)

        CDp = 0.

        fAVL.write(self.name + '\n')
        fAVL.write(str(mach)  + '\t\t\t| Mach\n')
        fAVL.write(str(iYsym) + ' ' +  str(iZsym) + ' ' +  str(Zsym) + '\t\t| iYsym  iZsym  Zsym\n')
        fAVL.write(str(Sref)  + ' ' +  str(Cref)  + ' ' +  str(Bref) + '\t\t| Sref   Cref   Bref\n')
        fAVL.write(str(Xref)  + ' ' +  str(Yref)  + ' ' +  str(Zref) + '\t\t| Xref   Yref   Zref\n')
        fAVL.write(str(CDp)   +                                        '\t\t\t| CDp  (optional)\n')
        
#===============================================================================
    def WriteAVLAircraft(self, filename, Xcg = None):
        """
        Writes out an AVL input file for the entire aircraft

        Inputs:
            filename - The name of the AVL file to be written
        """
        fAVL = open(filename, 'w')

        self._WriteAVLHeader(fAVL, Xcg)
        self._WriteAVLSurfaces(fAVL)

        fAVL.close()
        
                
#===============================================================================
    def EmptyCG(self):
        """
        Calculates the empty CG of the aircraft
        """
        if self.dirty:  self.Refresh()
        
        EmptyWeight = self.EmptyWeight
        CG = npy.array([0, 0, 0]) * IN*LBF

        for part in self.__dict__.itervalues():
            if hasattr(part, 'CG') and hasattr(part, 'Weight'):
                CG += part.CG()*part.Weight

        return CG/EmptyWeight

#===============================================================================
    def MOI(self):
        """
        Returns the moments of inertia for the aircraft.
        """
        if self.dirty: self.Refresh()
        
        return self.param.MOI

#===============================================================================
    def _UpdateMOI(self):
        """
        Updates the moments of inertia of the aircraft
        """
        CG = self.EmptyCG()
        
        MOI = npy.array([0,0,0])*LBM*IN**2
        
        def AddMOI(PCG, PM, PMOI):
            dCG = PCG - CG
            MOI[0] += PMOI[0] + PM*(dCG[1]**2 + dCG[2]**2) #Ixx
            MOI[1] += PMOI[1] + PM*(dCG[0]**2 + dCG[2]**2) #Iyy
            MOI[2] += PMOI[2] + PM*(dCG[1]**2 + dCG[0]**2) #Izz
             
        for part in self.__dict__.itervalues():
            if hasattr(part, 'CG') and hasattr(part, 'Weight') and hasattr(part, 'MOI'):
                AddMOI(part.CG(), part.Weight / g, part.MOI() )
        
        self.param.MOI = MOI

#===============================================================================
    def GetWeightTable(self):
        """
        Creates a weight table with all parts of the aircraft
        """
        if self.dirty: self.Refresh()
        
        CG = self.EmptyCG()
        
        WeightTable = ACWeightTable(CG)
                     
        for key, part in self.__dict__.iteritems():
            if hasattr(part, 'AddToWeightTable'):
                part.AddToWeightTable(key, WeightTable)
        
        return WeightTable
        
#===============================================================================
    def PayloadWeight(self):
        """
        Calculates the payload weight for the aircraft
        """
        if self.dirty: self.Refresh()
        
        return self.TotalWeight - self.EmptyWeight

#===============================================================================
    def _UpdateEmptyWeight(self):
        """
        Calculates the empty weight of the aircraft
        """
        EmptyWeight = 0*LBF

        for part in self.__dict__.itervalues():
            if hasattr(part, 'Weight'):
                EmptyWeight += part.Weight
            
        self.param.EmptyWeight = EmptyWeight
    
#===============================================================================
    def _CalcEmptyWeight(self):
        """
        Returns the empty weight of the aircraft
        """
        if self.dirty: self.Refresh()
        
        return self.param.EmptyWeight

#===============================================================================
    def __getattr__(self, key):
        """
        Attempt to calculate any variables set to None
        """
        ans = super(ACAircraft, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'EmptyWeight':
            return self._CalcEmptyWeight()
        #
        # If no calculation exist return None
        #
        return None

#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """
        Draw all the parts of the aircraft
        """
        if self.dirty: self.Refresh()
        
        #
        # Draw the CG (Draw first so the legend can be made)
        #
        pyl.figure(fig,figsize=(10,10))
        def DrawCG(cg, color):
            cgx = (cg[0]) / (IN)
            cgy = (cg[1]) / (IN)
            cgz = (cg[2]) / (IN)

            #pyl.subplot(top)
            #pyl.plot([cgx],[cgy],color)
            
            pyl.subplot(side)
            pyl.plot([cgx],[cgz],color)

        DrawCG(self.EmptyCG(), 'ro')          #Draw a red circle for the empty CG
        DrawCG(self.TotalCG(), 'mo')          #Draw a magenta circle for the total CG
        DCG = self.DesignCG()
        DrawCG(DCG,'go')                      #Draw a green circle for the deisgn CG
        DrawCG([self.Xnp(),0*IN,DCG[2]],'yo') #Draw a yellow circle for the neutral point
        DrawCG([self.Wing.Xac(),0*IN,DCG[2]],'bo') #Draw a blue circle for the Wing AC
        
        pyl.subplot(side)
        pyl.legend(['Empty CG', 'Total CG','Design CG','Neutral', 'AC'],loc='best',numpoints=1)
    
        #
        # Now draw the rest
        # 
        super(ACAircraft, self).Draw(fig, top, side, front)
        
        #
        # Go through all the attributes and call draw if available
        #
        for part in self.__dict__.itervalues():
            if hasattr(part, 'Draw'):
                part.Draw(fig, top, side, front)
                
                
        #
        # Add a table to the fourth subplot with some pertinent information
        #
        ECG = self.EmptyCG()[0] / (IN)
        ACG = self.DesignCG()[0] / (IN)
        dECG = ECG - ACG
        if self.Fuselage.Payload:
            PCG = self.Fuselage.Payload.CG()[0] / (IN)
            dPCG = PCG - ACG
        

        AsUnit = '%1.2f'
#        colLabels = []
        rowLabels = [r'Total Weight', 
                     r'Empty Weight', 
                     r'Payload Weight', 
                     r'Ground Roll',
                     r'Takeoff AoA', 
                     r'L.O. Rate of Climb', 
                     r'Max Rate of Climb', 
                     r'Overturn Angle < 63',
                     r'Empty $\Delta $CG$_x$']
        
        if self.Fuselage.Payload:
            rowLabels.append( r'Payload $\Delta $CG$_x$' )

        cellText  = [[(AsUnit + ' (lbf)'   ) % (self.TotalWeight / LBF)],
                     [(AsUnit + ' (lbf)'   ) % (self.EmptyWeight / LBF)],
                     [(AsUnit + ' (lbf)'   ) % (self.PayloadWeight() / LBF)],
                     [(AsUnit + ' (ft)'    ) % (self.Groundroll() / FT)],
                     [(AsUnit + ' (deg)'   ) % (self.GetAlphaFus_LO() / ARCDEG)],
                     [(AsUnit + ' (ft/min)') % (self.Rate_of_Climb() / (FT/MIN))],
                     [(AsUnit + ' (ft/min)') % (self.Rate_of_Climb( self.V_max_climb() ) / (FT/MIN))],
                     [(AsUnit + ' (deg)'   ) % (self.OverturnAngle() / ARCDEG)],
                     [(AsUnit + ' (in)'    ) % dECG]]
        
        if self.Fuselage.Payload:
            cellText.append( [(AsUnit + ' (in)'    ) % dPCG] )
                
        colWidths = self._GetTableColWidths(rowLabels, cellText)
        
        pyl.subplot(224)
        pyl.title("Aircraft Properties")
        pyl.axis('off')
        table = pyl.table(cellText  = cellText , cellLoc = 'center',
                          rowLabels = rowLabels, rowLoc  = 'left',
                          colWidths = colWidths,
                          #colLabels = colLabels, colLoc  = 'center',
                          loc = 'upper right')

        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1., 1.7)
        
                
#===============================================================================
    def Refresh(self):
        """
        Updates all the internal quantities.
        """
        super(ACAircraft, self).Refresh()
        self.param.refreshing = True

        self._CheckConsistent()

        self.param.ClSlopeAt = self.Wing.Alpha2D(self.CLSlopeAt)
        self.param.CmSlopeAt = self.Wing.Alpha2D(self.CMSlopeAt)

        self.param.refreshing = False
#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that the aircraft is not over or under specified
        """
        
        # There must be a total weight
        self._CheckEquation(['TotalWeight'], Need = 1)
        
        self._CheckEquation(['CLSlopeAt'], Need = 1)
        self._CheckEquation(['CMSlopeAt'], Need = 1)
        
        if self.CLSlopeAt[1] <= self.CLSlopeAt[0]:
            raise ACAircraftError(self.name + ": CLSlopeAt[1] must be grater than CLSlopeAt[0]")

        if self.CMSlopeAt[1] <= self.CMSlopeAt[0]:
            raise ACAircraftError(self.name + ": CMSlopeAt[1] must be grater than CMSlopeAt[0]")
        

#===============================================================================
    def _GetControlNames(self, LiftSurf):
        """
        Returns the names of the controls surfaces of a lifting surface.
        """
        ControlsNames = []
        for control in LiftSurf.Controls.__dict__.itervalues():
            ControlsNames.append(control.name)
        
        return ControlsNames       

#===============================================================================
    def GetIxx(self):
        """
        Returns Ixx for the aircraft.
        """
        if self.dirty: self.Refresh()
        
#        return 0.119*SLUG*FT**2
        
        return self.param.MOI[0]

#===============================================================================
    def GetIyy(self):
        """
        Returns Iyy for the aircraft.
        """
        if self.dirty: self.Refresh()
    
#        return 0.441*SLUG*FT**2
    
        return self.param.MOI[1]

#===============================================================================
    def GetIzz(self):
        """
        Returns Izz for the aircraft.
        """
        if self.dirty: self.Refresh()
    
#        return 0.518*SLUG*FT**2
    
        return self.param.MOI[2]

#===============================================================================
    def Cmadot(self):
        """
        TODO: Document what Cmadot is.
        
        It is for controls and relies on a horizontal tail.
        """
        assert(0)

#===============================================================================
    def Czq(self):
        """
        TODO: Document what Czq is.
        
        It is for controls
        """
        assert(0)
    
#===============================================================================
    def GetControlSurfaces(self):
        """
        Returns the names of all the control surfaces of this aircraft.
        """
        
        return ('')     

#===============================================================================
    def Lift(self, V, Alpha2d = None):
        """
        Calculates the lift of the aircraft
        
        Inputs: 
            V       - Velocity
            Alpha2d - Two dimensional angle of attack of the wing
        """
        h = self.GetAlt_LO()
        q = AeroUtil.q(h,V)
        S = self.GetS()
        a2d = self.AlphaTrim(V) if Alpha2d is None else Alpha2d                
        return self.CLTrim(a2d)*S*q

#===============================================================================
    def _LiftComponents(self, V, Alpha2d = None):
        """
        Calculates the lift components of the aircraft
        
        Inputs: 
            V       - Velocity
            Alpha2d - Two dimensional angle of attack of the wing
        """
        h = self.GetAlt_LO()
        q = AeroUtil.q(h,V)
        S = self.GetS()
        a2d = self.AlphaTrim(V) if Alpha2d is None else Alpha2d
        
        del_e = self.del_e_trim(a2d)
        CL, CLLegend = self._CLComponents(a2d, del_e)
                       
        CLCom = []
        for cl in CL:
            CLCom.append(cl*S*q)
        
        return CLCom, CLLegend

#===============================================================================
    def _DragComponents(self, V, Alpha2d = None):
        """
        Calculates the drag components of the aircraft
        
        Inputs: 
            V       - Velocity
            Alpha2d - Two dimensional angle of attack of the wing
        """
        h = self.GetAlt_LO()
        q = AeroUtil.q(h,V)
        S = self.GetS()
        a2d = self.AlphaTrim(V) if Alpha2d is None else Alpha2d
        
        del_e = self.del_e_trim(a2d)
        CD, CDLegend = self._CDComponents(a2d, del_e)
                       
        CDCom = []
        for cd in CD:
            CDCom.append(cd*S*q)
        
        return CDCom, CDLegend
        
#===============================================================================
    def AerodynamicDrag(self, V, Alpha2d = None):
        """
        Calculates the drag of the aircraft
        
        Inputs: 
            Alpha2d - Two dimensional angle of attack of the wing
            V       - Velocity
        """
        h = self.GetAlt_LO()
        q = AeroUtil.q(h,V)
        S = self.GetS()
        a2d = self.AlphaTrim(V) if Alpha2d is None else Alpha2d
        return self.CDTrim(a2d, V = V)*S*q

#===============================================================================
    def RollingDrag(self, V, Alpha2d = None):
        """
        Calculates the drag of the aircraft due to rolling friction on the ground
        
        Inputs: 
            V       - Velocity
            Alpha2d - Two dimensional angle of attack of the wing
        """
        a2d = self.AlphaTrim(V) if Alpha2d is None else Alpha2d
        return (self.TotalWeight - self.Lift(V, a2d))*self.Mu_r
    
#===============================================================================
    def Drag(self, V, Alpha2d = None):
        """
        Calculates the total drag of the aircraft
        
        Inputs: 
            V       - Velocity
            Alpha2d - Two dimensional angle of attack of the wing
        """
        V_LO = self.GetV_LO()
        
        V = self.MakeList(V)
        
        Drag = []
        for v in V:
            if v >= V_LO:
                Drag.append( self.AerodynamicDrag(v, Alpha2d) / (LBF))
            else:
                a2d = self.param.A_GR if Alpha2d is None else Alpha2d
                Drag.append( (self.RollingDrag(v, a2d) + self.AerodynamicDrag(v, a2d)) / (LBF) )
                
        return Drag[0]*LBF  if len(Drag) == 1 else npy.array(Drag)*LBF

#===============================================================================
    def VelocityRange(self, Vmin = None, Vmax = None, nV = 30):
        """
        Generates a range of velocity values.
        
        Inputs:
            Vmin - Minimum velocity
            Vmax - Maximum velocity
            nV   - Number of velocity points
        """
        if self.dirty: self.Refresh()
        
        if Vmin is None:
            Vmin = 0*FT/SEC
            
        if Vmax is None:
            Vmax = self.VmaxPltce
        
        return npy.linspace(Vmin / (FT/SEC), Vmax / (FT/SEC), nV)*FT/SEC

#===============================================================================
    def _UpdateAlphaTrim(self):
        """
        Calculates an array of Alpha Trim for the aircraft that will be interpolated
        in the method AlphaTrim
        """
        Lift  = self.Lift
        V_LO  = self.GetV_LO()
        Vmax  = self.VmaxPlt
        
        #TODO: For debugging
        # This is the stall speed
        #V_LO = V_LO/self.Wing.V_LOstall + 0.1*FT/SEC
        
        A_GR = self.param.A_GR / (ARCDEG)
        
        Vs = self.VelocityRange(0*FT/SEC, Vmax, 35)
                
        amin = self.Wing.AlphaClmin() / ARCDEG
        amax = self.Wing.AlphaClmax() / ARCDEG
        
        AlphaTrim = []
        V_Interp = []
        V_LO_i = 0

        for i in range(1,len(Vs)):
            if Vs[i] >= V_LO and Vs[i-1] < V_LO:
                Vs = npy.insert(Vs/(FT/SEC), i, [V_LO/(FT/SEC)] )*FT/SEC
                V_LO_i = i
                break
        
        for i in range(len(Vs)):
            if i < V_LO_i:
                AlphaTrim.append( A_GR )
                V_Interp.append( Vs[i] / (FT/SEC) )
                continue
                                    
            def y(Alpha2d):
                return (Lift(Vs[i], Alpha2d*ARCDEG)- self.TotalWeight) / (LBF)
            
            try:
                alpha = brentq(f = y, a = amin, b = amax, rtol = 1.e-2, maxiter = 1000)
                
                if i == V_LO_i:
                    #Add an additional point to make the interpolation here vertical
                    AlphaTrim.append( A_GR )
                    V_Interp.append( Vs[i] / (FT/SEC) )
                
                AlphaTrim.append( alpha )
                V_Interp.append( Vs[i] / (FT/SEC) )

                
            except ValueError:
                print "WARNING :: " + self.name + ": Failed to compute trimmed angle of attack for velocity " + AsUnit(Vs[i], 'ft/s')
                AlphaTrim.append( AlphaTrim[i-1] )
                V_Interp.append( Vs[i] / (FT/SEC) )
       
        self.param.AlphaTrim = interp1d(V_Interp, AlphaTrim, kind='linear')
    
#===============================================================================
    def AlphaTrim(self, V):
        """
        Calculates the Alpha Trim for the aircraft
        
        Inputs:
            V - Velocity
        """
        if self.dirty: self.Refresh()
        #
        # Return the interpolate 2D trimmed angle of attack if the aicraft has taken off
        # otherwise return the 2D ground roll angle of attack
        # 
        
        V_LO = self.GetV_LO()
        A_GR = self.param.A_GR / (ARCDEG)
        atrim = []
        
        V = self.MakeList(V)
        
        for v in V:
#            if v > V_LO:
            atrim.append(self.param.AlphaTrim(v / (FT/SEC)))
#            else:
#                atrim.append(A_GR)
        return atrim[0]*ARCDEG if len(atrim) == 1 else npy.array(atrim)*ARCDEG
    
#===============================================================================
    def V_flt(self, Alpha2d):
        """
        Determines flight velocity required for a given angle of attack
        
        Inputs:
            Alpha2d -  Two dimensional angle of attack for the wing
        """
        CLTrim  = self.CLTrim(Alpha2d)
        S       = self.GetS()
        h       = self.Wing.GetAlt_LO()
        Weight  = self.TotalWeight
        
        V_flt = (2 * Weight / (CLTrim * S * AeroUtil.DensAlt(h)))**0.5
        
        return V_flt

#===============================================================================
    def Vmax(self):
        """
        Returns the maximum flight velocity of the aircraft
        """
        if self.dirty: self.Refresh()
        
        return self.param.Vmax

#===============================================================================
    def V_max_climb(self):
        """
        Returns the velocity at the maximum climb of the aircraft
        """
        if self.dirty: self.Refresh()
        
        return self.param.V_max_climb
    
#===============================================================================
    def _UpdateCLTrim_max(self):
        """
        Calculates the lift of velocity i.e. when TotalWeight == Lift
        """

        amin   = self.Wing.AlphaClmin() / ARCDEG
        alphas = self.Wing.AlphaRange(nalpha=2) / ARCDEG

        da = 0.0001*ARCDEG
        def y(Alpha2d):
            return (( self.CLTrim(Alpha2d*ARCDEG + da) - self.CLTrim(Alpha2d*ARCDEG) )/da ) /(1./ARCDEG)
        
        try:
            alpha_CLmax = brentq(f = y, a = 0, b = alphas[1]-2*da/ARCDEG, rtol = 1.e-2, maxiter = 1000)*ARCDEG
        except ValueError:
            print "WARNING: Failed to compute trimmed CLmax"
            alpha_CLmax = self.Wing.AlphaClmax()

        self.param.CLTrim_max = self.CLTrim(alpha_CLmax)
        
#===============================================================================
    def _UpdateV_LO(self):
        """
        Calculates the lift of velocity i.e. V_LO = V_LOstall * VStall
        """

        #
        # Compute the stall speed of the aircraft
        #
        Alt = self.GetAlt_LO()
        rho = AeroUtil.DensAlt(Alt)
        Sw  = self.Wing.S
        CLmax = self.param.CLTrim_max
                                                   
        VStall = (2.0*self.TotalWeight/(rho*CLmax*Sw))**0.5

        #
        # Compute the lift of velocity of the aircraft
        #
        self.param.V_LO = self.Wing.V_LOstall * VStall

#===============================================================================
    def _UpdateVmax(self):
        """
        Calculates the maximum flight velocity of the aircraft
        """
        
        def y(Vmax):
            return (self.Thrust(Vmax*FT/SEC)- self.Drag(Vmax*FT/SEC)) / (LBF)
        
        V_LO = self.GetV_LO() / (FT/SEC)
        Vmax = self.VmaxPlt / (FT/SEC)
        #
        # Set the initial guesses and solve for thrust equal drag
        #
        x0 = V_LO
        x1 = x0 + 1
#        self.param.Vmax = SecantSolver(x0 = x0, x1 = x1, y = y, tol = 0.1, itmax = 500, xmin = V_LO, xmax = Vmax)*FT/SEC
        try:
            self.param.Vmax = brentq(f = y, a = V_LO / 1.2, b = Vmax, rtol = 0.1, maxiter = 500) * FT/SEC
        except ValueError:
            print "ERROR: Failed to compute the maximum velocity"
            V = npy.linspace(V_LO / 1.2, Vmax, 20)*FT/SEC
            T = self.Thrust(V) / LBF
            D = self.Drag(V) / LBF
            Vplt = V / (FT/SEC)
            pyl.figure(1)
            pyl.plot(Vplt,T)
            pyl.plot(Vplt,D)
            lgnd = ['Thrust', 'Drag']
            pyl.title('Thrust/Drag vs. Velocity')
            pyl.xlabel('Velocity (ft/s)') ; pyl.ylabel('Thrust/Drag (lbf)')
            pyl.legend(lgnd, loc = 'best')
            pyl.show()
            raise
        
#===============================================================================
    def _UpdateV_max_climb(self):
        """
        Calculates the flight velocity at the maximum rate of climb of the aircraft
        """
        Vmin = self.GetV_LO()*1.01
        Vmax = self.VmaxPlt
        
        V = self.VelocityRange(Vmin, Vmax, 30)
        
        RoC = self.Rate_of_Climb(V)
        
        V_max_climb = V[0]
        RocMax = RoC[0]
        for i in range(len(V)):
            if RoC[i] > RocMax:
                RocMax = RoC[i]
                V_max_climb = V[i]
            
        self.param.V_max_climb = V_max_climb

#===============================================================================
    def _da(self, Cfn, alph2dw, alph3dw):
        """
        Calculates the slope of the function Cfn
        
        Inputs:
            Cfn      - The function from which a slope is desired
            alpha2dw - The 2D wing angle of attack
            alpha3dw - The 3D wing angle of attack
        """
        if self.dirty:  self.Refresh()
                               
        #
        # Calculate the Cfn slope
        #
        dCfn_da = (Cfn(alph2dw[1]) - Cfn(alph2dw[0]))/(alph3dw[1] - alph3dw[0])

        return dCfn_da
             
#===============================================================================
    def Groundroll(self,nV = None):
        """
        Calculates ground roll for the aircraft at the given velocity.
        
        This uses Method 3 - From Raymer, pg. 565, eq. 17.101
        
        Inputs:
            nV - Used for creating a range of velocities
                 Returns GR, V
                 V[-3] = GR without rotation
                 V[-2] = GR with rotation
                 V[-1] = GR with rotation
              
        """
        if self.dirty: self.Refresh()
                            
        h       = self.GetAlt_LO()
        rho     = AeroUtil.DensAlt(h)
        S       = self.GetS()
        Weight  = self.TotalWeight
        Mu_r    = self.Mu_r
        A_GR    = self.param.A_GR
        
        #
        # Compute the lift and drag coefficients for the 2D ground roll angle of attack
        #
        CLTrim  = self.CLTrim(A_GR)
        CDTrim  = self.CDTrim(A_GR)
        
        if Mu_r < 0.00001:
            KA = 1.e-10*SEC**2/M**2  #Needed for Mu_r == 0 with catapults
        else: 
            KA = rho/(2*(Weight/S))*(Mu_r*CLTrim - CDTrim)
        
        #
        # This method (Method 3) takes into account
        # the decrease in thrust with increasing
        # velocity. F_avg is the average thrust per
        # velocity range. By feeding in the effective
        # thrust decrease, an accurate distance may
        # be calculated.
        #
        T = self.Thrust

        def SG(Vi, Vf):
            F_avg = (T(Vi) + T(Vf))/2
            KT    = (F_avg/Weight)  - Mu_r
            KAi = 1/(2*g*KA)
            Vr  = ((KT + KA*Vf**2)/(KT + KA*Vi**2))
            return KAi * math.log(Vr)


        def IntegrateGroundRoll(v):
            dV = 1*FT/SEC
            Vf = 0*FT/SEC
            Vold = 0*FT/SEC
            GR = 0*FT
            #
            # Integrate over unit steps of velocity
            #
            while(Vf + dV < v):
                Vf += dV
                GR += SG(Vold, Vf)
                Vold = Vf
            #
            # Add the last fraction of velocity
            #
            GR += SG(Vold, v)
            
            return GR

        #
        # Integrate (i.e. sum) SG now to get the total ground roll distance
        # 
        # When multiple velocities are desired, the ground roll is saved of
        # for each velocity as the integration is performed up to the maximum
        # velocity
        #
        
        GR = []
        
        
        V_LO = self.GetV_LO()
        
        if nV is None:
            #
            # Simply integrate up to the lift off velocity
            #
            GR = IntegrateGroundRoll(V_LO)
            #
            # Add rotation time
            #
            GR += V_LO*self.RotationTime
            return GR
        
        
        #
        # Sectond way of getting ground roll
        #
        
        V = npy.linspace(0, V_LO/(FT/SEC), nV)
        V = [v for v in V]
        
        #
        # Create the list of GR up to V_LO
        #
        for v in V:            
            GR.append( IntegrateGroundRoll(v*FT/SEC) /FT )
            
        
        GR.append(GR[-1] + (V_LO*self.RotationTime) / (FT) )
        V.append(V_LO/(FT/SEC)+0.00001)   

        #
        # Stop increasing the groundroll as the aircraft is now in the air
        #
        GR.append(GR[-1])
        V.append(1.1*V_LO/(FT/SEC)) 
                

        return npy.array(GR)*FT, npy.array(V)*FT/SEC

#===============================================================================
    def Rate_of_Climb(self, V = None):
        """
        Calculates the Rate of Climb for an Aircraft
        
        Inputs:
            V - Velocity of the aircraft (default V_LO)
        """
        if self.dirty: self.Refresh()

        if V is None:
            V = self.GetV_LO()+0.1*FT/SEC #For floating point comparison
        
        T    = self.Thrust
        Drag = self.Drag
        V_LO = self.GetV_LO()
        W    = self.TotalWeight
        
        V = self.MakeList(V)
        
        ROC = []
        for v in V:
            if v > V_LO:
                ROC.append((v*((T(v)-Drag(v))/W)) / (FT/MIN))
            else:
                ROC.append(0)
                
        return ROC[0]*FT/MIN if len(ROC) == 1 else npy.array(ROC)*FT/MIN

#===============================================================================
    def Climb_Angle(self, V):
        """
        Calculates the climb angle for an Aircraft
        
        Inputs:
            V - Velocity of the aircraft
        """
        T    = self.Thrust
        Drag = self.Drag
        V_LO = self.GetV_LO()
        W    = self.TotalWeight
        
        CA = []
        for v in V:
            if v > V_LO:
                CA.append(math.asin(((T(v) - Drag(v))/W) ))
            else:
                CA.append(0)
                
        if len(CA) == 1:
            return CA[0]*RAD
        else:
            return npy.array(CA)*RAD

#===============================================================================
    def PlotVNDiagram(self, fig=1, TotalWeights=None, nlimit = None):
        """
        Plots the V-N diagram for the aircraft.
        
        Inputs:
            fig          -- Figure where the V-N diagram is plotted
            TotalWeights -- A list of aircraft total weights
        """
        
        if TotalWeights is None:
            TotalWeights = [self.TotalWeight]
        
        h    = self.GetAlt_LO()
        Vmax = self.Vmax()
        qmax = AeroUtil.q(h,Vmax) 
        S    = self.GetS()
        Vne  = 1.10 * Vmax
        
        # Create ranges
        Alpha  = self.AlphaRange()
        CL     = self.CL(Alpha)
        Vrange = self.VelocityRange(Vmax=Vne)
        
        # Make velocities for plotting
        lVr = len(Vrange)
        Vplt = [0]*2*lVr
        for j in range(lVr):
            k = 2*lVr-j-1
            Vplt[j] = Vrange[j] / (FT/SEC)
            Vplt[k] = Vrange[j] / (FT/SEC)
        
        Clmax = max(CL)
        Clmin = min(CL)
        
        # Maximum G forces calculation method        
        def ncalc(V,TW,ClLimit):
            q = AeroUtil.q(h,V)
            n = (q * S * ClLimit / TW)
            n_lim = (qmax * S * ClLimit / TW) if nlimit is None else nlimit
            try:
                nlist = [0] * len(n)
                for i in range(len(n)):
                    if abs(n[i]) < abs(n_lim):
                        nlist[i] = max(n[i],0)
                    else:
                        nlist[i] = max(n_lim,0)
                return nlist
            except:
                return max(n,0)

        legend = []
        pyl.figure(fig)
        colors = 'bgrycmk'
        
        OrigTotalWeight = self.TotalWeight
        i=0;
        for TWeight in TotalWeights:
            self.TotalWeight = TWeight
            V_LO = self.GetV_LO()
            
            n     = ncalc(Vrange,TWeight,Clmax)  #Get max n
            n_min = ncalc(Vrange,TWeight,Clmin)  #Get min n
            nVLO  = ncalc(V_LO,TWeight,Clmax)
            
            # Append min n onto n array for plotting
            for j in range(len(n_min)):
                n.append(n_min[len(n_min)-j-1])
                         
            pyl.plot(Vplt,n,colors[i])
            pyl.plot([V_LO / (FT/SEC)],[nVLO],colors[i]+'o')
            #legend.append('W = %2.1f lbf' % (TWeight / LBF))
            #legend.append('Lift Off')
            i += 1
            
        pyl.legend(legend,loc='best',numpoints=1)
        pyl.title('Loading vs. Flight Velocity')
        pyl.xlabel("Aircraft Velocity (ft/sec)")
        pyl.ylabel("G loading")
        pyl.grid(b=True)
        
        #rest the total weight
        self.TotalWeight = OrigTotalWeight      
        
#===============================================================================
    def PlotPropulsionPerformance(self, fig = 1, Vmax = None, TotalWeights = None, LgndVar = 'W'):
        """
        Plots the thurst and drag vs. velocity along with rate of climb and climb angle.

        Inputs:
            fig          - The figure to draw to
            Vmax         - Maximum velocity to plot to
            TotalWeights - A list of weights to iterate over
            LgndVar      - The variable to use in the legend
        """
        
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)
        
        if Vmax is None:
            Vmax = self.VmaxPlt
            
        #
        # Create the velocities
        #
        V = self.VelocityRange( 0*FT/SEC, Vmax, int(Vmax / (FT/SEC)) )
        Vplt = V / (FT/SEC)
        
        
        #
        # Initialize the Legend Arrays
        #
        LegendD = []
        LegendROC = []
        LegendCA = []
        
        maxRoC = 0
        
        
        OldTW = self.TotalWeight
        #
        # Compute drag for possibly several total weights
        #
        if TotalWeights is None:
            TotalWeights = [self.TotalWeight]
        
        for TWeight in TotalWeights:
            if self.TotalWeight != TWeight:
                self.TotalWeight = TWeight
            
            if LgndVar == 'AR':
                self.Wing.Lift_LO = TWeight #change the lift at lift off of the wing
                
                lgndstr = 'Drag AR = %4.2f' % self.Wing.AR
            else:
                lgndstr = 'Drag' #'Drag W = %1.1f lbf' % (TWeight / LBF)
            
            Drag = self.Drag(V) / LBF

            AeroDrag = []
            for v in V:
                AeroDrag.append( self.AerodynamicDrag(v, self.AlphaTrim(v)) / LBF )

            pyl.subplot(131)
            pyl.plot(Vplt, Drag, label = lgndstr)
            pyl.plot(Vplt, AeroDrag, label = 'Drag (No Ground Roll)' )
            #LegendD.append(lgndstr)
            
            
            ROC = self.Rate_of_Climb(V)
            maxRoC = max(maxRoC, npy.max(ROC / (FT/MIN)) )
            pyl.subplot(132)
            pyl.plot(Vplt,ROC / (FT/MIN))
            LegendROC.append(lgndstr)
            

            CA = self.Climb_Angle(V)
            pyl.subplot(133)
            pyl.plot(Vplt,CA / (ARCDEG))
            LegendCA.append(lgndstr)
            
            
        #
        # Add the Thrust last so the legends for the 4 plots match.
        #
        T = self.Thrust(V)
        pyl.subplot(131)
        pyl.plot(Vplt, T / (LBF), label='Thrust')
        #LegendD.append('Thrust')
        pyl.grid(b=True)

        pyl.subplot(131)
        pyl.ylim(ymin = 0)
        pyl.xlabel('Velocity (ft/s)'); pyl.ylabel('Thrust/Drag (lbf)')
        #pyl.legend(LegendD, loc = 'best')
        pyl.legend(loc = 'best')
        pyl.grid(b=True)

        pyl.subplot(132)
        pyl.ylim(ymin = 0, ymax = maxRoC*1.1)
        pyl.xlabel('Velocity (ft/s)'); pyl.ylabel('Rate of Climb (ft/min)')
        #pyl.legend(LegendROC, loc = 'best')
        pyl.grid(b=True)
        
        pyl.subplot(133)
        pyl.ylim(ymin = 0)
        pyl.xlabel('Velocity (ft/s)'); pyl.ylabel('Climb Angle (deg)')
        #pyl.legend(LegendCA, loc = 'best')
        pyl.grid(b=True)

        pyl.annotate("Aircraft Propulsion Performance", xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)
        
        # Reset the total weight
        if OldTW != self.TotalWeight:
            self.TotalWeight = OldTW

#===============================================================================
    def PlotLiftVelocityGroundroll(self, fig = 1, GrndRollLimit = 190*FT, RulesRollLimit = 200*FT):
        """
        Draws the velocity and lift vs. groundroll

        Inputs:
            fig  - The figure to draw to
        """
        GR, V = self.Groundroll(nV = 30)
        
        LiftComp, CLLegend = self._LiftComponents(V)
        
        Lift = 0*LBF
        for lc in LiftComp:
            Lift += lc
        CLLegend.append('Total Lift')
        
        #V_LO = V_LO / (FT/SEC)
        V    = V / (FT/SEC)
        Lift = Lift / LBF
        GR   = GR / FT
        
        lines = []
        pyl.figure(fig)
        lines.append( pyl.plot(GR, V, 'm') )
        #pyl.axhline(y = V_LO, color = 'b')
        pyl.ylabel('Velocity (ft/s)')
        pyl.xlabel('Ground Roll (ft)')

        pyl.twinx()
        for i in xrange(len(LiftComp)):
            lines.append( pyl.plot(GR, LiftComp[i] / (LBF) ) )
        
        lines.append( pyl.plot(GR, Lift) )
        pyl.axvline(x = RulesRollLimit / FT, color = 'r')
        pyl.axvline(x = GrndRollLimit / FT, color = 'b')
        pyl.xlim(xmax = RulesRollLimit*1.05 / FT )
        pyl.ylabel('Lift (lbf)')
        pyl.legend(lines, ['Velocity'] + CLLegend, loc = 'best')
        pyl.title('Velocity and Lift vs. Groundroll')
        pyl.grid(b=True)

#===============================================================================
    def PlotWeightPrediction(self, TeamName, TeamNumber, fig = 1, TotalWeight=None, EmptyWeight=None, h=None, PayloadFraction = False, ShowDesign = True):
        """
        Draws the velocity and lift vs. groundroll

        Inputs:
            fig         - The figure to draw to
            TotalWeight - The total weight lifted at a given altitude
            EmptyWeight - The empty weight of the aircraft
            Alt         - The altitude at which the given total weight was lifted 
        """
        
        if TotalWeight is None:
            TotalWeight = self.TotalWeight
            
        if EmptyWeight is None:
            EmptyWeight = self.EmptyWeight
            
        if h is None:
            h = self.Wing.Alt_LO
        
        pyld = (TotalWeight-EmptyWeight) / LBF
        
        if PayloadFraction:
            pyld /= (TotalWeight / LBF)
            
        #
        # Calculate the density altitude for the given condition and create a linear
        #   curve for increased or decreased density altitude
        #
        rho1 = AeroUtil.DensAlt(h) / (LBM/FT**3)
        rho2 = AeroUtil.DensAlt(h + 1000*FT) / (LBM/FT**3)
        
        # slope = change in W / change in h
        # pyld = m*x + b  so   plyd(h1) - m*(h1) = b
        m = pyld * (rho2/rho1 - 1) / 1000
        b = pyld - m * h / FT        
        
        pyl.figure(fig, figsize=(10.85,8.35)) #Account for margin when saving as pdf
        if ShowDesign:
            pyl.plot([h / FT],[pyld],'ko')
            pyl.legend(['Design'],loc='best',numpoints = 1)

        pyl.plot([0,5000],[b,b+m*5000],'k')
        if PayloadFraction:
            pyl.title('Payload Fraction vs. Altitude')
            pyl.ylabel('Payload Fraction')
            string = 'Payload Fraction = '
            title = 'PAYLOAD FRACTION PREDICTION CHART'
        else:
            pyl.title('Payload Weight vs. Altitude')
            pyl.ylabel('Payload (lbf)')
            string = 'Payload = '
            title = 'PAYLOAD PREDICTION CHART'

        pyl.xlabel('Density Altitude (ft)'); 
        pyl.grid(b=True)
        
        string += '%4.2f - %7.6f * Density Altitude' % (b,-1*m)
        pyl.annotate(string, xy=(0.2, 0.15),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=15)

        pyl.annotate(title, xy=(.05, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=14)

        #Team and school name must be up top
        pyl.annotate(TeamName + ', UNIVERSITY OF CINCINNATI', xy=(0.95, .975), 
                xycoords='figure fraction',
                horizontalalignment='right', verticalalignment='top',
                fontsize=14)
        
        #The team number is required to be bottom right
        pyl.annotate('TEAM #%(#)03d'%{'#':TeamNumber}, xy=(0.95, .025), 
                xycoords='figure fraction',
                horizontalalignment='right', verticalalignment='bottom',
                fontsize=14)
         
#===============================================================================
#    def PlotPayloadvsGrossW(self, fig, WeightRange=None):
#        """
#        Plots the change in Payload with changing gross weight
#
#        Inputs:
#            fig         - The figure to draw to
#            WeightRange - A range of lifts for calculating aircraft weights
#        """
#        
#        EmptyW   = []
#        PayloadW = []
#        TotalW   = []
#        WingS    = []
#        
#        if WeightRange is None:
#            WeightRange = npy.linspace(20,40,5)*LBF
#        
#        Wplt = WeightRange / (LBF)
#        
#        for w in WeightRange:
#            self.Wing.Lift_LO = w*1.1
#            self.TotalWeight  = w
#            
#            WingS.append(self.Wing.S / (FT**2))
#            TotalW.append(w / (LBF))
#            EmptyW.append(self.EmptyWeight / (LBF))
#            PayloadW.append(self.PayloadWeight() / (LBF))
#            
#        pyl.figure(fig)
#        pyl.plot(Wplt,TotalW,'b',Wplt,EmptyW,'r',Wplt,PayloadW,'g',Wplt,WingS,'c')
#        pyl.legend(['Total Weight','Empty Weight','Payload Weight','Wing Area'],loc = 'best')
#        pyl.xlabel('Total Weight'); pyl.ylabel('Weight (lbf)')
        
#===============================================================================
#    def PlotRateOfClimbAndClimbAngle(self, fig = 1, Vmax = None):
#        """
#        Draws the rate of climb and climb angle vs. velocity
#
#        Inputs:
#            fig          - The figure to draw to
#            Vmax         - Maximum velocity to plot to
#
#        """
#        
#        if self.dirty: self.Refresh()
#        
#        pyl.figure(fig)
#        
#        if Vmax is None:
#            Vmax = self.VmaxPlt
#            
#        #
#        # Create the velocities
#        #
#        V = self.VelocityRange( 0*FT/SEC, Vmax, int(Vmax / (FT/SEC)) )
#        Vplt = V / (FT/SEC)
#        
#        #
#        # Create the plots.
#        #
#        ROC = self.Rate_of_Climb(V)
#        ClimbAngle = self.Climb_Angle(V) 
# 
#        pyl.subplot(121)
#        pyl.plot(Vplt, ROC / (FT/MIN))
#        pyl.ylim(ymin = 0)
#        pyl.xlabel('Velocity (ft/s)'); pyl.ylabel('Rate of Climb (ft/min)')
#
#        pyl.subplot(122)
#        pyl.plot(Vplt, ClimbAngle / (ARCDEG))
#        pyl.ylim(ymin = 0)
#        pyl.xlabel('Velocity (ft/s)'); pyl.ylabel('Climb Angle (deg)')
#
#        pyl.annotate("Aircraft Thurst and Lift/Drag", xy=(.025, .975),
#                xycoords='figure fraction',
#                horizontalalignment='left', verticalalignment='top',
#                fontsize=20)

        
################################################################################
class ACTailAircraft(ACAircraft):
    """
    A class for working with an aircraft with a tail

    Attributes:
        CLHTSlopeAt      - Fuselage angle of attack used to evaluate slope of the horizontal tail CL
        DWSlopeAt        - Fuselage angle of attack used to evaluate slope of the downwash from the main wing
        RotationAngle    - The angle the aircraft can rotate before the tail hits the ground
        NoseGearOffset   - Offset of the nose gear from the front bulkhead (default 0*in)
        VTailPos         - Horizontal location of the vertical tail in percentage of the horizonalt tail semi-span
        HTailPos         - Vertical location of the horizontal tail in percentage of the vertical tail semi-span
        EngineAlign      - Allignment of the engine (1 - top, 0 - middle, -1 - bottom)
        WingXOmega       - Relaxation factor for positioning the wing. Higher number relaxes the solution process. (default 1)
        WingXMaxIt       - Maximum number of iterations used to position the wing (default 10)
    """
#===============================================================================
    def __init__(self):
        super(ACTailAircraft, self).__init__()
        UnitList = self.UnitList

        #
        # Aircraft attributes
        #
        self.NoneList['CLHTSlopeAt']      = None        ; UnitList['CLHTSlopeAt']    = Angle
        self.NoneList['DWSlopeAt']        = None        ; UnitList['DWSlopeAt']      = Angle
        self.__dict__['RotationAngle']    = 0 * ARCDEG  ; UnitList['RotationAngle']  = Angle
        self.__dict__['NoseGearOffset']   = 0*IN        ; UnitList['NoseGearOffset'] = Length
        self.__dict__['VTailPos']         = 0           ; UnitList['VTailPos']       = Unitless
        self.__dict__['HTailPos']         = 0           ; UnitList['HTailPos']       = Unitless
        self.__dict__['EngineAlign']      = 0           ; UnitList['EngineAlign']    = Unitless
        self.__dict__['WingXOmega']       = 1           ; UnitList['WingXOmega']     = Unitless
        self.__dict__['WingXMaxIt']       = 10          ; UnitList['WingXMaxIt']     = Unitless
        
        #
        # Aircraft components
        #
        self.__dict__['Wing']       = None #ACMainWing(1)
        self.__dict__['HTail']      = ACHorizontalTailSurf(self.Wing, 2)
        self.__dict__['VTail']      = ACTailSurf(self.Wing, 3)
        self.__dict__['NoseGear']   = ACLandingGear()
        self.__dict__['MainGear']   = ACLandingGear()
        self.__dict__['Fuselage']   = None #ACFuselage()
        self.__dict__['Propulsion'] = None
        #
        # Make sure the aircraft gets dirtied if any of the following are dirtied
        #
        #self.Wing.Parents.append( self )
        self.HTail.Parents.append( self )
        self.VTail.Parents.append( self )
        self.NoseGear.Parents.append( self )
        self.MainGear.Parents.append( self )
        
        #
        # Give some approprate names
        #
        #self.Wing.name  = "Main Wing"
        self.HTail.name = "Horizontal Tail"
        self.VTail.name = "Vertical Tail"
        self.NoseGear.name   = "Nose Gear"
        self.MainGear.name   = "Main Gear"
        #self.Fuselage.name   = "Fuselage"
       
        #
        # Main Gear is always symmetric
        #
        self.MainGear.Symmetric = True 
        
        #ail = self.Wing.AddControl('Aileron')
        #ail.SgnDup = -1. #The aileron duplicate surface defects in the opposite direction

        elev = self.HTail.AddControl('Elevator')

        #
        # Elevator typically extends the entire horizontal tail
        #
        elev.Fb = 1.
        elev.Fc = 0.25
        elev.Ft = 0.

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
    def CM(self, alpha2dw, del_e = 0*ARCDEG, Xcg = None, Reht = None, i = None):
        """
        Calculates the total moment coefficient of the aircraft

        CM = CLw*(Xcg-Xac)/macw + CMw + eta*Sht/Sw*(CLht*(Xcg-Xht)/macw + CMht)

        Inputs:
            alhpah2dw - Angle of attack the Wing
            del_e     - Elevator deflection
            Xcg       - Center of gravity
            Reht      - Reynolds number of the horizontal tail
            i         - Incidence angle of the horizontal tail
        """
        if self.dirty:  self.Refresh()

        Wing  = self.Wing
        HTail = self.HTail

        if Xcg is None:
            Xcg = self.Xcg()

        eta  = HTail.eta
        Sht  = HTail.S
        CLht = HTail.CL(alpha2dw, del_e, Reht, i)
        CMht = HTail.CM(alpha2dw, del_e, Reht, i)
        Xht  = HTail.X[0]

        Sw   = Wing.S
        CLw  = Wing.CL(alpha2dw)
        CMw  = Wing.CM(alpha2dw)
        Xac  = Wing.Xac()
        macw = Wing.MAC()

        return (CLw*(Xcg-Xac)/macw + CMw + eta*Sht/Sw*( CMht + CLht*(Xcg-Xht)/macw ) ) 

#===============================================================================
    def dCM_da(self, del_e = 0*ARCDEG, Xcg = None, Reht = None, i = None):
        """
        Calculates the derivative of total moment coefficient w.r.t angle of attack of the aircraft

        Inputs:
            del_e     - Elevator deflection
            Xcg       - Center of gravity of the aircraft
        """
        if self.dirty:  self.Refresh()

        Wing  = self.Wing
        HTail = self.HTail
        
        a2dw = self.param.CmSlopeAt
        a3dw = self.CMSlopeAt

        if Reht is None:
            Reht  = HTail.Re()

        if Xcg is None:
            Xcg = self.Xcg()

        return ((self.CM(a2dw[1], del_e, Xcg = Xcg, Reht = Reht, i = i) - \
                 self.CM(a2dw[0], del_e, Xcg = Xcg, Reht = Reht, i = i))/(a3dw[1] - a3dw[0]))

#===============================================================================
    def Cmadot(self, Xcg = None):
        """
        TODO: Document what Cmadot is.
        
        It is for controls
        """
        if self.dirty:  self.Refresh()
       
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
        Xnp       = self.Xnp()
        Reht      = self.HTail.Re() #Just so it does not have to be recomputed every iteration
        a_Zero_CM = self.param.a_Zero_CM
        del_e     = 0. * ARCDEG     # Must have zero deflection on the elevator
        MAC       = self.Wing.MAC()
        
        #
        # The equation (Must have zero deflection on the elevator)
        #
        def y_WingX(WingX, iht):
            self.Wing.X[0] = WingX
            self.Wing.SetDirty()
            #Use Xnp here as we want dCm/da = 0 when Xcg == Xnp
            return self.dCM_da(del_e = del_e, Xcg = Xnp, Reht = Reht, i = iht) / (1./RAD)

        #
        # Setup the equation
        #
        def y_iht(WingX, iht):
            self.Wing.X[0] = WingX
            self.Wing.SetDirty()
            return self.CM(a_Zero_CM, del_e, Reht = Reht, i = iht)


        #
        # Solve the two equations simultaneously using Newtons method
        #

        dx = 0.01*IN     #Delta x location of the wing
        di = 0.01*ARCDEG #Delta insidence angle of the tail
        tol = 0.001
        maxit = self.WingXMaxIt
        omega = self.WingXOmega
        it = 0
        WingX = self.Wing.X[0]
        iht   = self.HTail.i
        y_WingX_v = y_WingX(WingX, iht)
        y_iht_v   = y_iht(WingX, iht)
        while(it < maxit and (abs(y_WingX_v) > tol or abs(y_iht_v) > tol)):            
            #
            # Newton's method
            #
            # [J][x] = [b] 
            #
            #[ a = dy_WingX/dx,  b = dy_WingX/di ][dWingX] = [-y_WingX]
            #[ c = dy_iht/dx  ,  c = dy_iht/di   ][diht]   = [-y_iht  ]
            #
            # Compute the elements of the Jacobian using numerical approximations
            #
            a = (y_WingX(WingX+dx, iht   ) - y_WingX_v)/dx
            b = (y_WingX(WingX   , iht+di) - y_WingX_v)/di
            c = (y_iht  (WingX+dx, iht   ) - y_iht_v)/dx
            d = (y_iht  (WingX   , iht+di) - y_iht_v)/di
            
            detJ = a*d - c*b
            
            #
            # J^-1 = 1/detJ [ d -b ]
            #               [ -c a ]
            #
            # Solve the systme of equations
            dWingX = ( d*(-y_WingX_v) - b*(-y_iht_v))/detJ
            diht   = (-c*(-y_WingX_v) + a*(-y_iht_v))/detJ
            
            #
            # Update the values
            #
            WingX += dWingX/omega
            iht   += diht/omega
            it += 1

            self.HTail.Refresh()
            
            y_WingX_v = y_WingX(WingX, iht)
            y_iht_v   = y_iht(WingX, iht)
        
        if abs(y_WingX_v) > tol or abs(y_iht_v) > tol:
#            iht = npy.linspace(-10, 10, 20)
#            CM = npy.empty(20)
#            for i in range(20):
#                CM[i] = y_iht(WingX, iht[i])
#
#            pyl.figure(1)
#            pyl.plot(iht,CM)
#            pyl.axhline(y = 0, color = 'r')
#            pyl.title('Aircraft CM vs. Horizontal Tail Incidence Angle')
#            pyl.xlabel('i (deg)'); pyl.ylabel('CM')
#            pyl.show()
            raise ACAircraftError("""Failed to position Wing and set HTail insidence angle.
                                     This is typically a result of a small tail or large static margine,
                                     or the polar slopes are poorly approximated.
                                     
                                     Alternatively, increasing WingXOmega and/or WingXMaxIt can lead to a solution.""" 
                                     + "\nReisudals : " + str(y_WingX_v) + ", " + str(y_iht_v))
        
        #
        # Save the result
        #
        self.Wing.X[0] = WingX
        self.Wing.SetDirty()
        self.HTail.i = iht
        
        self.dCM_da(Xcg=self.Xnp())

#===============================================================================
    def _Get_del_e_trim(self, Alpha2dw, Xcg = None):
        """
        Calculates the elevator deflection to trim the aircraft at a given angle of attack.

        Inputs:
            Alpha2dw - 2D angle of attack of the wing
            Xcg      - Center of gravity for the calculation
        """
        Reht = self.HTail.Re() #Just so it does not have to be recomputed every iteration
        
        if Xcg is None:
            Xcg = self.Xcg()
            
        #
        # Setup the equation
        #
        def y(del_e):
            return self.CM(a, Xcg = Xcg, del_e = del_e*ARCDEG, Reht = Reht)

        #
        # Give the initial guesses and solve the equation
        #
        del_e_trim = []
        for a in Alpha2dw:
            try:
                
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
        self._PositionTail()
        #
        # The tail must be positioned before the wing
        #
        self._PositionWing()
        self._PositionWingX()
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
        HTail    = self.HTail
        VTail    = self.VTail
        A_Rot    = self.RotationAngle / RAD
        
        # Update CG
        self.HTail.SetAircraftXcg(self.Xcg())
        self.VTail.SetAircraftXcg(self.Xcg())
                                    
        if Fuselage.Tail.Align is None:
            if VTail.FullWing:
                TE = VTail.TE(VTail.b/2)
            else:
                TE = max(VTail.MaxTE(),  HTail.MaxTE())

            Z_tail = (TE - MainGear.X[0])*math.tan(A_Rot)
    
            if VTail.FullWing:
                Z_tail += VTail.b/2
                
            Z_tail = max(Z_tail, Fuselage.Sections[-1].FrontBulk.Top() - Fuselage.Sections[-1].BackBulk.Height/2)
            #Z_tail = max(Z_tail, Fuselage.Sections[-1].BackBulk.Top() - Fuselage.Sections[-1].BackBulk.Height/2)
        else:
            Z_tail = Fuselage.Sections[-1].BackBulk.Top() - Fuselage.Sections[-1].BackBulk.Height/2
            
            
        HTail.X[2] = Z_tail
        HTail.SetDirty()

        Fuselage.TailBulk.X = self.HTail.GetX()
        Fuselage.SetDirty()
        
        VTail.X[1] = HTail.Tip()[1]*self.VTailPos
        
        if VTail.X[0] > Fuselage.TailBulk.X[0] or VTail.FullWing:
            VTail.X[2] = self.HTail.GetX()[2]
        else:
            VTail.X[2] = Fuselage.Top(VTail.X[0])
        
        VTail.SetDirty()

        if self.HTailPos > 0:
            #Now that everything is positioned, move the horizontal tail     
            HTail.X[2] = VTail.X[2] + VTail.b*self.HTailPos
            HTail.SetDirty()
                

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
        # This is not a great solution...
        #
        if isinstance(self.Propulsion.Engine, ACMotor):
            Motor = self.Propulsion.Engine
            Fuselage = self.Fuselage
            
            Fuselage.Nose.MotorBattery.Weight = Motor.Battery.Weight
            Fuselage.Nose.MotorBattery.LWH    = Motor.Battery.LWH
            
            Fuselage.Nose.SpeedController.Weight = Motor.SpeedController.Weight
            Fuselage.Nose.SpeedController.LWH    = Motor.SpeedController.LWH

            Fuselage.dirty = True
        
        #
        # Also update the Vmax and altitude for the propulsion system
        #
        if self.Propulsion.Vmax != self.VmaxPlt:
            self.Propulsion.Vmax = self.VmaxPlt

        if self.Propulsion.Alt != self.GetAlt_LO():
            self.Propulsion.Alt = self.GetAlt_LO()

#===============================================================================
    def _PositionPayload(self):
        """
        Positions the payload such that the CG is correct
        """
        if not self.Fuselage.Payload: return
        
        self.Fuselage.Payload.Weight = self.PayloadWeight()
        self.Fuselage.Payload.Refresh()
        deltaCG = 1*IN
        itmax = 20
        it = 0
        while( abs( deltaCG/IN ) > 0.01 and it < itmax):
        
            TotalCG = self.TotalCG()[0]
            DesignCG = self.DesignCG()[0]
            deltaCG = TotalCG - DesignCG
            incr = deltaCG/self.Fuselage.Payload.param.Section.Length
        
            self.Fuselage.Payload.XcgSecFrac -= incr
            self.Fuselage.Payload.Refresh()
            it += 1
                  
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
        TotalWeight = self.TotalWeight
        
        Payload = self.Fuselage.Payload
        
        #Only if there is a payload
        if Payload:
            TotalCG = (Payload.Weight*Payload.CG() + EmptyWeight*EmptyCG)/TotalWeight
        else:
            TotalCG = EmptyCG

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
        self.HTail.SetWing( Wing )
        self.VTail.SetWing( Wing )
        
        #
        # Make sure the aircraft is dirtied if the wing is dirtied
        #
        self.Wing.Parents.append( self )
        
        #
        # Give a better guess for the initial X value of the wing
        #
        if self.Fuselage is not None and self.Propulsion is not None:
            self.Fuselage.XOffset = self.Propulsion.Length()
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
        if self.Propulsion is not None:
            self.Fuselage.XOffset = self.Propulsion.Length()
            
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
        

        if self.Fuselage is not None:
            self.Fuselage.XOffset = self.Propulsion.Length()
            
            if self.Wing is not None:
                self.Wing.X[0] = self.Fuselage.AircraftCG()[0]
                self.Wing.SetDirty()
                
        #
        # This is not a great solution...
        #
        if self.Fuselage is not None and isinstance(self.Propulsion.Engine, ACMotor):
            Motor = self.Propulsion.Engine
            Fuselage = self.Fuselage
            
            Fuselage.Nose.MotorBattery.Weight = Motor.Battery.Weight
            Fuselage.Nose.MotorBattery.LWH    = Motor.Battery.LWH
            
            Fuselage.Nose.SpeedController.Weight = Motor.SpeedController.Weight
            Fuselage.Nose.SpeedController.LWH    = Motor.SpeedController.LWH

            Fuselage.SetDirty()
        
        
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
        super(ACTailAircraft, self)._CheckConsistent()

        self._CheckEquation(['CLHTSlopeAt'], Need = 1)
        self._CheckEquation(['DWSlopeAt'], Need = 1)

        if self.CLHTSlopeAt[1] <= self.CLHTSlopeAt[0]:
            raise ACAircraftError(self.name + ": CLHTSlopeAt[1] must be grater than CLHTSlopeAt[0]")

        if self.DWSlopeAt[1] <= self.DWSlopeAt[0]:
            raise ACAircraftError(self.name + ": DWSlopeAt[1] must be grater than DWSlopeAt[0]")
               
#===============================================================================
    def Refresh(self):
        """
        Updates all the internal quantities.
        """
        super(ACTailAircraft, self).Refresh()
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
        self.param.ClhtSlopeAt = self.Wing.Alpha2D(self.CLHTSlopeAt)
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
        self._UpdateCLTrim_max()
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
        self._PositionPayload()
                
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
        
        dx = (self.HTail.MaxTE() - WheelBase[0]) / (IN)
        
        #
        # Draw the rotation angle
        #
        RotAngle  = self.RotationAngle / (RAD)

        Rotation.xs = npy.empty(2)
        Rotation.zs = npy.empty(2)
                
        Rotation.xs[0] = WheelBase[0] / (IN)
        Rotation.zs[0] = 0
        
        Rotation.xs[1] = self.HTail.MaxTE() / (IN)
        Rotation.zs[1] = math.tan(RotAngle)*dx

        #
        # Draw the tipping angle
        #
        TipAngle  = self.TippingAngle / (RAD)

        Tipping.xs = npy.empty(2)
        Tipping.zs = npy.empty(2)

        Tipping.xs[0] = WheelBase[0] / (IN)
        Tipping.zs[0] = 0
        
        Tipping.xs[1] = self.HTail.MaxTE() / (IN)
        Tipping.zs[1] = math.tan(TipAngle)*dx
        
        #
        # Make sure the drawings are up to date
        #
        self.Wing.UpdateDrawing()
        self.HTail.UpdateDrawing()
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
        HTail = self.HTail

        #
        # Create a list of all lift components
        #
        CL = []
        CLLegend = []

        Wing = self.Wing
        HTail = self.HTail

        Sw = Wing.S
        CLw = Wing.CL(a2dw)
        CL.append(CLw)
        CLLegend.append('Wing')
        
        Sh  = HTail.S
        CLh = HTail.CL(a2dw, del_e = del_e)
        DHW = self.WingDownWash(a2dw) / (RAD)

        if not hasattr(DHW, '__len__'):
            DHW = [DHW]

        cosdw = []
        for dhw in DHW:
            cosdw.append(math.cos(dhw))
        cosdw = npy.array(cosdw)
            
        CL.append((CLh*Sh/Sw * cosdw)  )
        CLLegend.append('Horizontal Tail')
       
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
        HTail = self.HTail
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
        
        Sh  = HTail.S
        CDh = HTail.CD(a2dw, del_e)
        CD.append(CDh*Sh/Sw)
        CDLegend.append('Horizontal Tail')
                
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
        
        try:
            return CD[0] if len(CD) == 1 else CD
        except TypeError:
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
        HTail = self.HTail

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
        
        eta  = HTail.eta
        Sht  = HTail.S
        CLht = HTail.CL(alpha2dw, del_e, Reht, i)
        CMht = HTail.CM(alpha2dw, del_e, Reht, i)
        Xht  = HTail.X[0]

        CM.append( eta*Sht/Sw*(CMht + CLht*(Xcg-Xht)/macw) )
        CMLegend.append('Horizontal Tail')
                   
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
    def CDTrim(self, a2dw = None, Xcg = None, V = None):
        """
        Calculates the total trimmed  CD of the aircraft with the elevator deflected
        to produce a zero moment on the aircraft

        Inputs"
            a2dw - 2D angle of attack of the main wing
            Xcg  - A center of gravity position
            V    - Velocity
        """
        if self.dirty: self.Refresh()
        
        if V is not None and a2dw is None:
            a2dw = self.AlphaTrim(V)
        
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
        
        return self.HTail.DWF*self.Wing.DownWash(Alpha2dw)

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
        Returns the 2D angle of attack of the aircraft at take of conditions.
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
        
        a_Zero_CM = self.param.a_Zero_CM
        return self.HTail.CL(a_Zero_CM, del_e = self.del_e_trim(a_Zero_CM))

 #===============================================================================
    def GetHT_Design_CD(self):
        """
        Returns the drag coefficient of the horizonatl tail at zero CM. CD is relative to the tail area.
        """
        if self.dirty: self.Refresh()
        
        a_Zero_CM = self.param.a_Zero_CM
        return self.HTail.CD(a_Zero_CM, del_e = self.del_e_trim(a_Zero_CM))
    
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
        ControlNames += self._GetControlNames(self.HTail)
        ControlNames += self._GetControlNames(self.VTail)

        return ControlNames    



#===============================================================================
    def HTailLoad(self, V):
        """
        Returns the downward load on the horizontal tail at lift off

        Inputs:
            fig - The figure to draw to
        """
        
        if self.dirty: self.Refresh()
        
        a2dw = self.AlphaTrim(V)
        del_e = self.del_e_trim(a2dw)
        
        #
        # Calculate lift for the varying rudder deflections
        #
        CL = self.HTail.CL(a2dw, del_e)
        h  = self.GetAlt_LO()
        q  = AeroUtil.q(h,V)
        S  = self.HTail.S 
        L  = CL*q*S
        
        return L

#===============================================================================
    def PlotCMPolars(self, fig = 1, del_e = [0*ARCDEG], XcgOffsets = [0.05, 0.05]):
        """
        Plots CM with the CG at the neutral point, the nominal point, and some offset points

        Inputs:
            fig        - Figure number to draw the plots on
            del_e      - A list of elevator deflection
            XcgOffsets - A list of CG offsets in % MAC to display the range of stability
        """
        alphas = self.Wing.AlphaRange()

        Xnp = self.Xnp()
        Xcg = self.Xcg()

        CMnp  = self.CM(alphas, del_e = 0 * ARCDEG, Xcg = Xnp)
        CMnom = self.CM(alphas, del_e = 0 * ARCDEG, Xcg = Xcg)
        del_e_trim = self.del_e_trim(alphas)
        MAC = self.Wing.MAC()

        CMs = []
        del_e_trims = []
        legend = ['SM ' + str(self.StaticMargin*100) + '% MAC' ]
        for offset in XcgOffsets:
            CMs.append(self.CM(alphas, del_e = 0 * ARCDEG, Xcg = Xcg+offset*MAC))
            del_e_trims.append( self.del_e_trim(alphas, Xcg = Xcg+offset*MAC ) )
            sgn = '+' if offset > 0 else ' '
            #legend.append('$X_{cg}$ ' + sgn + str(offset) + '% MAC')
            legend.append(sgn + str(offset*100) + '% MAC')
        legend.append( '$X_{np}$' )

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

        pyl.plot(alphafus,CMnom)
        for i in range(len(XcgOffsets)):
            pyl.plot(alphafus,CMs[i])
        pyl.plot(alphafus,CMnp)

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

        pyl.axhline(color = 'k')

        pyl.title(r'$C_M$ Polar with Elevator Deflection')
        pyl.legend(legend, loc = 'best', frameon = False)
        pyl.xlabel(r'Aircraft Angle of Attack $\alpha^o$'); pyl.ylabel('$C_M$'); pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
 
 
        #
        # Plot the trimmed elevator deflections
        #
        
        del_e_trims = []
        legend = ['SM ' + str(self.StaticMargin*100) + '% MAC' ]
        for offset in XcgOffsets:
            del_e_trims.append( self.del_e_trim(alphas, Xcg = Xcg+offset*MAC ) )
            sgn = '+' if offset > 0 else ' '
            legend.append(sgn + str(offset*100) + '% MAC')
        legend.append( '$X_{np}$' )
        
        pyl.subplot(133)
        pyl.plot(alphafus, del_e_trim / ARCDEG)
        for i in range(len(XcgOffsets)):
            pyl.plot(alphafus,del_e_trims[i] / ARCDEG)
        pyl.plot( alphafus, self.del_e_trim(alphas, Xcg = Xnp) / ARCDEG )
        
        pyl.title('Trimmed Elevator Deflection')
        pyl.legend(legend, loc = 'best', frameon = False)
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
        CLht = self.HTail.CL(alpha)
        
        #
        # Moments produced by the wing
        #
        Xcg  = self.Xcg()
        Sw   = self.Wing.S
        Xac  = self.Wing.Xac()
        macw = self.Wing.MAC()
    
        CLw  = self.Wing.CL(alpha)
        CMw = self.Wing.CM(alpha)
        
        CMw_Cl = CLw*(Xcg-Xac)/macw
        CMw_total = CMw_Cl + CMw
        
        #
        # Moment for the tail
        #
        eta  = self.HTail.eta
        Sht  = self.HTail.S
        CMht = self.HTail.CM(alpha)
        Xht  = self.HTail.X[0]
        
        CLht = self.HTail.CL(alpha)
        CMht_Cl = eta*Sht/Sw*( CLht*(Xcg-Xht)/macw )
        CMht *= eta*Sht/Sw
        CMht_total = CMht + CMht_Cl
        
        alpha2dh = self.HTail._GetAlpha2DH(alpha, del_e=0, Re=self.HTail.Re(), iht=self.HTail.i) #The 2D angle of attack that the tail experiances
        
        dCL_da   = self.dCL_da() / (1/ARCDEG)
        dCM_da   = self.dCM_da() / (1/ARCDEG)
        dDW_da   = self.dDW_da() #This is already a number
        dCLht_da = self.dCLht_da() / (1/ARCDEG)
        
        #
        # Get the 2D wing angles of attack
        #        
        ClSlopeAt   = self.param.ClSlopeAt[0]
        CmSlopeAt   = self.param.CmSlopeAt[0]
        dwSlopeAt   = self.param.dwSlopeAt[0]
        ClhtSlopeAt = self.param.ClhtSlopeAt[0]

        CLSlopeAt   = self.CLSlopeAt[0]
        CMSlopeAt   = self.CMSlopeAt[0]
        CLHTSlopeAt = self.CLHTSlopeAt[0]
        DWSlopeAt   = self.DWSlopeAt[0]
        
        #
        # Make the slope go through the specified CLSlopeAt.
        #
        CLOffset   = self.CLTrim(ClSlopeAt)                  - dCL_da  *CLSlopeAt / (ARCDEG)
        CMOffset   = self.CM(CmSlopeAt)                      - dCM_da  *CMSlopeAt / (ARCDEG)
        CLhtOffset = self.HTail.CL(ClhtSlopeAt)              - dCLht_da*CLHTSlopeAt / (ARCDEG)
        DWOffset   = self.WingDownWash(dwSlopeAt) / (ARCDEG) - dDW_da*DWSlopeAt / (ARCDEG)
        
        alpha = self.Wing.AlphaFus(alpha) / ARCDEG
        alpha2dh /= ARCDEG
        #
        # These are in degrees without being a _scalar
        #
        alp_da = npy.array([min(alpha), max(alpha)])
        
        CL_da    = alp_da*dCL_da   + CLOffset
        CM_da    = alp_da*dCM_da   + CMOffset
        DW_da    = alp_da*dDW_da   + DWOffset
        CLht_da  = alp_da*dCLht_da + CLhtOffset

        #
        # Plot them
        #
        pyl.subplot(321)
        pyl.plot(alpha, CL)
        pyl.plot(alp_da, CL_da)
        pyl.xlabel(r'Aircraft $\alpha^o$'); pyl.ylabel(r'$C_l$')
        pyl.legend([r'$C_l$',r'$C_l$ Slope'], loc = 'best', frameon = False)
        pyl.grid()

        pyl.subplot(322)
        pyl.plot(alpha, DW)
        pyl.plot(alp_da, DW_da)
        pyl.xlabel(r'Aircraft $\alpha^o$'); pyl.ylabel(r'DownWash')
        pyl.legend([r'Down Wash',r'Down Wash Slope'], loc = 'best', frameon = False)
        pyl.grid()

        pyl.subplot(323)
        pyl.plot(alpha, CM)
        pyl.plot(alp_da, CM_da)
        pyl.xlabel(r'Aircraft $\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.legend([r'$C_M$ Aircraft',r'$C_M$ Slope'], loc = 'best', frameon = False)
        pyl.grid()

        pyl.subplot(324)
        pyl.plot(alpha, CLht)
        pyl.plot(alp_da, CLht_da)
        pyl.xlabel(r'Aircraft $\alpha^o$'); pyl.ylabel(r'$C_L$ Tail')
        pyl.legend([r'$C_L$ Tail',r'$C_L$ Tail Slope'], loc = 'best', frameon = False)
        pyl.grid()

        pyl.subplot(325)
        pyl.plot(alpha, CMw)
        pyl.plot(alpha, CMw_Cl)
        pyl.plot(alpha, CMw_total)
        pyl.plot(alpha, CM)
        pyl.xlabel(r'Aircraft $\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.legend([r'$C_M$ Wing', r'$C_M$($C_L$) Wing', r'$C_M$ Wing Total', r'$C_M$ Aircraft'], loc = 'best', frameon = False)
        pyl.grid()

        pyl.subplot(326)
        pyl.plot(alpha, CMht)
        pyl.plot(alpha, CMht_Cl)
        pyl.plot(alpha, CMht_total)
        pyl.plot(alpha, CM)
        pyl.xlabel(r'Aircraft $\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.legend([r'$C_M$ Tail', r'$C_M$($C_L$) Tail', r'$C_M$ Tail Total', r'$C_M$ Aircraft'], loc = 'best', frameon = False)
        pyl.grid()
        
#        pyl.subplot(326)
#        pyl.plot(alpha, alpha2dh)
#        pyl.xlabel(r'Aircraft $\alpha^o$'); pyl.ylabel(r'2D Tail AoA ($^o$)')
#        pyl.grid()


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
        
        pyl.figure(fig)
        
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
        HTail = self.HTail
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
            return eta*Sht/Sw*(CMht + CLht*(Xcg-Xht)/macw)
        
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
        pyl.legend(legend, loc = 'best', labelspacing = 0.0, frameon = False)

        pyl.subplot(132)
        for i in range(len(CMs)):
            pyl.plot(alphafus, CMs[i])

        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.legend(legend, loc = 'best', labelspacing = 0.0, frameon = False)

        pyl.subplot(133)
        for i in range(len(del_es)):
#            if del_es[i] == 0*ARCDEG
            pyl.plot(alphafus, CMs[i])

        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_M$')
        pyl.legend(legend, loc = 'best', labelspacing = 0.0, frameon = False)

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
        S = self.HTail.S
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
 
        pyl.annotate("Horizontal Tail Characteristics", xy=(.025, .975),
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
   
################################################################################
if __name__ == '__main__':
    from TestFuselage import Fuselage
    from TestPropulsion import Propulsion

    aircraft = ACTailAircraft()
    
    aircraft.SetFuselage(Fuselage)
    aircraft.SetPropulsion(Propulsion)
    
    aircraft.TippingAngle  = 10*ARCDEG
    aircraft.RotationAngle = 10*ARCDEG
    
    aircraft.Alpha_Groundroll = 0*ARCDEG

    aircraft.TotalWeight    = 30*LBF
    aircraft.CMSlopeAt      = 0 * ARCDEG 
    aircraft.CLSlopeAt      = 6 * ARCDEG
    aircraft.CLHTSlopeAt    = 8 * ARCDEG 
    aircraft.DWSlopeAt      = 7 * ARCDEG
    
    aircraft.StaticMargin   = 0.1
    aircraft.Alpha_Zero_CM  = 3 *ARCDEG
    aircraft.Xnp_Init       = 40.5 *IN

    aircraft.Wing.X         = [1. * M, 0 * IN, 10. * IN]

#    aircraft.Wing.Lift_LO       = 26.799 * LBF
    aircraft.Wing.V_Stall       = 39.2 * FT/SEC
    aircraft.Wing.Alt_LO        = 600 * FT
    aircraft.Wing.ClSlopeAt     = 6 * ARCDEG

    aircraft.Wing.S             = 900 * IN**2

#    aircraft.Wing.AR            = 16
    aircraft.Wing.b             = 10. * FT
    aircraft.Wing.TR            = [1, 0.7, 0.7]
    aircraft.Wing.Fb            = [0.3, 0.8, 1]
#    aircraft.Wing.Gam           = [0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG]
    aircraft.Wing.CEdge         = 'LE'
    aircraft.Wing.ConstUpper    = True
    aircraft.Wing.Airfoil       = ('S1223.dat','S1223')
    aircraft.Wing.FullWing      = True
    aircraft.Wing.LLFile        = 'LLTPolar.txt'

    aircraft.Wing.Controls.Aileron.Fc = 0.25
    aircraft.Wing.Controls.Aileron.Fb = 0.42
    aircraft.Wing.Controls.Aileron.Ft = 0.2

    aircraft.HTail.Airfoil   = ('NACA0012.dat','NACA0012')
    aircraft.HTail.AR        = 3
    aircraft.HTail.TR        = 0.5
    aircraft.HTail.S         = 100 * IN**2
#    aircraft.HTail.L         = 30 * IN
    aircraft.HTail.VC        = 0.75
    aircraft.HTail.DWF       = 2.0 #Main wing Down wash factor
    aircraft.HTail.FullWing  = True
    aircraft.HTail.ClSlopeAt = None

    aircraft.HTail.Controls.Elevator.Fc = 0.3

    #Set the sweep about the elevator hinge
    aircraft.HTail.SweepFc  = 1.0 - aircraft.HTail.Controls.Elevator.Fc


    aircraft.VTail.Airfoil      = ('NACA0012.dat','NACA0012')
    aircraft.VTail.VC           = 0.06
    aircraft.VTail.AR           = 1.6
    aircraft.VTail.TR           = 0.7
    aircraft.VTail.Axis         = (0, 1)
#    aircraft.VTail.L            = 51.572 * IN
    aircraft.VTail.S            = 69 * IN**2

    aircraft.VTail.Controls.Rudder.Fc = 0.4

    #Set the sweep about the rudder hinge
    aircraft.VTail.SweepFc      = 1.0 - aircraft.VTail.Controls.Rudder.Fc

    # Aircraft Landing Gear
    MainGear = aircraft.MainGear
    MainGear.Theta = 45*ARCDEG
    MainGear.GearHeight, MainGear.StrutW, MainGear.StrutH = 6 * IN, 0.2 * IN, 0.1 * IN
    MainGear.X = [5 * IN, 2 *IN, MainGear.GearHeight]
    
    NoseGear = aircraft.NoseGear
    NoseGear.Theta = 0*ARCDEG
    NoseGear.GearHeight, NoseGear.StrutW, NoseGear.StrutH = 5 * IN, 0.2 * IN, 0.1 * IN
    NoseGear.X = [0 * IN, 0 * IN, MainGear.GearHeight]

#    Xcg = aircraft.Wing.Xac() + 2. * IN
#    del_e = 0 * ARCDEG
#    alpha2dw = 5 * ARCDEG
#    print 'CM', aircraft.CM(alpha2dw, del_e)
#    print 'dCM_da', aircraft.dCM_da(del_e)
#    print 'Wing AC', aircraft.wing.Xac()
#    print 'Neutral point', aircraft.Xnp()
#    print 'Center of Gravity', aircraft.Xcg()
#    print 'HTail i', aircraft.HTail.i

    V_LO = aircraft.GetV_LO()

    print 'Aircraft Xnp   :', AsUnit( aircraft.Xnp(), 'in' )
    print 'Aircraft Xcg   :', AsUnit( aircraft.Xcg(), 'in' )
    print 'Aircraft CM    :', aircraft.CM(15*ARCDEG, del_e = 10*ARCDEG)
    print 'Aircraft dCM_da:', aircraft.dCM_da(0*ARCDEG, AsUnit(aircraft.Xnp()), '1/rad') 
    print
    print 'Wing Area      :', AsUnit( aircraft.Wing.S, 'in**2' )
    print 'Wing MAC       :', AsUnit( aircraft.Wing.MAC(), 'in' )
    print 'Wing dCl_da    :', AsUnit( aircraft.Wing.dCL_da(), '1/rad' )
    print 'Wing Xac       :', AsUnit( aircraft.Wing.Xac(), 'in' )
    print 'Wing X         :', AsUnit( aircraft.Wing.X[0], 'in' )
    print
    print 'Horiz Area     :', AsUnit( aircraft.HTail.S, 'in**2' )
    print 'Horiz Length   :', AsUnit( aircraft.HTail.L, 'in' )
    print 'Horiz iht      :', AsUnit( aircraft.HTail.i, 'deg' )
    print 'Horiz dCl_da   :', AsUnit( aircraft.HTail.dCL_da(), '1/rad' )
    print 'Horiz Xac      :', AsUnit( aircraft.HTail.Xac(), 'in' )
    print 'Horiz X        :', AsUnit( aircraft.HTail.X[0], 'in' )
    print
    print 'Vert Area      :', AsUnit( aircraft.VTail.S, 'in**2' )
    print 'Vert Span      :', AsUnit( aircraft.VTail.b, 'in' )
    print
    print 'Groundroll     :', AsUnit( aircraft.Groundroll(V_LO), 'ft' )
    print 'A_GR           :', aircraft.param.A_GR
    print 'Drag V_LO      :', aircraft.Drag(V_LO)
    print 'Drag V_LO+1    :', aircraft.Drag(V_LO+1*FT/SEC)
    print
    print 'Propulsion X   :', aircraft.Propulsion.Engine.X


#    alphas = aircraft.Wing.AlphaRange()
#    del_e_trim = aircraft.del_e_trim(alphas)

#    pyl.figure(3)
#    pyl.plot(alphas / (ARCDEG),del_e_trim / (ARCDEG))

#    print
#    print 'Aircraft del_e_trim:', aircraft.del_e_trim(5*ARCDEG)

#    aircraft.WriteAVLMainWing('AVLMainWing.avl')

#    aircraft.Wing.LLTPlots(3)
#    aircraft.Wing.Draw3DWingPolars(5)
#    aircraft.HTail.Draw3DWingPolars(4)
#    aircraft.HTail.Draw2DAirfoilPolars(3)
    
    aircraft.PlotVNDiagram(fig=7,TotalWeights=[20*LBF,30*LBF,40*LBF])
    
    pyl.figure(6)
    V = aircraft.VelocityRange(aircraft.Wing.GetV_LO(), 100*FT/SEC)
    Atrim = aircraft.Wing.AlphaFus(aircraft.AlphaTrim(V))
    pyl.plot(V / (FT/SEC),Atrim / (ARCDEG))

    aircraft.PlotPropulsionPerformance(5)
    aircraft.PlotDragBuildup(4)
    aircraft.PlotTrimmedPolars(3)
    aircraft.PlotCMPolars(2, (-10*ARCDEG, -5*ARCDEG, 0*ARCDEG, +5*ARCDEG, +10 * ARCDEG), (+0.5 * IN, -0.5 * IN))
    aircraft.Draw()
    pyl.show()