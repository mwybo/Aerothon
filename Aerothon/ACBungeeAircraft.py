from __future__ import division # let 5/2 = 2.5 rather than 2
from ACAircraft import ACTailAircraft
from ACBase import LinearForceDensity, Length, Angle, Pressure, Area
from scalar.units import IN, LBF, ARCDEG, gacc, FT, SEC, RAD, M, Pa, N, MIN
from scalar.units import AsUnit
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import cmath as math
import numpy as npy
import pylab as pyl

#
# Integrate the 2nd order ODE for computing the final velocity for a bungee launched airplanes
#
def IntegrateDistanceODE(x, ODE, History = False, tol = 1e-10, dt = 0.05*SEC):
    
    t0 = 0*SEC
    
    int_beta=0.5
    int_lambda=0.25
    
    y0 = [0*FT]
    y0d = [0*FT/SEC]
    y0dd = [0*FT/SEC**2]

    y1 = [0*FT]
    y1d = [0*FT/SEC]
    y1dd = [0*FT/SEC**2]
    
    int_ts = [t0]
    while y1[-1] <= x:
        int_ts.append( int_ts[-1] + dt )
        y1.append(y0[-1]) ; y1_old = y1[-1] + 1*FT
        y1d.append(y0d[-1]) ; y1d_old = y1d[-1] + 1*FT/SEC
        y1dd.append(y0dd[-1]) ; y1dd_old = y1dd[-1] + 1*FT/SEC**2
        while(tol*FT < abs(y1[-1] - y1_old) and tol*FT/SEC < abs(y1d[-1] - y1d_old) and tol*FT/SEC**2 < abs(y1dd[-1] - y1dd_old)):
            
            y1_old = y1[-1]; y1d_old = y1d[-1]; y1dd_old = y1dd[-1]
            
            y1dd[-1] = ODE(x-y1[-1], y1d[-1]) #(k*y1[-1]+self.Thrust(y1d[-1])-self.Drag(y1d[-1])-Fg_comp)/m

            y1[-1] = y0[-1]+y0d[-1]*dt+((0.5-int_beta)*y0dd[-1]+int_beta*y1dd[-1])*dt**2
            y1d[-1] = y0d[-1]+((1-int_lambda)*y0dd[-1]+int_lambda*y1dd[-1])*dt

#        print "Distance traveled         : ", AsUnit( y1[-1], 'ft' )
#        print "Velocity at end of ramp   : ", AsUnit( y1d[-1], 'ft/s' )
        if y1[-1] < y1[-2]:
            print "WARNING: Reversed motion detected. Halting launch velocity calculation."
            y1 = [0*FT]
            y1d = [0*FT/SEC]
            y1dd = [0*FT/SEC**2]
            break
        
        y0.append(y1[-1])
        y0d.append(y1d[-1])
        y0dd.append(y1dd[-1])

    #Interpolate to the last desired value
    if y1[-1] > x:
        y1d[-1]  = (y1d[-1]  - y1d[-2] )/(y1[-1] - y1[-2])*(y1[-1] - x) + y1d[-2]
        y1dd[-1] = (y1dd[-1] - y1dd[-2])/(y1[-1] - y1[-2])*(y1[-1] - x) + y1dd[-2]
        y1[-1] = x

    if History:
        return int_ts, y1, y1d, y1dd
    else:
        return int_ts[-1], y1[-1], y1d[-1], y1dd[-1]
#
# Create an aircraft class to impose the total length restriction
#
class ACBungeeAircraft(ACTailAircraft):
    """
    An aircraft that is launched with a bunge coord
    
    Attributes:
        Bungee_K     - Spring constant on the bungee coord
        Bungee_X     - Length the bungee coord is streched
        Bungee_Alpha - Angle of the launch system
        
    """
#===============================================================================
    def __init__(self):
        
        super(ACBungeeAircraft, self).__init__()
        self.__dict__['Bungee_K']     = 0*LBF/IN ; self.UnitList['Bungee_K']     = LinearForceDensity
        self.__dict__['Bungee_X']     = 0*IN     ; self.UnitList['Bungee_X']     = Length
        self.__dict__['Bungee_Alpha'] = 0*ARCDEG ; self.UnitList['Bungee_Alpha'] = Angle
        self.__dict__['Bungee_E']     = 0*N/M**2 ; 
#        self.Unitlist['Bungee_E']     = Pressure
        self.__dict__['Bungee_Lo']    = 0*IN     ; self.UnitList['Bungee_Lo']    = Length
        self.__dict__['Bungee_Lf']    = 0*IN     ; self.UnitList['Bungee_Lf']    = Length
        self.__dict__['Bungee_OD']    = 0*IN     ; self.UnitList['Bungee_OD']    = Length
        self.__dict__['Bungee_ID']    = 0*IN     ; self.UnitList['Bungee_ID']    = Length
        self.__dict__['Bungee_Area']  = 0*IN**2  ; self.UnitList['Bungee_Area']  = Area
        
        self.Mu_r = 0
        
        del self.__dict__['NoseGear']
        del self.__dict__['MainGear']
                
#===============================================================================
    def _PositionParts(self):
        """
        Positions all the parts on the aircraft
        """
        self._PositionFuselage()
        self._PositionTail()
        #
        # The tail must be positioned before the wing
        #
        self._PositionWing()
        self._PositionWingX()
        self._PositionPropulsion()

#===============================================================================
    def _PositionFuselage(self):
        """
        Positions the nose and main landing gear
        """
        Fuselage   = self.Fuselage
        TAngle     = self.TippingAngle / RAD
        A_GR       = self.Alpha_Groundroll / RAD
        
        Fuselage.XOffset = self.Propulsion.Length()
                
        #
        # The fuselage is generated relative to the foremost bulkhead
        #
        NoseBulk    = Fuselage.GetNoseBulk()
        NoseSection = Fuselage.GetNoseSection()
        
        NoseBulk.X[2] -= Fuselage.XcgSection.FrontBulk.Bottom()
                
        Fuselage.Refresh()
        
#===============================================================================
    def _PositionTail(self):
        """
        Positions the fuselage on the horizontal tail
        """
        Fuselage = self.Fuselage
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

            Z_htail = (TE - Fuselage.Sections[-1].FrontBulk.X[0])*math.tan(A_Rot).real
    
            if VTail.FullWing:
                Z_htail += VTail.b/2
                
            Z_htail = max(Z_htail, Fuselage.Sections[-1].FrontBulk.Top())
            Z_htail = max(Z_htail, Fuselage.Sections[-1].BackBulk.Top() - Fuselage.Sections[-1].BackBulk.Height/2)
        else:
            Z_htail = Fuselage.Sections[-1].BackBulk.Top() - Fuselage.Sections[-1].BackBulk.Height/2
            
            
        HTail.X[2] = Z_htail
        HTail.SetDirty()

        Fuselage.TailBulk.X = self.HTail.GetX()
        Fuselage.SetDirty()
        
        VTail.X[1] = HTail.Tip()[1]*self.VTailPos
        
        if VTail.X[0] > Fuselage.TailBulk.X[0] or VTail.FullWing:
            VTail.X[2] = self.HTail.GetX()[2]
        else:
            VTail.X[2] = Fuselage.Top(VTail.X[0])
        
        VTail.SetDirty()
                       
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
        
        Dfus = Fuselage.CD(a2dw, V, self.GetAlt_LO(), self.Wing)
        CD.append( (Dfus)*self.DragKFactor/Sw )
        CDLegend.append('Fuselage')

        return CD, CDLegend
    
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
        WheelBase = self.Fuselage.Sections[-1].FrontBulk.X
        
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
        Rotation.zs[1] = (math.tan(RotAngle)*dx).real

        #
        # Draw the tipping angle
        #
        TipAngle  = self.TippingAngle / (RAD)

        Tipping.xs = npy.empty(2)
        Tipping.zs = npy.empty(2)

        Tipping.xs[0] = WheelBase[0] / (IN)
        Tipping.zs[0] = 0
        
        Tipping.xs[1] = self.HTail.MaxTE() / (IN)
        Tipping.zs[1] = (math.tan(TipAngle)*dx).real
        
        #
        # Make sure the drawings are up to date
        #
        self.Wing.UpdateDrawing()
        self.HTail.UpdateDrawing()
        self.VTail.UpdateDrawing()
        self.Fuselage.UpdateDrawing()
        self.Propulsion.UpdateDrawing()
            
#===============================================================================
    def OverturnAngle(self):
        """
        Analyze Overturn Angle and throw warning if criteria not met
        Reymer p265 Fig 11.5, p266 paragraph 2
        """ 
        if self.dirty: self.Refresh()
        return 0*ARCDEG
    
#===============================================================================
    def RollingDrag(self, V, Alpha2d = None):
        """
        Calculates the drag of the aircraft due to rolling friction on the ground
        
        Inputs: 
            V       - Velocity
            Alpha2d - Two dimensional angle of attack of the wing
        """
        a2d = self.AlphaTrim(V) if Alpha2d is None else Alpha2d
        #return (self.TotalWeight*math.cos(self.Bungee_Alpha/RAD).real - self.Lift(V, a2d))*self.Mu_r
        return 0*LBF #(self.TotalWeight*math.cos(self.Bungee_Alpha/RAD).real - self.Lift(V, a2d))*self.Mu_r

#===============================================================================
    def BungeeForce(self, x):
        """
        Calculates the force produced by the catapult based on the deformation of the bungee coord
        
        Inputs: 
            x - Current Deflection
        """
        k=self.Bungee_K
        return k*x

#===============================================================================
    def DeadLaunch_V_LO(self):
        """
        Computes the takeoff velocity with the dead engine
        """
        if self.dirty: self.Refresh()
        
        x=self.Bungee_X
        k=self.Bungee_K
        
        Weight = self.TotalWeight
        m = Weight/gacc
        Fg = Weight*math.sin(self.Bungee_Alpha/RAD).real

        A_GR = self.param.A_GR

        def ODE(x, V):
            return (self.BungeeForce(x) - self.AerodynamicDrag(V, Alpha2d = A_GR) - self.RollingDrag(V, Alpha2d = A_GR) - Fg)/m
        
        int_ts, y1, y1d, y1dd = IntegrateDistanceODE(x, ODE)
        
        return y1d
        
#===============================================================================
    def DeadLaunch_Acceleration(self):
        """
        Computes the acceleration for a dead launch
        """
        if self.dirty: self.Refresh()
        
        x=self.Bungee_X
        
        Weight = self.TotalWeight
        m = Weight/gacc
        Fg = Weight*math.sin(self.Bungee_Alpha/RAD).real

        A_GR = self.param.A_GR

        def ODE(x, V):
            return (self.BungeeForce(x) - self.AerodynamicDrag(V, Alpha2d = A_GR) - self.RollingDrag(V, Alpha2d = A_GR) - Fg)/m
        
        int_ts, y1, y1d, y1dd = IntegrateDistanceODE(x, ODE, History=True)
        
        return self.ToUnumList(int_ts,SEC), self.ToUnumList(y1dd, M/SEC**2)
    
#===============================================================================
    def LiveLaunch_V_LO(self, Result = 'v', History = False):
        """
        Computes the takeoff velocity with the live engine
        """
        if self.dirty: self.Refresh()

        angle = self.Bungee_Alpha
        x = self.Bungee_X
        
        Weight = self.TotalWeight
        m = Weight/gacc
        Fg = Weight*math.sin(self.Bungee_Alpha/RAD).real
        
        A_GR = self.param.A_GR

        def ODE(x, V):
            return (self.BungeeForce(x) + self.Thrust(V) - self.AerodynamicDrag(V, Alpha2d = A_GR) - self.RollingDrag(V, Alpha2d = A_GR) - Fg)/m
        
        int_ts, y1, y1d, y1dd = IntegrateDistanceODE(x, ODE, History)
        
        ret = []
        
        if 'x' in Result:
            ret.append( self.ToUnumList(y1, M) )
        if 'v' in Result:
            ret.append( self.ToUnumList(y1d, M/SEC) )
        if 'a' in Result:
            ret.append( self.ToUnumList(y1dd, M/SEC**2) )
        
        return ret[0] if len(ret) == 1 else ret
        
#===============================================================================
    def LiveLaunch_Acceleration(self):
        """
        Computes the takeoff velocity with the live engine
        """
        if self.dirty: self.Refresh()

        angle = self.Bungee_Alpha
        x = self.Bungee_X
        
        Weight = self.TotalWeight
        m = Weight/gacc
        Fg = Weight*math.sin(self.Bungee_Alpha/RAD).real
        
        A_GR = self.param.A_GR

        def ODE(x, V):
            return (self.BungeeForce(x) + self.Thrust(V) - self.AerodynamicDrag(V, Alpha2d = A_GR) - self.RollingDrag(V, Alpha2d = A_GR) - Fg)/m
        
        int_ts, y1, y1d, y1dd = IntegrateDistanceODE(x, ODE)
        
        return int_ts, self.ToUnumList(y1dd, M/SEC**2)
                
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
            a2d = self.AlphaTrim(v) if Alpha2d is None else Alpha2d                    
            Drag.append( self.AerodynamicDrag(v, a2d) / (LBF))
                
        return Drag[0]*LBF  if len(Drag) == 1 else npy.array(Drag)*LBF
        
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

            pyl.subplot(131)
            pyl.plot(Vplt, Drag, label = lgndstr)
            #LegendD.append(lgndstr)
            
            
            ROC = self.Rate_of_Climb(V)
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
#        Xs, Vs = self.LiveLaunch_V_LO(Result = 'xv', History = True)
#        
#        
#        BengeeLength = interp1d(Vs/(FT/SEC), Xs/IN, kind='linear')
#        for i in xrange(len(V)):
#            if V[i] < Vs[-1]:
#                T[i] += self.BungeeForce(self.Bungee_X - BengeeLength(V[i]/(FT/SEC))*IN)
#            else:
#                RealseThrust = T[i]
#                break
            
        pyl.subplot(131)
#        pyl.plot(Vs[-1]/(FT/SEC), RealseThrust/LBF, 'ro', label='Release')
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
        pyl.ylim(ymin = 0)
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