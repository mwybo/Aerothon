"""
This is the class for analyzing the propulsion system of the aircraft
Uses the ACPropeller and ACEngine classes for the anlysis
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACStructural, Length, Velocity, Unitless, Force, Frequency, g
from scalar.units import RPM, FT, SEC, LBF, LBM, HP, HR, TSFC, IN, ARCDEG, PSFC, OZF, A, W, MIN
from scalar.units import AsUnit
from RootSolvers import SecantSolver
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import AeroUtil
import cmath as math
import numpy as npy
import pylab as pyl

################################################################################
class ACPropulsionError(Exception):
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg

################################################################################
class ACPropulsion(ACStructural):
    """
    Propulsion class
    
    Inputs:
        Prop   - Propeller class to be used in the propulsions analysis
        Engine - Engine class to be used in the propulsions analysis  
        
    Attributes:
        Alt    - Altitude of the aircraft
        Vmax   - Maximum velocity of the aircraft
        nV     - Number of velocity interpolation points for the RPMMatch (default = 50)
        Weight - Weight of the total propulsion system
        PRMInit - Initial RPM guess when mathcing power of the engine and propeller (default 10,000)
    """
    
#===============================================================================    
    def __init__(self,Prop,Engine):
        super(ACPropulsion,self).__init__()
        UnitList = self.UnitList
        
        self.name = "Propulsion"
        self.__dict__['Prop']    = Prop
        self.__dict__['Engine']  = Engine
        
        self.__dict__['Alt']     = 0*FT      ; UnitList['Alt']     = Length
        self.__dict__['Vmax']    = 0*FT/SEC  ; UnitList['Vmax']    = Velocity
        self.__dict__['nV']      = 50        ; UnitList['nV']      = Unitless
        self.__dict__['RPMInit'] = 10000*RPM ; UnitList['RPMInit'] = Frequency
        
        #
        # Make sure the propulsion gets refreshed if the prop or engine get dirtied
        #
        self.Prop.Parents.append( self )
        self.Engine.Parents.append( self )

#===============================================================================
    def _CalcWeight(self):
        """
        Calculates the weight of the propulsion as the sum of the prop and engine
        """
        fac = 2 if self.Symmetric else 1
        return fac*(self.Prop.Weight + self.Engine.Weight)

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the propulsion and its parts to the Weight Table
        
        Input:
            PartName    - The name give to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WeightTable[PartName] = self
        
        self.Prop.AddToWeightTable('Prop', WeightTable[PartName])
        self.Engine.AddToWeightTable('Engine', WeightTable[PartName])
        
#===============================================================================
    def Length(self):
        """
        Computes the length of the propulsion system
        """
        return self.Prop.Thickness + self.Engine.Length()
    
#===============================================================================
    def CG(self):
        """
        Computes the center of gravity of the propulsion system
        """
        Prop = self.Prop
        Engine = self.Engine
        return (Prop.CG() * Prop.Weight + Engine.CG() * Engine.Weight) / (Prop.Weight + Engine.Weight)

#===============================================================================
    def MOI(self):
        """
        Computes the moments of inertia of the propulsion system
        """
        Prop = self.Prop
        Engine = self.Engine
        
        CG  = self.CG()
        MOI = npy.array([0,0,0]) * LBM*IN**2

        def AddMOI(PCG, PM, PMOI):
            dCG = PCG - CG
            MOI[0] += PMOI[0] + PM*(dCG[1]**2 + dCG[2]**2) #Ixx
            MOI[1] += PMOI[1] + PM*(dCG[0]**2 + dCG[2]**2) #Iyy
            MOI[2] += PMOI[2] + PM*(dCG[1]**2 + dCG[0]**2) #Izz
        
        AddMOI(Prop.CG()  , Prop.Weight / g  , Prop.MOI())
        AddMOI(Engine.CG(), Engine.Weight / g, Engine.MOI())
        
        return MOI

        
#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACPropulsion,self).Refresh()
        self.param.refreshing = True
       
        self._UpdateRPMMatch()

        self.Prop.Symmetric = self.Symmetric
        self.Engine.Symmetric = self.Symmetric
        
        self.Prop.X = self.Engine.PropX()
        
        self.param.refreshing = False

#===============================================================================
    def UpdateDrawing(self):
        """
        Updates all drawings of the propulsion systemd
        """
        self.Prop.UpdateDrawing()            
        self.Engine.UpdateDrawing()
        
#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """
        Draws this propulsion system

        Inputs:
            fig   - Integer number of figure to draw into
            top   - Subplot for top view
            side  - Subplot for side view
            front - Subplot for front view
        """
        self.UpdateDrawing()
        
        self.Prop.Draw(fig, top, side, front)
        self.Engine.Draw(fig, top, side, front)
       
#===============================================================================    
    def RPMMatch(self,V):
        """
        Matches the RPM of the engine and propeller
            - Outputs the matched RPM
                
        Inputs:
            V      - forward velocity of the aircraft 
        """
        
        if self.dirty:
            self.Refresh()

        V = self.MakeList(V)

        Nbp = []
        for v in V:
            if v > self.Vmax:
                raise ACPropulsionError('ERROR: Maximum velocity exceeded. Increase Vmax for the propulsion.')
            N = self.param.Nbp(v / (FT/SEC))
            Nbp.append(N)
            
        return Nbp[0]*RPM if len(Nbp) == 1 else npy.array(Nbp)*RPM

#===============================================================================    
    def _UpdateRPMMatch(self):
        """
        Updates the interpolating array for the RPMMatch
        """
        
        V  = npy.linspace(0, self.Vmax /(FT/SEC), self.nV)*FT/SEC
        h  = self.Alt
        Pp = self.Prop.P
        Pb = self.Engine.Pb
        dn = 1000
        Nval = self.RPMInit / RPM
        
        # Use the secant solver to find where the propeller power and engine power intersect
        Nbp = []
        for v in V:
            def y(N):
                return (Pp(N*RPM,v,h)- Pb(h,N*RPM)) / HP
#            Nval = SecantSolver(x0 = Nval, x1 = Nval + dn , y = y, tol = 1.e-2, itmax = 100)
            Nval = brentq(f = y, a = 1000, b = 30000, rtol = 1.e-5, maxiter = 1000)
            Nbp.append(Nval)
            
        Nbp = npy.array(Nbp)
        self.param.Nbp = interp1d(V / (FT/SEC), Nbp, kind='linear')
        
#===============================================================================    
    def T(self,V):
        """
        Calculates the tip velocity of the propeller blades
                
        Inputs:
            V      - forward velocity of the aircraft 
        """
        
        Tp = self.Prop.T
        h  = self.Alt
        
        V = self.MakeList(V)
        
        T = []
        for v in V:
            Nbp = self.RPMMatch(v)
            T.append(Tp(Nbp,v,h) / LBF)
        
        fac = 2 if self.Symmetric else 1
        
        return fac*T[0]*LBF if len(T) == 1 else fac*npy.array(T)*LBF 

#===============================================================================    
    def P(self,V):
        """
        Power that the propeller consumes at the matched RPM
                
        Inputs:
            V      - forward velocity of the aircraft 
        """
        
        Pp = self.Prop.P
        h  = self.Alt
        
        P = []
        for v in V:
            Nbp = self.RPMMatch(v)
            P.append(Pp(Nbp,v,h) / HP)
        
        return P[0]*HP if len(P) == 1 else npy.array(P)*HP 
        
#===============================================================================    
    def Vtip(self,V):
        """
        Calculates the tip velocity of the propeller blades
                
        Inputs:
            V      - forward velocity of the aircraft 
        """
        
        d  = self.Prop.D
        pi = math.pi
        
        Nbp = self.RPMMatch(V)
        
        return pi*d*Nbp

#===============================================================================    
    def Efficiency(self,V):
        """
        Calculates the propulsion system useful power and efficiency at the matched RPM
                
        Inputs:
            V      - forward velocity of the aircraft 
        """
        
        Tp = self.Prop.T
        Pp = self.Prop.P
        h  = self.Alt
        
        Pu  = []; Eta = []
        for v in V:
            Nbp = self.RPMMatch(v)
            Pbp = Tp(Nbp,v,h)
            Pbp = Pp(Nbp,v,h)
            Puval = Tbp*v
            
            Pu.append(Puval / HP)
            Eta.append(Puval/Pbp)
        
        if len(Pu) == 1:
            return Eta[0], Pu[0]
        else:
            return npy.array(Eta), npy.array(Pu)*HP

#===============================================================================    
    def Mdotfv(self,V):
        """
        Calculates the fuel mass flow and the TSFC at the matched RPM
                
        Inputs:
            V      - forward velocity of the aircraft 
        """
        
        h     = self.Alt
        Tp    = self.Prop.T
        mdotf = self.Engine.Mdotf
        ps    = self.Engine.PS
         
        mdotfv = []
        tsfc = []
        for v in V:
            Nbp = self.RPMMatch(v)
            #
            # Note that TSFC = LBM/(HR*LBF)
            #
            Tbp = Tp(Nbp,v,h) / LBF
            mdotfvval = mdotf(h,Nbp) / (LBM/HR)
            
            mdotfv.append(mdotfvval)
            tsfc.append(mdotfvval/Tbp)
        
        if len(mdotfv) == 1:
            return mdotfv[0]*LBM/HR, tsfc[0]*TSFC
        else:
            return npy.array(mdotfv)*LBM/HR, npy.array(tsfc)*TSFC
            
#===============================================================================    
    def PlotTPvsN(self, N, V, fig = 1):
        """
        Plots the power of the propeller and engine and the propeller thrust
                
        Inputs:
            N   - A list of RPMs
            V   - A list of velocities
            fig - Figure number
        """
        Engine = self.Engine
        h = self.Alt
        
        pyl.figure(fig)
        
        Pb = self.Engine.Pb(h,N)
        Pb = Pb / Engine.PowerUnit
        
        Nplt  = N / RPM
        
        legend1 = []
        legend2 = []
        
        V = self.MakeList(V)
       
        for v in V:
            Pp = self.Prop.P(N,v,h)
            Tp = self.Prop.T(N,v,h)
            
            Pp = Pp / Engine.PowerUnit
            Tp = Tp / Engine.ThrustUnit
            
            Vplt = AsUnit(v, "ft/s", "%3.0f")
            
            pyl.subplot(121)
            pyl.plot(Nplt,Pp)

            legend1.append('Prop Power Req. ' + Vplt)
            
            pyl.subplot(122)
            pyl.plot(Nplt,Tp)
            
            legend2.append(Vplt)

        pyl.subplot(121)
        pyl.plot(Nplt,Pb, 'b')
        legend1.append('Engine Power Provided')                        
        pyl.xlabel('RPM'); pyl.ylabel('Power (' + Engine.PowerUnitName + ')')
        pyl.legend(legend1,loc ='best')

        pyl.subplot(122)
        pyl.xlabel('RPM'); pyl.ylabel('Propeller Thrust (' + Engine.ThrustUnitName + ')')
        pyl.legend(legend2,loc ='best')
        
        pyl.annotate("Power and Propeller Thrust", xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)


#===============================================================================    
    def PlotElectricMatched(self, V, Vprop, fig = 1, ylimits = [None,None,None,None], PlotProp = True):
        """
        Plots the matched RPM, Power and Thrust for a range of velocities
                
        Inputs:
            V           - A list of velocities at which to match the RPMs 
            N           - A list of PRM values for the engine/propeller power match
            Vprop       - A list of velocities for plotting propeller power usage
            fig         - Figure number
            ylimits     - Y axis ranges for each plot [ (ymin1, ymax1), (ymin2, ymax2), (ymin3, ymax3)]
            thrust_unit - Unit on the thrust
        """
        Motor = self.Engine
        Nm = self.RPMMatch(V)
        
        N  = Motor.NRange()
        Ib = Motor.Ib(Nm)/A

        self.PlotMatched(V, N, Vprop, fig, ylimits, PlotProp)
        
        Battery = Motor.Battery
        SpeedController = Motor.SpeedController
        
        Imax = min(SpeedController.Imax/A, Battery.Imax/A)
        Vmax = self.Vmax / (FT/SEC)
        
        pyl.figure(fig+1)
        
        pyl.subplot(221)
        pyl.title(self.Engine.name + ', ' + self.Prop.name)
        pyl.plot(V/(FT/SEC), Ib)
        pyl.plot([0, Vmax], [Imax, Imax], 'r-.', linewidth=2 )
        ylim = pyl.ylim()
        pyl.ylim( min(ylim[0], Imax*0.95), max( ylim[1], Imax*1.05 ) )
        pyl.xlabel('V (ft/s)'); pyl.ylabel('Current (A)') 
        
        pyl.subplot(222)
        pyl.plot(V/(FT/SEC), Motor.Duration(N=Nm)/MIN)
        pyl.xlabel('V (ft/s)'); pyl.ylabel('Duration (min)')
    
        pyl.subplot(223)
        pyl.plot([0, Vmax], [Motor.Wmax/W, Motor.Wmax/W], 'r-.', linewidth=2 )
        pyl.plot(V/(FT/SEC), Motor.Pin(Nm)/W)
        ylim = pyl.ylim()
        pyl.ylim( min(ylim[0], Motor.Wmax/W*0.95), max( ylim[1], Motor.Wmax/W*1.05 ) )
        pyl.xlabel('V (ft/s)'); pyl.ylabel('Power In (watt)') 
        
        pyl.subplot(224)
        pyl.plot(V/(FT/SEC), Motor.Efficiency(Nm))
        pyl.xlabel('V (ft/s)'); pyl.ylabel('Efficiency (%)') 
        
        
#===============================================================================    
    def PlotMatched(self, V, N, Vprop, fig = 1, ylimits = [None,None,None,None], PlotProp = True):
        """
        Plots the matched RPM, Power and Thrust for a range of velocities
                
        Inputs:
            V           - A list of velocities at which to match the RPMs 
            N           - A list of PRM values for the engine/propeller power match
            Vprop       - A list of velocities for plotting propeller power usage
            fig         - Figure number
            ylimits     - Y axis ranges for each plot [ (ymin1, ymax1), (ymin2, ymax2), (ymin3, ymax3)]
            thrust_unit - Unit on the thrust
        """
        Engine = self.Engine
                
        Nbp = self.RPMMatch(V)
        Tbp = self.T(V)
        Pbp = self.P(V)
                
        V   = V / (FT/SEC)
        Nbp = Nbp / RPM
        Tbp = Tbp / Engine.ThrustUnit
        Pbp = Pbp / Engine.PowerUnit
                
        pyl.figure(fig)
        pyl.subplot(222)
        pyl.plot(V,Nbp)
        pyl.grid(b=True)
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Matched RPM')
        if ylimits[1] is not None:
            pyl.ylim(ylimits[0])
#        else:
#            pyl.ylim( (0, max(Nbp)*1.1) )
        
        pyl.subplot(223)
        pyl.plot(V,Tbp)
        pyl.grid(b=True)
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Thrust at matched RPM (' + Engine.ThrustUnitName + ')')
        if ylimits[2] is not None:
            pyl.ylim(ylimits[1])
#        else:
#            pyl.ylim( (0, max(Tbp)*1.1) )
                    
        pyl.subplot(224)
        pyl.plot(V,Pbp)
        pyl.grid(b=True)
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Power at matched RPM (' + Engine.PowerUnitName + ')')
        if ylimits[3] is not None:
            pyl.ylim(ylimits[2])
#        else:
#            pyl.ylim( (0, max(Pbp)*1.1) )
            
            
        legend = ['Engine']
        
        Vprop = self.MakeList(Vprop)
        h = self.Alt
        Nplt = N / RPM
        
        Pb = self.Engine.Pb(h,N)
        Pb = Pb / Engine.PowerUnit
            
        ax = pyl.subplot(221)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))

        pyl.plot(Nplt,Pb)
        pyl.grid(b=True)
        pyl.xlabel('RPM'); pyl.ylabel('Power (' + Engine.PowerUnitName + ')')

        if ylimits[3] is not None:
            ymax = ylimits[2]
        else:
            ymax = max(max(Pb)*1.5, pyl.ylim()[1])

        if PlotProp:
            for v in Vprop:
                Pp = self.Prop.P(N,v,h)
                
                Pp = Pp / Engine.PowerUnit
                
                Vplt = AsUnit(v, "ft/s", "%3.0f")
                
                pyl.subplot(221)
                pyl.plot(Nplt,Pp)
    
                legend.append('Prop ' + Vplt)
                    
        legend_fontsize = pyl.rcParams['legend.fontsize']
        pyl.rcParams.update({'legend.fontsize': 10})
        
        if PlotProp:
            pyl.legend(legend,loc ='upper left')
                    
        pyl.ylim( ymin = 0, ymax = ymax)
        
#        pyl.annotate("Matched RPM, Thrust and Power", xy=(.025, .975),
#                xycoords='figure fraction',
#                horizontalalignment='left', verticalalignment='top',
#                fontsize=20)
        pyl.annotate(self.Engine.name + ", " + self.Prop.name, xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)
        
        pyl.rcParams.update({'legend.fontsize': legend_fontsize})

#===============================================================================    
    def PlotTipSpeed(self, V, fig = 1):
        """
        Plots tip speed vs. aircraft velocities
                
        Inputs:
            V   - A list of velocities at which to match the RPMs 
            fig - Figure number
        """
        
        Nbp = self.RPMMatch(V)
        Vtip = self.Vtip(V)
        TipMach = AeroUtil.Mach(self.Alt,Vtip)
        
        V   = V / (FT/SEC)
        Vtip = Vtip / (FT/SEC)
       
        pyl.figure(fig)
        pyl.subplot(121)
        pyl.plot(V,Vtip)
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Tip Speed (ft/sec)')
        
        #
        # Plot the tip mach number here
        #
        pyl.subplot(122)
        pyl.plot(V,TipMach)
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Tip Mach Number')
        
                
#===============================================================================    
    def PlotFuelFlow(self, V, fig=1):
        """
        Plots the matched fuel flow and tsfc for a range of velocities
                
        Inputs:
            V   - A list of velocities at which to match the RPMs 
            fig - Figure number
        """
        
        mdotfv, tsfc = Propulsion.Mdotfv(V);
        
        V      = V / (FT/SEC)
        mdotfv = mdotfv / (LBM/HR)
        tsfc   = tsfc / TSFC
        
        pyl.figure(fig)
        pyl.subplot(121)
        pyl.plot(V,mdotfv)
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('Fuel flow at matched RPM (lbm/hr)')
        
        pyl.subplot(122)
        pyl.plot(V,tsfc)
        pyl.xlabel('Aircraft Velocity (ft/sec)'); pyl.ylabel('TSFC at matched RPM (lbm/(hr*lbf))')

#===============================================================================
    def PlotTestData(self, fig = 1):
        """
        Plots measured data for both the propeller and engine models
            
        Inputs:
            fig - Figure number
        """        
        TestData = npy.array(self.Engine.TestData)

        Nmax = max(TestData[:,0]) / RPM    
        Nmin = 0
             
        N  = npy.linspace(Nmin,Nmax,20)*RPM
        Tb = self.Engine.Tb(0*FT, N)
        Pb = self.Engine.Pb(0*FT, N)
        
        N  = N / RPM
        Tb = Tb / self.Prop.TorqueUnit
        Pb = Pb / self.Prop.PowerUnit
                
        Ntest = self.ToUnumList(TestData[:,0],RPM)
        Ttest = self.ToUnumList(TestData[:,1],IN*OZF)
  
        Ptest = (Ntest*Ttest) / self.Prop.PowerUnit
        Ntest = Ntest / RPM
        Ttest = Ttest / self.Prop.TorqueUnit

        subplot4 = None
        if self.Prop.ThrustData is not None and self.Prop.TorqueData is not None and self.Prop.DynThrustData is not None:
            subplot1 = 221; subplot2 = 222;  subplot3 = 223;  subplot4 = 224;
        elif self.Prop.ThrustData is not None and self.Prop.TorqueData is not None and self.Prop.DynThrustData is None:
            subplot1 = 131; subplot2 = 132;  subplot3 = 133;
        else:
            subplot1 = 111; subplot2 = 121; subplot3 = 122;

        pyl.figure(fig)
        ax = pyl.subplot(subplot2)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
        
        pyl.plot(N,Tb,'k')
        pyl.plot(Ntest,Ttest,'om')
        #pyl.legend(['Theoretical', 'Measured', 'Max Torque', 'Max Power'], loc = 'best')
        pyl.ylim( ymin = 0 )
        pyl.legend( ['Engine Model', 'Test Data'], loc = 'best' )
        
        pyl.subplot(subplot3)
        pyl.plot(N,Pb,'k')
        pyl.plot(Ntest,Ptest,'om')
        #pyl.legend(['Theoretical', 'Measured', "Max Torque", 'Max : ' + AsUnit(self.Pr(), self.PowerUnitName, "%3.1f") ], loc = 'best')
        pyl.ylim( ymin = 0 )            

        self.Prop.PlotTestData(fig)
        
        if subplot4 is not None:
            DynThrustData = npy.array(self.Prop.DynThrustData)
            
            VelocityData = self.ToUnumList(DynThrustData[:,2], FT/SEC) / (FT/SEC)
            V = npy.linspace(0,VelocityData[-1], 30)*(FT/SEC)
            
            Tbp = self.T(V)
            V   = V / (FT/SEC)
            Tbp = Tbp / self.Engine.ThrustUnit
            
            pyl.subplot(subplot4)
            pyl.plot(V,Tbp,'k')
            #pyl.legend(['Test Data', 'Propeller Model', "Propulsion Model" ], loc = 'best')
            pyl.ylim( ymin = 0, ymax = Tbp[0]*1.1 )            
            
################################################################################
if __name__ == '__main__':
    from ACPropeller import ACPropeller
    from ACEngine import ACEngine
    
    # Set Propeller properties
    Prop = ACPropeller()
    Prop.D          = 14.2*IN
    Prop.LenDi      = [0*IN, Prop.D]
    Prop.PitchAngle = 12*ARCDEG
    Prop.dAlpha     = 0*ARCDEG
    Prop.Solidity   = 0.0136
    Prop.RD         = 3/8
    Prop.AlphaStall = 14*ARCDEG
    Prop.Weight     = 0.3*LBF
    h=0*FT
    
    # Set Engine properties
    Engine  = ACEngine()
    Engine.Rbs          = 1.1
    Engine.Rla          = 3.5
    Engine.NumCyl       = 1
    Engine.NumRev       = 1
    Engine.CompRatio    = 9
    Engine.Vd           = 0.607*IN**3
    Engine.PistonSpeedR = 38.27*FT/SEC
    Engine.BMEPR        = 50.038*LBF/IN**2
    Engine.Rnmt         = 0.01
    Engine.Rtmt         = 1.5
    Engine.MEPtlmt      = 10.1526*LBF/IN**2
    Engine.SFCmt        = 1*PSFC
    Engine.Rnsfc        = 0.8
    Engine.A_F          = 16
    Engine.PS           = 1
    Engine.LWH          = [3*IN, 2*IN, 2*IN]
    Engine.Xat          = [1, 0.5, 0.5]
    Engine.Weight       = 1.5*LBF
    
    # Set Propulsion properties
    Propulsion = ACPropulsion(Prop,Engine)
    Propulsion.Alt  = 0*FT
    Propulsion.Vmax = 100*FT/SEC
    Propulsion.nV   = 20
    
    print 'Propulsion Weight :', AsUnit( Propulsion.Weight, "lbf")
    
    V1  = npy.array([0, 50])*FT/SEC
    V2  = npy.linspace(0,100,11)*FT/SEC
    N   = npy.linspace(1000,17000,17)*RPM
    
    Propulsion.PlotTPvsN(N, V1, 1)
    #Propulsion.PlotMatched(V2, 2)
    Propulsion.PlotFuelFlow(V2, 3)
    Propulsion.Draw(fig = 4)
    
    pyl.show()