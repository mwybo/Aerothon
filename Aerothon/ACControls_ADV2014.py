"""
A class for creating control derivatives
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACBase, ParamKeyError, Velocity, Force, Length, Area, Angle, MomentOfInertia, g, Unitless
from scalar.units import ARCDEG, SEC, RAD, FT, LBM, IN
from scalar.units import AsUnit
from ACTableBase import ACTableBase
from RootSolvers import BisectionSolver
from ACAVL import ACAVL
import AeroUtil
import cmath as math
import pylab as pyl
import numpy as npy
from scipy.optimize import brentq

################################################################################
class ACControlsError(Exception):
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg

################################################################################
class ACControls(ACBase, ACAVL):
    """
    A class for creating control derivatives
    
    Unless otherwise noted, referenced equations are from:
    Nelson, R. C., "Flight Stability and Automatic Control", 2nd ed., McGraw Hill, 1989.
    
    
    Attributes:
        Alpha2d       - Reference two dimensional angle of attack
        u0            - Flight Velocity
        Alt           - Flight altitude
        S             - Area of the main wing
        b             - Main wing span
        MAC           - Main wing mean aerodynamic chord
        Weight        - The weight of the aircraft
        Ixx, Iyy, Izz - Moments of inertia of the aircraft
    """

#===============================================================================
    def __init__(self, Aircraft):
        """
        Read in a AVL output file with control derivatives

        Inputs:
            Aircraft  - Aircraft to generate control derivatives for
        """
        super(ACControls, self).__init__()
        ACAVL.__init__(self)
        UnitList = self.UnitList
        
        self.__dict__['Alpha2d'] = 0*FT/SEC ;  UnitList['Alpha2d'] = Angle
        self.__dict__['u0']      = 0*FT/SEC ;  UnitList['u0']      = Velocity
        self.__dict__['Alt']     = 0*FT     ;  UnitList['Alt']     = Length
        self.__dict__['Weight']  = 0*FT     ;  UnitList['Weight']  = Force
        
        self.NoneList['StaticMargin']  = None;  UnitList['StaticMargin']  = Unitless

        self.__dict__['S']   = 0*IN**2     ;  UnitList['S']       = Area
        self.__dict__['b']   = 0*IN        ;  UnitList['b']       = Length
        self.__dict__['MAC'] = 0*IN        ;  UnitList['MAC']     = Length
        self.__dict__['Ixx'] = 0*LBM*IN**2 ;  UnitList['Ixx']     = MomentOfInertia
        self.__dict__['Iyy'] = 0*LBM*IN**2 ;  UnitList['Iyy']     = MomentOfInertia
        self.__dict__['Izz'] = 0*LBM*IN**2 ;  UnitList['Izz']     = MomentOfInertia
        
        self.__dict__['Deriv'] = []
                
        self.param.Aircraft = Aircraft

        #
        # Set initial values
        #
        self.SetToLO()
        self.SetToAircraft()

#===============================================================================
    def SetToLO(self):
        """
        Set all the attributes to the lift off conditions
        """
        Aircraft = self.param.Aircraft
       
        self.Alpha2d = Aircraft.GetAlpha2d_LO()
        self.u0      = Aircraft.GetV_LO()
        self.Alt     = Aircraft.GetAlt_LO()
        
#===============================================================================
    def SetToAircraft(self):
        """
        Set all the attributes to those of the aircraft
        """
        Aircraft = self.param.Aircraft
        
        self.S   = Aircraft.GetS()
        self.b   = Aircraft.Getb()
        self.MAC = Aircraft.GetMAC()
        self.Ixx = Aircraft.GetIxx()
        self.Iyy = Aircraft.GetIyy()
        self.Izz = Aircraft.GetIzz()

        self.Weight = Aircraft.EmptyWeight

#===============================================================================
    def q(self):
        """
        Returns the dynamic pressure for the controls analysis
        """
        return AeroUtil.q(self.Alt, self.u0)

#===============================================================================
    def rho(self):
        """
        Returns the density for the controls analysis
        """
        return AeroUtil.DensAlt(self.Alt)
    
#===============================================================================
    class _ACVarInit(ACBase):
        """
        Initializes the list of variables in _Coefficients
        """
        def __init__(self):
            super(ACControls._ACVarInit, self).__init__()            
            UnitList = self.UnitList
            for C in self._Coefficients:
                self.__dict__[C] = 0  ; UnitList[C] = Unitless
        
#===============================================================================
    class _ACAircraftDeriv(_ACVarInit, ACTableBase):
        """
        All control derivatives for the airplane that will get read in
        """
        _Coefficients = ('Cma', 'Cla', 'Clp', 'CYb', 'CYr', 'Cnp', 'Clb', 'Cnb', 'Cnr', 'Clr', 'Cmq', 'Xnp', 'Cmtot')
        
#=============================================================================== 
        def __init__(self, Control):
            """
            Inputs:
                Control - An ACControl instance this belongs to
            """
            super(ACControls._ACAircraftDeriv, self).__init__()
            
            Aircraft = Control.param.Aircraft
            UnitList = self.UnitList
            
            #
            # Creates an instance of each control surface in ControlSurfaces
            #
            self.param.ControlSurfaces = Aircraft.GetControlSurfaces()
            for ControlSurf in self.param.ControlSurfaces:
                self.__dict__[ControlSurf] = ACControls._ACControlSurfDeriv()
            
            #
            # Save off the Control for future use
            #
            self.param.Control = Control

#===============================================================================
        def Cza(self):
            """
            TODO: Documents Cza
            """
            Control  = self.param.Control
            Aircraft = Control.param.Aircraft
            Alpha2d  = Control.Alpha2d
            u0       = Control.u0
            return -(Aircraft.dCL_da() / (1/RAD) + Aircraft.CDTrim(Alpha2d, V = u0))
                
#===============================================================================
        def Yb(self):
            """
            Units: ft/s^2
            """
            Control = self.param.Control
            q    = Control.q()
            S    = Control.S
            CYb  = self.CYb
            Mass = Control.Weight / g
            return q * S * CYb / Mass
    
#===============================================================================
        def Yr(self):
            """
            Units: ft/s            
            """
            Control = self.param.Control
            q    = Control.q()
            S    = Control.S
            b    = Control.b
            u0   = Control.u0
            CYr  = self.CYr
            Mass = Control.Weight / g  
            return q * b * S * CYr / (2 * Mass * u0)
    
#===============================================================================
        def Nb(self):
            """
            Directional stability, or weathercock stability.
            
            Units: 1/s^2
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            b   = Control.b
            Cnb = self.Cnb
            Izz = Control.Izz
            return q * b * S * Cnb / Izz
    
#===============================================================================
        def Nr(self):
            """
            Yaw rate damping.
            
            Units: 1/s
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            b   = Control.b
            u0  = Control.u0
            Cnr = self.Cnr
            Izz = Control.Izz
            return q * b**2 * S * Cnr / (2 * Izz * u0)
    
#===============================================================================
        def Np(self):
            """
            Units: 1/s
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            b   = Control.b
            u0  = Control.u0
            Cnp = self.Cnp
            Izz = Control.Izz
            return q * b**2 * S * Cnp / (2 * Izz * u0)
    
#===============================================================================
        def Ndc(self, dc, ControlSurfName):
            """
            Units: 1/s^2
            """
            Control = self.param.Control
            q    = Control.q()
            S    = Control.S
            b    = Control.b
            Cndc = self.__dict__[ControlSurfName].Cn
            Izz  = Control.Izz
            return q * b * S * Cndc / Izz
    
#===============================================================================
        def Lb(self):
            """
            Dihedral effect
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            b   = Control.b
            Clb = self.Clb
            Ixx = Control.Ixx
            return q * b * S * Clb / Ixx

#===============================================================================
        def Lp(self):
            """
            Roll moment due to roll rate.
            
            Units: 1/s
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            b   = Control.b
            u0  = Control.u0
            Clp = self.Clp
            Ixx = Control.Ixx
            
            return q * S * b**2 * Clp / (2 * Ixx * u0)

#===============================================================================
        def Lr(self):
            """
            Roll moment due to yaw rate.
            
            Units: 1/s
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            b   = Control.b
            u0  = Control.u0
            Clr = self.Clr
            Ixx = Control.Ixx
            return q * b**2 * S * self.Clr / (2 * Ixx * u0)
    
#===============================================================================
        def Ldc(self, dc, ControlSurfName):
            """
            Calculates the roll moment due to a control surface deflection.
                 
            Units: 1/s^2
            
            Inputs:
                dc              - Control surface deflection
                ControlSurfName - Name of the control surface being deflected
            """
            Control = self.param.Control
            q    = Control.q()
            S    = Control.S
            b    = Control.b
            Cldc = self.__dict__[ControlSurfName].Cl
            Ixx  = Control.Ixx
            return q * b * S * Cldc / Ixx
    
#===============================================================================
        def Ma(self):
            """
            Units: 1/s^2
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            Cma = self.Cma
            Iyy = Control.Iyy
            MAC = Control.MAC
            
            if Cma >= 0.:
                raise ACControlsError("Positive Cm slope detected. Please check the Polar Slope approximations.")
            
            return Cma * q * S * MAC / Iyy
    
#===============================================================================
        def Madot(self):
            """
            Units: 1/s
            """
            Control = self.param.Control
            Aircraft = Control.param.Aircraft
            q      = Control.q()
            S      = Control.S
            u0     = Control.u0
            Iyy    = Control.Iyy
            MAC    = Control.MAC
            
            Xcg = None
            if Control.StaticMargin is not None:
                Xcg = Aircraft.Xnp() - MAC*Control.StaticMargin
            
            Cmadot = Aircraft.Cmadot(Xcg)

            return Cmadot * MAC**2 * q * S / (2 * u0 * Iyy)
    
#===============================================================================
        def Mq(self):
            """
            Units: 1/s
            """
            Control = self.param.Control
            q   = Control.q()
            S   = Control.S
            u0  = Control.u0
            Cmq = self.Cmq
            Iyy = Control.Iyy
            MAC = Control.MAC
            return Cmq * MAC**2 * q * S / (2 *u0 * Iyy)

#===============================================================================
        def Mdc(self, dc, ControlSurfName):
            """
            Returns the change in pitching moment due to a control surface deflection.
            
            Units: 1/s^2
            """
            Control = self.param.Control
            Cmdc = self.__dict__[ControlSurfName].Cm
            q    = Control.q()
            S    = Control.S
            MAC  = Control.MAC
            Iyy  = Control.Iyy
            return Cmdc * q * S * MAC / Iyy
                         
#===============================================================================
        def Zq(self):
            """n
            Units: ft/s
            """
            Control = self.param.Control
            Aircraft = Control.param.Aircraft
            q    = Control.q()
            S    = Control.S
            u0   = Control.u0
            Czq  = Aircraft.Czq()
            MAC  = Control.MAC
            Mass = Control.Weight / g 
            return -Czq * MAC * q * S / (2 * u0 * Mass)
    
#===============================================================================
        def Za(self):
            """
            Units: ft/s^2
            """
            Control = self.param.Control
            q    = Control.q()
            S    = Control.S
            Cza  = self.Cza()
            Mass = Control.Weight / g
            return Cza * q * S / Mass

#===============================================================================
        def Zu(self):
            """
            Units: 1/s
            """
            Control = self.param.Control
            Aircraft = Control.param.Aircraft
            q    = Control.q()
            S    = Control.S
            u0   = Control.u0
            CL0  = Aircraft.CLTrim(Control.Alpha2d)
            Cza  = self.Cza()
            Mass = Control.Weight / g
            return -2 * CL0 * q * S / (u0 * Mass)
                        
#===============================================================================
        def Xu(self):
            """
            Units: 1/s
            """
            Control = self.param.Control
            Aircraft = Control.param.Aircraft
            q    = Control.q()
            S    = Control.S
            u0   = Control.u0
            CD0  = Aircraft.CDTrim(Control.Alpha2d, V = u0)
            Mass = Control.Weight / g
            return -2 * CD0 * q * S / (u0 * Mass)
        
#===============================================================================
        def Pdc(self, dc, ControlSurfName):
            """   
            Returns the roll rate from control surface deflection.
            
            Eq. 5.7, pg. 182
            
            Inputs:
                dc              - Control surface deflection
                ControlSurfName - Name of the control surface being deflected
            """
            Control = self.param.Control
            b    = Control.b
            u0   = Control.u0
            Clp  = self.Clp
            Mass = Control.Weight / g
            Cldc = self.__dict__[ControlSurfName].Cl
            return -2 * u0 / b * Cldc / Clp * dc
                    
#===============================================================================
        def DutchRollFreq(self):
            """
            Returns the natural frequency of the Dutch roll motion.
            
            Eq. 5.47, pg. 198
            """
            Control = self.param.Control
            Yb = self.Yb()
            Yr = self.Yr()
            Nb = self.Nb()
            Nr = self.Nr()
            u0 = Control.u0
            return ( (Yb * Nr - Nb * Yr + u0 * Nb) / u0 )**0.5
    
#===============================================================================
        def DutchRollDamp(self):
            """
            Returns the damping ratio of the Dutch roll motion.
            
            Eq. 5.48, pg. 198
            """
            Control = self.param.Control
            Yb  = self.Yb()
            Nr  = self.Nr()
            u0  = Control.u0
            DRf = self.DutchRollFreq()
            return -((Yb + u0 * Nr) / (2 * DRf * u0))
    
#===============================================================================
        def ShortPeriodFreq(self):
            """
            Returns the frequency of the short period motion.
                        
            Table 4.3, pg. 155
            """
            Control = self.param.Control
            Ma = self.Ma()
            Mq = self.Mq()
            Za = self.Za()
            u0 = Control.u0
            return (Za * Mq / u0 - Ma)**0.5
    
#===============================================================================
        def ShortPeriodDamp(self):
            """
            Returns the damping ratio of the short period motion.
            
            Table 4.3, pg. 155
            """
            Control = self.param.Control
            Mq    = self.Mq()
            Madot = self.Madot()
            Za    = self.Za()
            u0    = Control.u0
            SPf   = self.ShortPeriodFreq()
            return -((Mq + Madot + Za / u0) / (2 * SPf))
    
#===============================================================================
        def LongPeriodFreq(self):
            """
            Returns the frequency of the long period (phugoid) motion.
            
            Table 4.3, pg. 155
            """
            Control = self.param.Control
            Zu = self.Zu()
            u0 = Control.u0
            return (-Zu * g / u0)**0.5
    
#===============================================================================
        def LongPeriodDamp(self):
            """
            Returns the damping ratio of the long period (phugoid) motion.
            
            Table 4.3, pg. 155
            """
            Xu  = self.Xu()
            LPf = self.LongPeriodFreq()
            return -(Xu / (2 * LPf))

#===============================================================================
        def SpiralMode(self):
            """
            Returns the frequency of the spiral mode.
            
            Units: 1/s
            
            Eq. 5.41, pg. 197
            """
            Lb = self.Lb()
            Lr = self.Lr()
            Nb = self.Nb()
            Nr = self.Nr()
            
#            if (Lb * Nr - Lr * Nb) / (1/(SEC**3)) > 0:
#                print "\n Spiral Mode Numerator is positive. (Probably Stable)"
            return (Lb * Nr - Lr * Nb) / Lb

#===============================================================================
        def RollReversal(self, AileronName = 'Aileron'):
            """
            Roll reversal occurs when this value drops below zero, Nicolai, pg 23-9.
            
            Units: Unitless
            """
            if not AileronName in self.__dict__:
                return 0
                #raise TypeError("Please provide the correct AileronName to RollReversal")

            Cnb = self.Cnb
            Clb = self.Clb 
            Cnda = self.__dict__[AileronName].Cn
            Clda = self.__dict__[AileronName].Cl
            
            return Cnb - Clb * Cnda/Clda
 
#===============================================================================
        def _CrossWind(self, Beta, dr, RudderName = 'Rudder'):
            """
            Cross wind calculation, Nicolai, pg 23.10.
            
            Beta       - Side slip angle
            dr         - Rudder deflection
            RudderName - Name of the rudder control surface
            
            Units: ARCDEG
            """
            if not RudderName in self.__dict__:
                raise TypeError("Please provide the correct RudderName to CrossWind")

            Cnb = self.Cnb
            Cndr = self.__dict__[RudderName].Cn
            
            return Cnb*Beta + Cndr*dr

#===============================================================================
        def RudderDeflection(self, Beta, RudderName = 'Rudder'):
            """
            Required rudder deflection to maintain straight flight in cross winds, Nicolai, pg 23.10.
            
            Beta       - Side slip angle
            RudderName - Name of the rudder control surface
            
            Units: ARCDEG
            """
            if not RudderName in self.__dict__:
                raise TypeError("Please provide the correct RudderName to RudderDeflection")

            def y(dr):
                return self._CrossWind(b, dr*ARCDEG, RudderName) / ARCDEG
            
            try:
                dr = []
                for b in Beta:
                    try:
                        dr.append( brentq(f = y, a = -30, b = 30, rtol = 1.e-2, maxiter = 500) )
                    except ValueError:
                        dr.append(-90)
                return npy.array(dr)*ARCDEG
            except:
                b = Beta
                return brentq(f = y, a = -30, b = 30, rtol = 1.e-2, maxiter = 100)*ARCDEG

#===============================================================================
        def CnbDynamic(self, Alpha3D):
            """
            If this term drops below zero, the aircraft will be susceptible to spin in hard nonrolling turns
            
            Alpha3D - Aircraft angle of attack
            
            Units: Unitless
            """
            Cnb = self.Cnb
            Clb = self.Clb
            Izz = self.param.Control.Izz
            Ixx = self.param.Control.Ixx

            def y(a):
                return Cnb - Clb*Izz/Ixx*math.tan(a/RAD).real

            try:
                CnbD = []
                for a in Alpha3D:
                    CnbD.append( y(a/RAD)  )
                return npy.array(CnbD)
            except:
                return y(Alpha3D/RAD)

                
#===============================================================================
        def SpiralDoubleTime(self):
            """
            Returns the time for the spiral motion to double in amplitude.
            
            Should be less than 12s or negative.
            """
            SM = self.SpiralMode()
            return math.log(2).real / SM
        
#===============================================================================
        def PitchFreq(self):
            """
            Returns the undamped natural frequency of the pitching motion.
            
            Units: 1/s
            
            Eq. 4.38, pg. 140
            """
            Ma = self.Ma()
            Mq = self.Mq()
            Za    = self.Za()
            u0   = self.param.Control.u0
            
            if Ma > 0.0/SEC**2:
                return 0.00001/SEC
            else:            
                return (-Ma + Za*Mq/u0)**0.5

#===============================================================================
        def PitchDamp(self):
            """
            Returns the damping ratio of the pitching motion.
            
            Eq. 4.39, pg. 140
            """
            Ma    = self.Ma()
            Mq    = self.Mq()
            Madot = self.Madot()
            Za    = self.Za()
            u0   = self.param.Control.u0
                        
            wn = self.PitchFreq()
            
            #Madot = 0
            
            if Ma > 0.0/SEC**2:
                return 0.00001
            else:
                return -((Mq + Madot + Za/u0) / (2 * wn ))

#===============================================================================
        def PitchPhaseAng(self):
            """
            Returns the phase angle of the pitching motion.
            
            Eq. 4.45, pg. 141
            """
            Pd = self.PitchDamp()
            return math.atan(math.sqrt(1 - Pd**2)/ Pd)

#===============================================================================
        def PitchDoubleTime(self):
            """
            Returns the time for the pitch motion to double in amplitude.
            
            Eq. 4.74, pg. 143
            """
            Pf = self.PitchFreq()
            Pd = self.PitchDamp()
            return math.log(2).real / (-Pd * Pf)
        
#===============================================================================
        def PitchResponseDueToElevator(self, de, T, ElevatorName):
            """
            Returns the normalized change in alpha due to a step input of the elevator.
            
            Eq. 4.45, pg. 141
            """
            PPa      = self.PitchPhaseAng()
            PF       = self.PitchFreq()
            PD       = self.PitchDamp()
            
#            PD = 0.5
#            PF = 8.98/SEC
#            PPa = math.atan(math.sqrt(1 - PD**2)/ PD)
            
            # Split the equation for aesthetic purposes
            dA = []
            for t in T:
                Num   = math.e**( -PD * PF * t )
                Denom = math.sqrt(1. - PD**2)
                Angle = math.sqrt(1. - PD**2) * PF * t + PPa
                SinA  = math.sin(Angle)
                dA.append( max(0, (1. - Num / Denom * SinA).real) )
            
            if len(dA)==1:
                return dA[0]
            else:
                return npy.array(dA)
            
#===============================================================================
        def DAlphaFromElevator(self, de, t, ElevatorName):
            """
            Returns the time to reach a unit deflection of the elevator.
            """
            Ma         = self.Ma()
            Mde        = self.Mdc(de, ElevatorName)
            DAlphaTrim = -(Mde * de) / Ma

            return self.PitchResponseDueToElevator(de, t, ElevatorName) * DAlphaTrim

#===============================================================================
        def PlotPitchResponseDueToElevator(self, de, t, ElevatorName, fig=1):
            """
            Plots the time to reach a unit deflection of the elevator.
            """
            
            dA = self.PitchResponseDueToElevator(de, t, ElevatorName)
            pyl.figure(fig)
            pyl.plot(t / SEC, dA)
            pyl.xlabel('Time (sec)')
            pyl.ylabel('Elevator step input response')

#===============================================================================
        def PlotStability(self, fig=1, RudderName = 'Rudder'):
            """
            Plots the time to reach a unit deflection of the elevator.
            """
            
            Alpha3D = npy.linspace(-20, 25, 30)*ARCDEG
            CnbDy = self.CnbDynamic(Alpha3D)
            Beta = npy.linspace(0, 15, 30)*ARCDEG
            dr = self.RudderDeflection(Beta, RudderName)
            
            pyl.figure(fig)
            pyl.subplot(121)
            pyl.plot(Alpha3D/ARCDEG, CnbDy)
            pyl.xlabel('Aircraft Angle of Attack (deg)')
            pyl.ylabel(r'$C_{n\beta}$ Dynamic')
            pyl.title('If below zero, aircraft susceptible\n to spin in hard nonrolling turns')
            pyl.grid(b=True)

            pyl.subplot(122)
            pyl.plot(Beta/ARCDEG, dr/ARCDEG)
            pyl.xlabel(r'Aircraft Sideslip $\beta$ (deg)')
            pyl.ylabel('Rudder Deflection, $\delta_r$ (deg)')
            pyl.title('Required rudder deflection\n for straight flight')
            pyl.grid(b=True)
            pyl.gca().yaxis.set_label_coords(-0.1,0.5)
                            
#===============================================================================
        def TimeForPitchChange(self, de, ElevatorName):
            """
            Returns the time to reach 90% of the defined alpha change.
            """
            def DAlpha(t):
                return self.PitchResponseDueToElevator(de, t * SEC, ElevatorName) - 0.9
            t = BisectionSolver(0, 10, DAlpha, 0.0001 , 100)

            return t * SEC
                
#===============================================================================
        def RollDueToAileron(self, da, AileronName, V = None):
            """
            Returns the steady state roll rate due to an aileron deflection of da.
            
            Pg. 184
            
            Inputs:
                da          -- Aileron deflection angle
                AileronName -- Name of the aileron
                V           -- Velocity used in the calculation
            """
            Control = self.param.Control
            
            if V is None:
                V  = Control.u0
            
            b    = Control.b
            Clda = self.__dict__[AileronName].Cl
            Clp  = self.Clp
            
            return -2*V*Clda*da/(b*Clp)*RAD

#===============================================================================
        def StabilityTable(self, fig = 1):
            """
            Displays the table of stability derivatives and their appropriate
            stable ranges.
            """
            
            pyl.figure(fig)
            
            AsUnit = '%1.2f'
            
            colLabels = ('Phugoid', 'Short Period', 'Pure Pitching', 'Spiral Mode', 'Dutch Roll')
            rowLabels = ('Required',r'$\omega$ (1/s)',r'$\zeta$')
            #cellText  = [[r'0.04 < $\zeta$',r'0.35 < $\zeta$ < 1.3',r'1.0 < $\zeta$',r'$\zeta$ 12s or negative',r'0.08 < $\zeta$; 0.4 < $\omega$']]
            cellText  = [[r'0.04 < $\zeta$',r'0.35 < $\zeta$ < 1.3',r'1.0 < $\zeta$',r'"$\zeta$" >20s or <0',r'0.19 < $\zeta$; 0.4 < $\omega$']]
                         

            Freq = [self.LongPeriodFreq() / (1/SEC),
                    self.ShortPeriodFreq() / (1/SEC),
                    self.PitchFreq() / (1/SEC),
                    self.SpiralMode() / (1/SEC),
                    self.DutchRollFreq() / (1/SEC)]
            Damp = [self.LongPeriodDamp(),
                    self.ShortPeriodDamp(),
                    self.PitchDamp(),
                    self.SpiralDoubleTime() / (SEC),
                    self.DutchRollDamp()]

            cellText.append( [AsUnit % Freq[i] for i in range(len(Freq))] )
            cellText.append( [AsUnit % Damp[i] for i in range(len(Damp))] )
            colWidths = self._GetTableColWidths(rowLabels, cellText, colLabels)
           
            #
            # Mark any values outside of the appropriate range
            #
            cellColours = []
            for r in cellText:
                cellColours.append(['w' for i in range(len(r))])
            
            #
            # Frequencies
            #
#            cellColours[1][3] = 'w' if 0.0 < Freq[3] else 'r'
            cellColours[1][4] = 'w' if 0.4 < Freq[4] else 'r'
            
            #
            # Damping rates
            #
            cellColours[2][0] = 'w' if 0.04 < Damp[0] else 'r'
            cellColours[2][1] = 'w' if 0.35 < Damp[1] and Damp[1] < 1.3 else 'r'
            cellColours[2][2] = 'w' if 1.0  < Damp[2] else 'y' if 0.8 < Damp[2] else 'r'
            cellColours[2][3] = 'w' if 12.0 < Damp[3] or Damp[3] < 0.0 else 'r'
            cellColours[2][4] = 'w' if 0.08 < Damp[4] else 'r'

            
            pyl.axes([0.01, 0.01, 0.99, 0.95])
            pyl.axis('off')
            table = pyl.table(cellText = cellText, cellColours=cellColours,
                              rowLabels = rowLabels, #rowColours=colours,
                              colWidths = colWidths,
                              colLabels = colLabels, colLoc = 'center',
                              loc = 'upper right')
            
            table.auto_set_font_size(False)
            table.set_fontsize(11)
            table.scale(1., 1.7)
            
            #------------------------------------------------------------------
            # Second table
            #------------------------------------------------------------------
            
            AsUnit = '%1.4f'
            
            colLabels = ['Roll Reversal', r'$C_{n\beta}$ Dynamic']
            rowLabels = ['Required',r'Value']
            cellText  = [[r'0 < ', r'0 < ']]
                         
            Values = [self.RollReversal(),
                      (self.CnbDynamic(-20*ARCDEG),self.CnbDynamic(25*ARCDEG)) ]

            cellText.append( [AsUnit % v for v in Values[:-1]] )
            cellText[-1].append( (AsUnit % Values[-1][0]) + ', '  + (AsUnit % Values[-1][1]) )
            colWidths = self._GetTableColWidths(rowLabels, cellText, colLabels)
           
            #
            # Mark any values outside of the appropriate range
            #
            cellColours = []
            for r in cellText:
                cellColours.append(['w' for i in range(len(r))])
                        
            #
            # Mark them as red
            #
            cellColours[1][0] = 'w' if 0.00 < Values[0] else 'r'
            cellColours[1][1] = 'w' if 0.00 < Values[-1][0] and 0.00 < Values[-1][1] else 'r'

            table = pyl.table(cellText = cellText, cellColours=cellColours,
                              rowLabels = rowLabels, #rowColours=colours,
                              colWidths = colWidths,
                              colLabels = colLabels, colLoc = 'center',
                              loc = 'center right')
            
            table.auto_set_font_size(False)
            table.set_fontsize(11)
            table.scale(1., 1.7)            

            #------------------------------------------------------------------
            # Third table
            #------------------------------------------------------------------
            
            AsUnit = '%1.3f'
            
            colLabels = [r'$C_{m\alpha}$', r'$C_{mq}$', r'$C_{m \delta e}$', r'$C_{n\beta}$', r'$C_{n \delta r}$', r'$C_{l \beta}$', r'$C_{l \delta a}$']
            rowLabels = [' ']
            cellText  = []
            
            try:
                Aileron_Cl = self.Aileron.Cl
            except ParamKeyError:
                Aileron_Cl = 0
            
            Values = [self.Cma,
                      self.Cmq,
                      self.Elevator.Cm,
                      self.Cnb,
                      self.Rudder.Cn,
                      self.Clb,
                      Aileron_Cl ]

            cellText.append( [AsUnit % v for v in Values] )
            colWidths = self._GetTableColWidths(rowLabels, cellText, colLabels)
           
            #
            # Just set all colors to white
            #
            cellColours = []
            for r in cellText:
                cellColours.append(['w' for i in range(len(r))])
                        
            table = pyl.table(cellText = cellText, cellColours=cellColours,
                              rowLabels = rowLabels, #rowColours=colours,
                              colWidths = colWidths,
                              colLabels = colLabels, colLoc = 'center',
                              loc = 'lower right')
            
            table.auto_set_font_size(False)
            table.set_fontsize(11)
            table.scale(1., 1.7)            

#===============================================================================
        def ControlsTables(self, fig = 1):
            """
            Display the control derivatives for the control surfaces
            """
            pyl.figure(fig)
            
            pyl.axes([0.01, 0.01, 0.99, 0.95])
            pyl.axis('off')
            
            loc = ['upper right', 'center right', 'lower right']
            
            n = 0
            for ControlName in self.param.ControlSurfaces:
                
                Control = self.__dict__[ControlName]
                
                AsUnit = '%1.6f'
                
                colLabels = ACControls._ACControlSurfDeriv._Coefficients
                rowLabels = [ControlName]
                cellText  = []

                cellText.append( [AsUnit % (Control.__dict__[Deriv]/(1/ARCDEG)) for Deriv in ACControls._ACControlSurfDeriv._Coefficients] )
                colWidths = self._GetTableColWidths(rowLabels, cellText, colLabels)
               
                #
                # Mark any values outside of the appropriate range
                #
                cellColours = []
                for r in cellText:
                    cellColours.append(['w' for i in range(len(r))])
                
                #
                # Frequencies
                #
    #            cellColours[1][3] = 'w' if 0.0 < Freq[3] else 'r'
                #cellColours[1][4] = 'w' if 0.4 < Freq[4] else 'r'
                
                #
                # Damping rates
                #
                #cellColours[2][0] = 'w' if 0.04 < Damp[0] else 'r'
                #cellColours[2][1] = 'w' if 0.35 < Damp[1] and Damp[1] < 1.3 else 'r'
                #cellColours[2][2] = 'w' if 1.0  < Damp[2] else 'y' if 0.8 < Damp[2] else 'r'
                #cellColours[2][3] = 'w' if 12.0 < Damp[3] or Damp[3] < 0.0 else 'r'
                #cellColours[2][4] = 'w' if 0.08 < Damp[4] else 'r'
    
                
                table = pyl.table(cellText = cellText, cellColours=cellColours,
                                  rowLabels = rowLabels, #rowColours=colours,
                                  colWidths = colWidths,
                                  colLabels = colLabels, colLoc = 'center',
                                  loc = loc[n])
                
                table.label = ControlName
                table.auto_set_font_size(False)
                table.set_fontsize(11)
                table.scale(1., 1.7)
                
                fig += 1
                n += 1
            
#===============================================================================
    class _ACControlSurfDeriv(_ACVarInit):
        """
        All control derivatives for the control surfaces that will get read in
        """
        _Coefficients = ('CL', 'CY', 'Cl', 'Cm', 'Cn')
                                                                               
#===============================================================================
    def _ReadAVLControlSurfOutput(self, ControlName, ControlDeriv, AVLOutputFile):
        """
        Reads the controls surface derivatives from the AVL output file

        Inputs
            ControlDeriv - ACDerivaties class with a list of derivatives to read in
        """

        try:

            f = open(AVLOutputFile,'r')
            #
            # Find the line containing the derivative
            #
#            line = f.readline()
#            while line.find('derivatives') == -1:
#                line = f.readline()
#                if len(line) == 0:
#                    break #Hit the end of the file 
            #
            # Find the second occurrence of the control name
            #
            for n in range(2):
                line = f.readline()
                while line.find(ControlName) == -1:
                    line = f.readline()
                    if len(line) == 0:
                        break #Hit the end of the file 

            if line.find(ControlName) == -1:
                raise "Could not find " + ControlName + " in " + AVLOutputFile
            
            line  = line[line.find(ControlName)+len(ControlName):]
            line  = line.lstrip(' ')
            ControlID  = line[:line.find(' ')]
            f.close()
        except IOError, (errno, strerror):
            print "While reading '", AVLOutputFile, "'"
            raise

        class CDeriv():
            pass

        Derivs = CDeriv()
        Derivs._Coefficients = []
        for key in ControlDeriv._Coefficients:
            Derivs._Coefficients.append(key + ControlID)

        ACAVL._ReadAVLOutput(self, Derivs, AVLOutputFile)

        for key in ControlDeriv._Coefficients:
            #All control derivatives are in units of 1/ARCDEG
            ControlDeriv.__dict__[key] = Derivs.__dict__[key + ControlID] * 1/ARCDEG 

#===============================================================================
    def AddRun(self, RunName, AVLInputFile, WriteAVLInput = True):
        """
        
        """
        ACAVL.AddRun(self, RunName, AVLInputFile)
        
        if WriteAVLInput:
            Xcg = None
            if self.StaticMargin is not None:
                Aircraft = self.param.Aircraft
                Xcg = Aircraft.Xnp() - Aircraft.Wing.MAC()*self.StaticMargin
        
            self.param.Aircraft.WriteAVLAircraft(self.RunDir + AVLInputFile, Xcg)

#===============================================================================
    def ReadAVLFiles(self):
        """
        
        """
        self.Deriv = ACAVL.ReadAVLFiles(self, ACControls._ACAircraftDeriv, Control = self)
        
#===============================================================================
    def _ReadAVLOutput(self, Derivatives, AVLOutputFile):
        """
        Reads the AVL output file

        Inputs
            Derivatives - ACDerivaties class with a list of derivatives to read in
        """
        ACAVL._ReadAVLOutput(self, Derivatives, AVLOutputFile)
        
        for ControlName in Derivatives.param.ControlSurfaces:
            self._ReadAVLControlSurfOutput(ControlName, Derivatives.__dict__[ControlName], AVLOutputFile)

################################################################################
if __name__ == '__main__':
    from TestAircraft import aircraft
    
#    Deriv = ACControls._ACAircraftDeriv(aircraft)
#    me = ACControls(aircraft)
#    me._ReadAVLOutput(Deriv, 'AVLDeriv.txt')
#    for key in Deriv._Coefficients:
#        print key, ' = ', Deriv.__dict__[key]
#
#    Aileron = ACControls._ACControlSurfDeriv()
#    me._ReadAVLControlSurfOutput('Aileron', Aileron, 'AVLDeriv.txt')
#    for key in Aileron._Coefficients:
#        print key, ' = ', Aileron.__dict__[key]

    Controls = ACControls(aircraft)
    Controls.AddRun('Test', 'AVLInputFile.dat', WriteAVLInput = False)
    Controls.Test.DumpStability('AVLDeriv.txt')
    Controls.Test.Exit()
    #Controls.ExecuteAVL()
    
    Controls.ReadAVLFiles()
    
    Deriv = Controls.Deriv[0]

    print
    print
    print 'Yb = ', Deriv.Yb()
    print 'Yr = ', Deriv.Yr()
    print 'Nb = ', Deriv.Nb()
    print 'Nr = ', Deriv.Nr()
    print 'Np = ', Deriv.Np()
    print 'Lb = ', Deriv.Lb()
    print 'Lr = ', Deriv.Lr()
    print 'Ma = ', Deriv.Ma()
    print 'Mq = ', Deriv.Mq()
    print 'Za = ', Deriv.Za()
    print 'Xu = ', Deriv.Xu()
    print 'Ndc = ', Deriv.Ndc(3 * ARCDEG ,'Aileron')
    print 'Ldc = ', Deriv.Ldc(3 * ARCDEG ,'Aileron')
    print 'Pdc = ', Deriv.Pdc(3 * ARCDEG ,'Rudder')
    print
    print 'Dutch Roll Freq: ', Deriv.DutchRollFreq()
    print 'Dutch Roll Damp: ', Deriv.DutchRollDamp()    
    print
    print 'Short Period Freq: ', Deriv.ShortPeriodFreq()
    print 'Short Period Damp: ', Deriv.ShortPeriodDamp()
    print
    print 'Long Period Freq: ', Deriv.LongPeriodFreq()
    print 'Long Period Damp: ', Deriv.LongPeriodDamp()    
    print
    print 'Spiral Mode: ', Deriv.SpiralMode()
    print 'Spiral double time: ', Deriv.SpiralDoubleTime()
    print
    print 'Pitch Freq: ', Deriv.PitchFreq()
    print 'Pitch Damp: ', Deriv.PitchDamp()
    print 'Pitch Phase Ang: ', Deriv.PitchPhaseAng()
    print 'Pitch Double Time: ', Deriv.PitchDoubleTime()
    print
    print 'Madot = ', Deriv.Madot()
    print 'Zq = ', Deriv.Zq()
    print
    print 'RollReversal = ', Deriv.RollReversal()
    print 'Rudder Deflection = ', AsUnit(Deriv.RudderDeflection(11.5*ARCDEG), 'deg')
    
#    print 'DAlphaFromElevator = ', Deriv.DAlphaFromElevator(3, 5 * SEC,'Elevator')

#    print 'TimeForPitchChange = ', Deriv.TimeForPitchChange(3 * ARCDEG, 'Elevator')
    
    print 'Steady state roll rate = ', AsUnit( Deriv.RollDueToAileron(15 * ARCDEG, 'Aileron'), "deg/s" )
    
    t = npy.linspace(0, 5, 20)*SEC
    Deriv.StabilityTable(2)
    Deriv.PlotStability(fig=3)
    Deriv.PlotPitchResponseDueToElevator(fig=4, de=5*ARCDEG, t=t, ElevatorName = 'Elevator')
    
    
#    Controls = ACControls(aircraft)
#    Controls.AddRun('Run1','AVLAircraft1.txt')
#    Controls.Run1.DumpStability('STFile1.txt')
#    Controls.Run1.Exit()
#    
#    Controls.ExecuteAVL()
#
#    Controls.ReadAVLFiles()
#    
#    Controls.Deriv[0].StabilityTable(2)

    pyl.show()
    
    print 'End of script'
