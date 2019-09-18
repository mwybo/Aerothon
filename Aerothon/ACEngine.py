"""
This is the class for analyzing the internal combustion engine
Used for Propulsion analysis
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACBase, Unitless, Volume, Velocity, Pressure, Power_SFC, PerAngle, Force, Power, Torque
from ACMass import ACMassBox, ACMassCyl
from scalar.units import RPM, IN, ARCDEG, RAD, ATM, FT, SEC, LBF, LBM, PSFC, PSI, GAL, HP, HR, OZF, M
from scalar.units import AsUnit
import cmath as math
import numpy as npy
import pylab as pyl
import AeroUtil as AUtil

################################################################################
class ACEngineError(Exception):
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg

################################################################################
class ACEngine(ACMassBox):
    """
    Engine class
    
    Attributes:
        Rbs          - Ratio of engine bore diameter to stroke length (default 1.1)
        Rla          - Ratio of engine connecting rod length to crank radius (default 3.5)
        NumCyl       - Number of Cylinders (default 1)
        NumRev       - Number of Revelutions per Power Stroke (Two Stroke=1, Four Stroke=2) (default 1)
        CompRatio    - Compression Ratio
        Vd           - Displacement Volume (default 0.607 in^3)
        
        -----------------------------------------------------------
        Engine Family Correlating parameters
        - Use these to correlate the model to match a particular engine family
        -----------------------------------------------------------
        PistonSpeedR - Average Rated Piston Speed (Sp)                   : Controls Rated RPM
        BMEPR        - Brake Mean Effective Pressure (Rated)             : Controls Rated Torque
        Rnmt         - Ratio of RPM at max torque to RPM at at max power : Controls Torque curve shape
        Rtmt         - Ratio of Torque at max torque to rated torque     : Controls Torque curve shape
        -----------------------------------------------------------
        
        MEPtlmt      - Total loss for mean effective pressure at max torque (default 10.1526 psi)
        SFCmt        - Specific Fuel Consumption (SFC) at Max Torque 
        Rnsfc        - Ratio RPM at min SFC to RPM at rated torque (max power) (default 0.8)
        A_F          - Air to Fuel Ratio (default 16)
        PS           - Power setting (Throttle setting) (default 1)
        
        TorqueUnit     - Unit for torque (default in*ozf)
        TorqueUnitName - Name of unit for thrust (default in*ozf)
        ThrustUnit     - Unit for thrust (default lbf)
        ThrustUnitName - Name of unit for thrust (default lbf)
        PowerUnit      - Unit for power (default hp)
        PowerUnitName  - Name of unit for power (default hp)
        
        GearRatio    - Gear ratio for a gear box (default 1)
        GearEff      - Gear Box efficiency (default 1.0)
        
        NoseCone     - An ACMassCyl representing the nose cone (for drawing purposes)
        Muffler      - An ACMassCyl representing the muffler (for drawing purposes)
        
        TestData - An array of test data
                    #            RPM,        Torque
                    TestData = [(7000  *RPM, 107.35 *IN*OZF),
                                (8500  *RPM, 104.24 *IN*OZF),
                                (9500  *RPM, 101.13 *IN*OZF)]
    """
    
#===============================================================================    
    def __init__(self):
        super(ACEngine,self).__init__()
        
        self.name = 'Engine'
        
        UnitList = self.UnitList
        self.__dict__['Rbs']            = 1.1               ; UnitList['Rbs']          = Unitless
        self.__dict__['Rla']            = 3.5               ; UnitList['Rla']          = Unitless
        self.__dict__['NumCyl']         = 1                 ; UnitList['NumCyl']       = Unitless
        self.__dict__['NumRev']         = 1                 ; UnitList['NumRev']       = Unitless
        self.__dict__['CompRatio']      = 9                 ; UnitList['CompRatio']    = Unitless
        self.__dict__['Vd']             = 0.607*IN**3       ; UnitList['Vd']           = Volume
        self.__dict__['PS']             = 1                 ; UnitList['PS']           = Unitless
        self.__dict__['PistonSpeedR']   = 0*FT/SEC          ; UnitList['PistonSpeedR'] = Velocity
        self.__dict__['BMEPR']          = 0*LBF/IN**2       ; UnitList['BMEPR']        = Pressure
        self.__dict__['Rnmt']           = 0.01              ; UnitList['Rnmt']         = Unitless
        self.__dict__['Rtmt']           = 0                 ; UnitList['Rtmt']         = Unitless
        self.__dict__['MEPtlmt']        = 10.1526*LBF/IN**2 ; UnitList['MEPtlmt']      = Pressure
        self.__dict__['SFCmt']          = 0*PSFC            ; UnitList['SFCmt']        = Power_SFC
        self.__dict__['Rnsfc']          = 0.8               ; UnitList['Rnsfc']        = Unitless
        self.__dict__['A_F']            = 16                ; UnitList['A_F']          = Unitless
        self.__dict__['PS']             = 1                 ; UnitList['PS']           = Unitless

        self.__dict__['TorqueUnit'] = IN*OZF                ; UnitList['ThrustUnit'] = Torque
        self.__dict__['TorqueUnitName'] = 'in*ozf'

        self.__dict__['ThrustUnit'] = LBF                   ; UnitList['ThrustUnit'] = Force
        self.__dict__['ThrustUnitName'] = 'lbf'

        self.__dict__['PowerUnit'] = HP                     ; UnitList['PowerUnit'] = Power
        self.__dict__['PowerUnitName'] = 'hp'          

        self.__dict__['GearRatio']      = 1                 ; UnitList['GearRatio']    = Unitless
        self.__dict__['GearEff']        = 1                 ; UnitList['GearEff']      = Unitless

        self.__dict__['NoseCone']       = ACMassCyl()
        self.__dict__['Muffler']        = ACMassCyl()
        self.__dict__['Mount']          = ACMassBox()
        self.__dict__['Cylinder']       = ACMassCyl()
        self.__dict__['Crankcase']      = ACMassCyl()
                        
        self.__dict__['TestData']       = None
        
        self.__dict__['DrawDetail']     = False
        
        self.NoseCone.LenDi = [1.5*IN,0.75*IN] #TODO: This is just a guess
        self.Muffler.LenDi = [6*IN,1.5*IN] #TODO: This is just a guess

#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """
        Draws the propeller hub and muffler

        Inputs:
            fig - Integer number of figure to draw into
            top - Subplot for top view
            side - Subplot for side view
            front - Subplot for front view
        """
        if self.DrawDetail is False:
            super(ACEngine,self).Draw(fig, top, side, front)
        else:
            self.Muffler.Draw(fig, top, side, front)
            self.Mount.Draw(fig, top, side, front)
            self.Cylinder.Draw(fig, top, side, front)
            self.Crankcase.Draw(fig, top, side, front)
        
        self.NoseCone.Draw(fig, top, side, front)
        
#===============================================================================
    def Length(self):
        """
        Computes the length of the engine
        """
        return self.LWH[0] + self.NoseCone.LenDi[0]

#===============================================================================
    def PropX(self):
        """
        Computes the position of the propeller on the engine
        """
        X     = self.GetX().copy()
        X[0] -= self.LWH[0]*self.Axis[0]*self.Xat[0]
        X[0] -= self.NoseCone.LenDi[0] * 0.5
        
        return X

#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACEngine, self).Refresh()

        self.param.refreshing = True
        
        self.NoseCone.Symmetric = self.Symmetric
        self.NoseCone.Refresh()
        
        #
        #  Draw component parts of engine....
        #
        L = self.LWH[0]
        W = self.LWH[1]
        H = self.LWH[2]
        
        self.Mount.Axis = [-1,0,0]
        self.Mount.X    = self.GetX()
        self.Mount.LWH  = [L * 0.6 , W * 1.0 , H * 0.1]
        
        self.Crankcase.Axis  = [-1,0,0]
        self.Crankcase.X     = self.GetX() + npy.array([-0.2 * L / IN , 0 , 0 ])*IN
        self.Crankcase.LenDi = [L * (1 - 0.2), W * 0.75]
        self.Crankcase.Refresh()
        
        self.Cylinder.Axis   = [0,0,1]
        self.Cylinder.X      = self.GetX() + npy.array([-0.5 * L / IN , 0 , 0 ])*IN
        self.Cylinder.LenDi  = [H * 0.75, W * 0.8]
        self.Cylinder.Refresh()
        
        self.NoseCone.Axis   = [-1,0,0]
        self.NoseCone.X      = self.GetX() + npy.array([-L / IN,0,0])*IN
        self.NoseCone.Refresh()
        
        self.Muffler.X       = self.GetX() + npy.array([-0.90 * L / IN, H / IN * 0.75 , H / IN * 0.6])*IN
        self.Muffler.Refresh()
        
        self.param.refreshing = False

#===============================================================================
    def _CalcRla(self,l,a):
        """
        Inputs:
            l - connecting rod length
            a - crank radius
        """
        newRla = l/a
        
        self.Rla = newRla
    
#===============================================================================
    def B(self):
        """
        Bore calculation from displacement volume and bore to stroke ratio
        """
        Vd  = self.Vd
        rbs = self.Rbs
        nc  = self.NumCyl
        
        return (Vd*4/math.pi*rbs/nc)**(1/3)
        
#===============================================================================
    def L(self):
        """
        Stroke length calculation
        """
        rbs = self.Rbs  
        B   = self.B()
        
        return B/rbs
    
#===============================================================================
    def Nr(self):
        """
        Rated crankshaft rotational speed
        """
        spr = self.PistonSpeedR 
        L   = self.L()
        
        return spr/(2*L) * 2.*math.pi / self.GearRatio #Multiply by 2 pi to put it in RPM
    
#===============================================================================
    def Tr(self):
        """ 
        Torque at rated max power function of mean effective pressure at rated speed and displacement volume
        """
        bmepr = self.BMEPR
        vd    = self.Vd
        nr    = self.NumRev
        
        return self.GearRatio*bmepr*vd/(2*math.pi*nr)
    
#===============================================================================
    def Pr(self):
        """ 
        Power at rated max power function of mean effective pressure at rated speed and displacement volume
        """  
        Nr  = self.Nr()
        Tr  = self.Tr()
        
        return Nr*Tr * self.GearEff
    
#===============================================================================
    def Nmt(self):
        """ 
        RPM at max torque
        """
        Nr      = self.Nr()
        rnmt    = self.Rnmt
          
        return Nr*rnmt / self.GearRatio

#===============================================================================
    def Tm(self):
        """ 
        Max torque
        """
        rtmt    = self.Rtmt
        Tr      = self.Tr()
        
        return Tr*rtmt*self.GearRatio

#===============================================================================
    def BMEPMT(self):
        """ 
        Max brake mean effective pressure
        """
        Tm  = self.Tm()
        nr  = self.NumRev
        vd  = self.Vd
        
        return Tm*2*math.pi*nr/vd
    
#===============================================================================
    def BMEPsl(self,N):
        """
        Brake mean effective pressure at sea level
        
        Input:
            N - RPM
        """
        BMEPMT  = self.BMEPMT()
        bmepr   = self.BMEPR
        Nmt     = self.Nmt()
        rnmt    = self.Rnmt
        
        return BMEPMT+(bmepr-BMEPMT)*((N/Nmt-1)**2/(1/rnmt-1)**2)
 
#===============================================================================
    def IMEPMT(self):
        """
        Max indicated mean effective pressure at max torque
        """
        BMEPMT  = self.BMEPMT()
        meptlmt = self.MEPtlmt
        
        return BMEPMT + meptlmt
    
#===============================================================================
    def MEPtlsl(self,N):
        """
        Assumed total losses at max rated torque
        
        Inputs:
            N - RPM
        """
        IMEPMT  = self.IMEPMT()
        BMEPsl  = self.BMEPsl(N)
        
        return IMEPMT - BMEPsl
        
 
 #===============================================================================
    def MEPtf(self,N):
        """
        Total friction losses in terms of mean effective pressure
        
        Inputs:
            N - RPM
        """       
        Nr  = self.Nr()
                
        return 0.8*(0.97+0.27*(N/Nr)+0.19*(N/Nr)**2)*ATM
 
#===============================================================================   
    def MEPtp(self,alt,N):
        """
        Total pumping losses in terms of mean effective pressure
        
        Inputs:
            N   - RPM
            alt - altitude
        """
        ps      = self.PS
        MEPtlsl = self.MEPtlsl(N)
        MEPtf   = self.MEPtf(N)
        
        alt = self.MakeList(alt)
        
        sigma = []
        for h in alt:
            sigma.append(AUtil.DensityRatioAlt(h))
            
        if len(sigma) == 1:
            sigma = sigma[0]
        else:        
            sigma = npy.array(sigma)
            
        return (MEPtlsl-MEPtf)*ps*sigma
    
#===============================================================================   
    def IMEP(self, alt, N):
        """
        Indicated(gas)mean effective pressure
        
        Inputs:
            alt - altitude
            N   - RPM
        """
        BMEP    = self.BMEP(alt, N)
        MEPtf   = self.MEPtf(N)
        MEPtp   = self.MEPtp(alt, N)
        
        return BMEP + MEPtp + MEPtf
    
#===============================================================================   
    def BMEP(self, alt, N):
        """
        Brake mean effective pressure
        
        Inputs:
            alt - altitude
            N   - RPM
        """
        ps      = self.PS
        BMEPsl = self.BMEPsl(N)
        
        alt = self.MakeList(alt)
        
        sigma = []
        for h in alt:
            sigma.append(AUtil.DensityRatioAlt(h))
        
        if len(sigma) == 1:
            sigma = sigma[0]
        else:        
            sigma = npy.array(sigma)
        
        return BMEPsl*ps*sigma

#===============================================================================   
    def Ti(self,alt, N):
        """
        Indicated Torque
        
        Inputs:
            alt - altitude
            N   - RPM
        """
        IMEP  = self.IMEP(alt, N)
        MEPtf = self.MEPtf(N)
        MEPtp = self.MEPtp(alt, N)
        vd   = self.Vd
        nr   = self.NumRev;
        
        return self.GearRatio*IMEP*vd/(2. * math.pi*nr)

#===============================================================================   
    def Pi(self,alt,N):
        """
        Indicated Power
        
        Inputs:
            N   - RPM
            alt - altitude
        """
          
        return N*self.Ti(alt, N) * self.GearEff
    
#===============================================================================   
    def Tb(self,alt,N):
        """
        Brake Torque
        
        Inputs:
            N   - RPM
            alt - altitude
        """
        BMEP    = self.BMEP(alt,N)
        vd      = self.Vd
        nr      = self.NumRev
                
        return self.GearRatio*BMEP*vd/(2.*math.pi*nr)
    
#===============================================================================   
    def Pb(self,alt,N):
        """
        Brake Power
        
        Inputs:
            N   - RPM
            alt - altitude
        """
        Tb      = self.Tb(alt,N)
        
        return N*Tb * self.GearEff

#===============================================================================
    def Mdotfr(self):
        """
        Rated Fuel Flow at min specific fuel consumption
        """
        Pb      = self.Pb
        sfcmt   = self.SFCmt
        Nr      = self.Nr()
        rnsfc   = self.Rnsfc 
        
        return sfcmt*Pb(0*FT,Nr*rnsfc)
 
#===============================================================================
    def Mdotf(self,alt,N):
        """
        Fuel Flow
        """
        Pi      = self.Pi
        Mdotfr  = self.Mdotfr()
        Nr      = self.Nr()
        rnsfc   = self.Rnsfc 
        
        return Mdotfr*Pi(alt,N)/Pi(0*FT,rnsfc*Nr)
    
#===============================================================================
    def Mdota(self,alt,N):
        """
        Air Flow (a=air)
        
        Inputs:
            N   - RPM
            alt - altitude
        """
        a_f     = self.A_F
        Mdotf   = self.Mdotf(alt,N)
 
        
        return Mdotf*a_f
    
#===============================================================================
    def FuelFlow(self,alt,N):
        """
        Volume Fuel Flow
        
        Inputs:
            N   - RPM
            alt - altitude
        """
        Mdotf   = self.Mdotf(alt,N)
 
        return Mdotf/(6*LBM/GAL)
    
#===============================================================================
    def SFC(self,alt,N):
        """
        Specific Fuel Consumption
        
        Inputs:
            N   - RPM
            alt - altitude
        """
        Mdotf   = self.Mdotf(alt,N)
        Pb      = self.Pb(alt, N)
 
        return Mdotf/Pb
  
#===============================================================================
    def PlotTestData(self, fig = 1, RPMmin = None, RPMmax = None):
        """
        Plots Torque and Power along with test data
        
        Inputs:
            fig - Figure number
        """
        
        if self.TestData is None:
            raise ACEngineError("No test data has been given")
  
        TestData = npy.array(self.TestData)

        Nmt = self.Nmt()
        Tm  = self.Tm() / self.TorqueUnit
        Pmt = self.Pb(0*FT,Nmt) /self.PowerUnit
        Nmt /= RPM

        Pr = self.Pr() / self.PowerUnit
        Tr = self.Tr() / self.TorqueUnit
        Nr = self.Nr() / RPM

        if RPMmax is None:
            Nmax = max( max(TestData[:,0]) / RPM, Nmt, Nr)
        else:
            Nmax = RPMmax / RPM
            
        Nmin = 0
        if RPMmin is not None:
            Nmin = RPMmin /RPM
             
        N  = npy.linspace(Nmin,Nmax,20)*RPM
        Tb = self.Tb(0*FT, N)
        Ti = self.Ti(0*FT, N)
        Pb = self.Pb(0*FT, N)
        
        N  = N / RPM
        Tb = Tb / self.TorqueUnit
        Pb = Pb / self.PowerUnit
                
        Ntest = self.ToUnumList(TestData[:,0],RPM)
        Ttest = self.ToUnumList(TestData[:,1],IN*OZF)
  
        Ptest = (Ntest*Ttest) / self.PowerUnit
        Ntest = Ntest / RPM
        Ttest = Ttest / self.TorqueUnit
        
        pyl.figure(fig)
        ax = pyl.subplot(121)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))        
        pyl.plot(N,Tb)
        pyl.plot(Ntest,Ttest,'o')
        pyl.plot(Nmt,Tm,'o')
        pyl.plot(Nr,Tr,'o')
        pyl.xlabel('RPM'); pyl.ylabel('Torque (' + self.TorqueUnitName + ')')
        pyl.legend(['Theoretical', 'Measured', 'Max Torque', 'Max Power'], loc = 'best')
        pyl.ylim( ymin = 0 )
        
        ax = pyl.subplot(122)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))        
        pyl.plot(N,Pb)
        pyl.plot(Ntest,Ptest,'o')
        pyl.plot(Nmt,Pmt,'o')
        pyl.plot(Nr,Pr,'o')
        pyl.xlabel('RPM'); pyl.ylabel('Power (' + self.PowerUnitName + ')') 
        pyl.legend(['Theoretical', 'Measured', "Max Torque", 'Max : ' + AsUnit(self.Pr(), self.PowerUnitName, "%3.1f") ], loc = 'best')
        pyl.ylim( ymin = 0 )
        
        pyl.annotate("Engine Torque and Power Test Matching", xy=(.025, .975),
                xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='top',
                fontsize=20)

################################################################################
if __name__ == '__main__':
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
    
    Tb = Engine.Tb(0*FT, 15997*RPM)
    Pb = Engine.Pb(0*FT, 15997*RPM)
    
    N  = npy.linspace(6000,18000,13)*RPM
    h  = npy.linspace(0,30000,7)*FT
    Nr = Engine.Nr()
    
    Tb    = Engine.Tb(0*FT, N)
    Pb    = Engine.Pb(0*FT, N)
    mdotf = Engine.Mdotf(h, Nr)
    SFC   = Engine.SFC(0*FT, N)
    
    N     = N / RPM
    h     = h / FT
    mdotf = mdotf / (LBM/HR)
    SFC   = SFC / PSFC
    Tb    = Tb / (FT*LBF)
    Pb    = Pb / HP
    
    pyl.figure(1)
    pyl.subplot(121)
    pyl.plot(N,Tb)
    pyl.xlabel('RPM'); pyl.ylabel('Tb (FT*LBF)')
    
    pyl.subplot(122)
    pyl.plot(N,Pb)
    pyl.xlabel('RPM'); pyl.ylabel('Pb (HP)')
    
    pyl.figure(2)
    pyl.subplot(121)
    pyl.plot(h,mdotf)
    pyl.xlabel('Altitude (FT)'); pyl.ylabel('Fuel Flow Rate (LBM/HR)')
    
    pyl.subplot(122)
    pyl.plot(N,SFC)
    pyl.xlabel('RPM'); pyl.ylabel('Specific Fuel Consumption (PSFC)')
    
    pyl.show()
    
#    print 'Bore   = ',Engine.B()
#    print 'Stroke = ',Engine.L()
#    print ''
#    Nr = Engine.Nr()
#    print 'Rated crankshaft rotational speed (Nr)   = ',Nr
#    print 'Torque at Nr                             = ',Engine.Tr()
#    print 'Brake power at Nr                        = ',Engine.Pr()
#    print 'RPM at max torque                        = ',Engine.Nmt()
#    print 'Max Torque                               = ',Engine.Tm()
#    print 'Max brake mean effective pressure        = ',Engine.BMEPMT()
#    print 'Quadratic fit of BMEP at sea level at Nr = ',Engine.BMEPsl(Nr)
#    print 'Max indicated mean effective pressure    = ',Engine.IMEPMT()
#    print 'Assumed total losses                     = ',Engine.MEPtlsl(Nr)
#    print 'Total friction losses                    = ',Engine.MEPtf(Nr)
#    print 'Total pumping losses                     = ',Engine.MEPtp(0*FT, Nr)
#    print 'Indicated (gas) mean effective pressure  = ',Engine.IMEP(0*FT)
#    print 'Brake mean effective pressure            = ',Engine.BMEP(0*FT, Nr)
#    print ''
#    print 'Indicated torque = ',Engine.Ti(0*FT)
#    print 'Indicated power  = ',Engine.Pi(0*FT, Nr)
#    print 'Brake torque     = ',Engine.Tb(0*FT, Nr)
#    print 'Brake power      = ',Engine.Pb(0*FT, Nr)
#    print ''
#    print 'Rated fuel flow at minimum sfc = ',Engine.Mdotfr()
#    print 'Fuel flow (mass/hr)            = ',Engine.Mdotf(0*FT, Nr)
#    print 'Air flow  (mass/hr)            = ',Engine.Mdota(0*FT, Nr)
#    print 'Fuel flow (volume/hr)          = ',Engine.FuelFlow(0*FT, Nr)
#    print 'Specific fuel consumption      = ',Engine.SFC(0*FT, Nr)
#    print ''
    print 'Test Completed'