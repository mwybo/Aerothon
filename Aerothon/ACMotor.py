"""
This is the class for analyzing an electric motor
Used for Propulsion analysis
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACBase, Current, Resistance, ElectricPotential, ElectricCapacity, RPMPerVolt, Unitless, ACUnitCheck, Force, Power, Frequency
from ACMass import ACMassCyl, ACMassBox
from scalar.units import RPM, OHM, A, V, IN, OZF, HP, W, HR, mAh, GRAM, gacc
from scalar.units import AsUnit
import cmath as math
import numpy as npy
import pylab as pyl
import AeroUtil as AUtil

################################################################################
class ACMotorError(Exception):
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg
    
################################################################################
class ACSpeedController(ACMassBox):
    """
    Speed Controller class
        
    Attributes:
        Imax  - Maximum Current
    """
    def __init__(self):
        super(ACSpeedController,self).__init__()
        
        self.name = 'SpeedController'
        
        UnitList = self.UnitList
        self.NoneList['Imax']  = None  ; UnitList['Imax'] = Current

################################################################################
class ACBattery(ACMassBox):
    """
    Battery class
        
    Attributes:
        --Battery Properties
        Voltage    - Battery Voltage
        Capacity   - Battery Capacity (Amp-Hour)
        C_Rating   - Discharge rating of the battery in 'C'
        Cells      - Number of Cells
        Cells_V    - Voltage of a single cell (default 3.7V)
        Imax       - Maximum Current Draw
        
    """
    def __init__(self):
        super(ACBattery,self).__init__()
        
        self.name = 'Battery'
        
        UnitList = self.UnitList
        self.NoneList['Voltage']   = None  ; UnitList['Voltage']   = ElectricPotential
        self.NoneList['Capacity']  = None  ; UnitList['Capacity']  = ElectricCapacity
        self.NoneList['C_Rating']  = None  ; UnitList['C_Rating']  = Unitless
        self.NoneList['Cells']     = None  ; UnitList['Cells']     = Unitless
        self.__dict__['Cells_V']   = 3.7*V ; UnitList['Cells_V']   = ElectricPotential
        self.NoneList['Imax']      = None  ; UnitList['Imax']      = Current
       
#===============================================================================   
    def Power(self):
        """
        Computes the power of the battery
        """
        return self.Capacity*self.Voltage
       
#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACBattery, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'Voltage':
            return self._CalcVoltage()
        elif key == 'Imax':
            return self._CalcImax()
        elif key == 'C_Rating':
            return self._Calc_C_Rating()
        #
        # If no calculation exist return None
        #
        return None
    
#===============================================================================   
    def _CalcVoltage(self):
        """
        Computes the voltage based on the cell count
        """
        return self.Cells*self.Cells_V
           
#===============================================================================   
    def _CalcImax(self):
        """
        Computes the maximum allowed current
        """
        return self.Capacity*self.C_Rating/HR

#===============================================================================   
    def _Calc_C_Rating(self):
        """
        Computes the maximum allowed current
        """
        return self.Imax/self.Capacity
       
       
################################################################################
class ACMotor(ACMassCyl):
    """
    Motor class
    
    Much of the theory is based on this webpage
    http://www.innovatia.com/Design_Center/Power_Electronics.htm
    
    Another useful guide can be found at
    http://www.adamone.rchomepage.com/guide5.htm
    
    Attributes:
        --Battery Properties
        Vb        - Battery Voltage
        Cb        - Battery Capacity (Amp-Hour)
        
        --Motor Properties
        Ri        - Coil resistance (Ohms)
        Io        - Idle current (Ampere)
        Kv        - Voltage Constant (Ratio of the RPM to the Voltage at the motor's terminals minus 
                    the Voltage loss inside the motor due to the coil's resistance Ri)
        Kt        - Torque Constant
        aKt       - Constant to relate Kt to Kv ( Kt = aKt/Kv)
        
        Wmax      - Maximum continuous watts
        Imax      - Maximum continuous current
        Vmax      - Maximum voltage
        Nmax      - Maximum RPM for the voltage Vmax

        --Losses
        xRm       - Unitless multiplier such that Rm = xRm * Ri where Rm is the shunt resistance for magnetic loss
        Pz        - Ratio of drive power loss to maximum input power
        
        --Units
        ThrustUnit     - Unit for thrust (default ozf)
        ThrustUnitName - Name of unit for thrust (default ozf)
        PowerUnit      - Unit for power (default watt)
        PowerUnitName  - Name of unit for power (default watt)

        --Gear box
        GearRatio - Gear ratio for a gear box (default 1)
        GearEff   - Gear Box efficiency (default 1.0)
               
        TestData  - An array of test data
                    #            RPM,        Torque           Current     Voltage
                    TestData = [(7000  *RPM, 107.35 *IN*OZF,   12*A,      6.4*V),
                                (8500  *RPM, 104.24 *IN*OZF,   11*A,      6.4*V),
                                (9500  *RPM, 101.13 *IN*OZF,   10*A,      6.4*V)]
    """
    
#===============================================================================    
    def __init__(self):
        super(ACMotor,self).__init__()
        
        self.name = 'Motor'
        
        UnitList = self.UnitList
        #self.NoneList['Vb']        = None                  ; UnitList['Vb']        = ElectricPotential
        #self.NoneList['Cb']        = None                  ; UnitList['Cb']        = ElectricCapacity
        self.NoneList['Ri']        = None                  ; UnitList['Ri']        = Resistance
        self.NoneList['Io']        = None                  ; UnitList['Io']        = Current
        self.NoneList['Kv']        = None                  ; UnitList['Kv']        = RPMPerVolt
        self.NoneList['Kt']        = None                  ; UnitList['Kt']        = ACUnitCheck(IN*OZF/A, 'Electric Motor Torque Slope')
        self.__dict__['Kta']       = 1352*IN*OZF*RPM/(A*V) ; UnitList['Kta']       = ACUnitCheck(IN*OZF*RPM/(A*V), 'Electric Motor Torque Conversion')
        
        self.NoneList['Wmax']      = None                  ; UnitList['Wmax']      = Power
        self.NoneList['Imax']      = None                  ; UnitList['Imax']      = Current
        self.NoneList['Vmax']      = None                  ; UnitList['Vmax']      = ElectricPotential
        self.NoneList['Nmax']      = None                  ; UnitList['Nmax']      = Frequency

        self.NoneList['xRm']       = None                  ; UnitList['xRm']       = Unitless
        self.NoneList['Pz']        = None                  ; UnitList['Pz']        = Unitless

        self.__dict__['Battery']         = None
        self.__dict__['SpeedController'] = None

        self.__dict__['ThrustUnit'] = OZF                  ; UnitList['ThrustUnit'] = Force
        self.__dict__['ThrustUnitName'] = 'ozf'

        self.__dict__['PowerUnit'] = W                     ; UnitList['PowerUnit'] = Power
        self.__dict__['PowerUnitName'] = 'watt'          
 
        self.__dict__['GearRatio'] = 1                    ; UnitList['GearRatio'] = Unitless
        self.__dict__['GearEff']   = 1                    ; UnitList['GearEff']   = Unitless
                        
        self.__dict__['TestData']  = None       

#===============================================================================
#    def AddToWeightTable(self, PartName, WeightTable):
#        """
#        Adds the lift surface and its parts to the Weight Table
#        
#        Input:
#            PartName    - The name give to this part
#            WeightTable - The table to add self to 
#        """
#        if self.dirty: self.Refresh()
#        
#        WeightTable[PartName] = self
#                
#        self.Battery.AddToWeightTable(Battery.name, WeightTable[PartName])
#        self.SpeedController.AddToWeightTable(SpeedController.name, WeightTable[PartName])
        
#===============================================================================
#    def _CalcWeight(self):
#        """
#        Calculates the total weight
#        """
#        if self.dirty: self.Refresh()
#        
#        return super(ACMotor,self).Weight + Battery.Weight + SpeedController.Weight
                
#===============================================================================
    def Length(self):
        """
        Computes the length of the motor
        """
        return self.LenDi[0]
    
#===============================================================================
    def PropX(self):
        """
        Computes the position of the propeller on the motor
        """
        X     = self.GetX().copy()
        X[0] -= self.LenDi[0]*self.Axis[0]*self.Xat[0]
        
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
        super(ACMotor, self).Refresh()

        self.param.refreshing = True
                
        self.param.refreshing = False

#===============================================================================   
    def _CalcKt(self):
        """
        Computes the torque slope
        """
        return self.Kta/self.Kv

#===============================================================================   
    def _CalcWmax(self):
        """
        Computes the maximum allowed power
        """
        return self.Vmax*self.Imax

#===============================================================================   
    def _CalcImax(self):
        """
        Computes the maximum allowed current
        """
        return self.Wmax/self.Vmax

#===============================================================================   
    def _CalcVmax(self):
        """
        Computes the maximum allowed voltage
        """
        return self.Wmax/self.Imax
    
#===============================================================================
    def _CalcNmax(self):
        """ 
        Calculates Max RPM
        """
        Kv = self.Kv
        Vb = self.Battery.Voltage
        Ri = self.Ri
        Io = self.Io
                
        N  = Kv*(Vb - Io*Ri)
        
        return N / self.GearRatio

#===============================================================================
    def _CalcRi(self):
        """ 
        Calculates the motor resistance based on Nmax and Vb
        """
        Kv = self.Kv
        Vmax = self.Vmax
        Nmax = self.Nmax
        Io = self.Io
                
        Ri = (Vmax - Nmax/Kv)/Io
        
        return Ri

#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
#        if key == 'Weight':
#            return self._CalcWeight()

        ans = super(ACMotor, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'Kt':
            return self._CalcKt()
        elif key == 'Ri':
            return self._CalcRi()
        elif key == 'Wmax':
            return self._CalcWmax()
        elif key == 'Imax':
            return self._CalcImax()
        elif key == 'Vmax':
            return self._CalcVmax()
        elif key == 'Nmax':
            return self._CalcNmax()
        #
        # If no calculation exist return None
        #
        return None
    
#===============================================================================
    def NRange(self, npt = 30):
        """ 
        Range of RPM
        """
        Nmax  = self.Nmax/RPM
        return npy.linspace(0, Nmax, npt)*RPM

#===============================================================================
    def Tmax(self):
        """ 
        Max torque
        """
        return self.Tb(0, 0*RPM)

#===============================================================================
    def Pmax(self):
        """ 
        Max power
        """
        Kv = self.Kv
        Vb = self.Battery.Voltage
        Ri = self.Ri
        Io = self.Io
        
        N = (Vb/Ri - Io)*(Kv*Ri)/2
        
        return self.Pb(0,N)

#===============================================================================   
    def I_Effmax(self):
        """
        Current for maximum efficiency
        """        
        return self.Ib( self.N_Effmax() )
    
#===============================================================================   
    def N_Effmax(self):
        """
        RPM for maximum efficiency
        
        Pout = Tb*N - Pm = Kt*(Iin - Io)*N - Pm = 
             = Kt*((Vb - Kt*N)/Ri - Io)*N - (Kt*N)**2/Rm
        Pin = Vb*Iin = Vb/Ri*(Vb - Kt*N)
        
        eta = (N/Nmax)*(1 - (Ri/Rm)*(N/(Nmax + N)))
        Solve( deta/dN = 0 )
        """
        xRmp1 = self.xRm + 1
        Pz    = self.Pz
        Vb    = self.Battery.Voltage
        Io    = self.Io
        Rm    = self.Ri * self.xRm
        Kt    = self.Kt
                
        return ((Pz+1)*Vb*xRmp1 - V*math.sqrt( ((Pz+1)*Vb*xRmp1*(Io*Rm + Vb*(Pz*xRmp1 + 1))) / V**2).real)/(Kt*xRmp1) / self.GearRatio


#===============================================================================   
    def Effmax(self):
        """
        Maximum efficiency
        """
        return self.Efficiency(self.N_Effmax())

#===============================================================================   
    def MagneticLoss(self, N=None, Ib=None):
        """
        Computes the magnetic loss
        """
        if N is None and Ib is None:
            raise ACMotorError("Must specify N or Ib.")

        Ri = self.Ri
        Rm = self.xRm * Ri
        Kt = self.Kt
        if N is None:
            N = self.Kv*(self.Battery.Voltage - Ib*Ri)
            
        return (Kt * N*self.GearRatio)**2/Rm
             

#===============================================================================   
    def Ib(self,N):
        """
        Brake Current (Current used by motor to achieve the desired RPM)
        
        Inputs:
            N   - RPM
        """
        Kv = self.Kv
        Vb = self.Battery.Voltage
        Ri = self.Ri
        
        Vloss = Vb - N*self.GearRatio/Kv
        
        return Vloss/Ri
    
#===============================================================================   
    def Tb(self,alt,N=None,Ib=None):
        """
        Brake Torque
        
        Inputs:
            N   - RPM
            alt - Altitude (not used)
        """
        
        Kt = self.Kt
        if N is not None:
            Iin = self.Ib(N)
        elif Ib is not None:
            Iin = Ib
        else:
            raise ACMotorError("Must specify N or Ib.")
                
        return self.GearRatio*Kt*(Iin - self.Io)
    
#===============================================================================   
    def Pb(self,alt,N=None,Ib=None):
        """
        Brake Power
        
        Inputs:
            N   - RPM
            alt - Altitude (not used)
        """
        if N is None and Ib is None:
            raise ACMotorError("Must specify N or Ib.")

        Pm = self.MagneticLoss(N,Ib)
        
        if N is not None:
            Tb = self.Tb(alt,N,Ib)
            return (N*self.GearRatio*Tb -  Pm) * self.GearEff / self.GearRatio
        elif Ib is not None:
            Vb = self.Battery.Voltage
            Ri = self.Ri
            Io = self.Io
            return ((Vb - Ib*Ri)*(Ib - Io) -  Pm) * self.GearEff

#===============================================================================   
    def Pin(self,N=None,Ib=None):
        """
        Input Power
        
        Inputs:
            N   - RPM
            alt - Altitude (not used)
        """
        if N is None and Ib is None:
            raise ACMotorError("Must specify N or Ib.")

        if Ib is None:
            Ib = self.Ib(N)
        
        Vb = self.Battery.Voltage
        
        return Vb*Ib + self.Pz*Vb**2/self.Ri
        
#===============================================================================   
    def Efficiency(self, N=None, Ib=None):
        """
        Power efficiency
        
        Inputs:
            N  - RPM
            Ib - Input Current
        """
        if N is None and Ib is None:
            raise ACMotorError("Must specify N or Ib.")
        
        Pin = self.Pin(N,Ib)
        Pout = self.Pb(0,N,Ib)
        
        return Pout/Pin*100
 
 #===============================================================================   
    def Duration(self, N = None, Ib = None):
        """
        The duration at max throttle
        """
        if N is None and Ib is None:
            raise ACMotorError("Must specify N or Ib.")

        if N is not None:
            Ib = self.Ib(N)
        
        Cb = self.Battery.Capacity
        
        return Cb/Ib

#===============================================================================   
    def Vmotor(self, Ib):
        """
        Voltage that the motor sees as a result of losses
        """
        return self.Battery.Voltage - Ib*self.Ri

#===============================================================================
    def PlotProperties(self, fig = 1):
        """
        Plots Torque, Power, Current and Efficiency
        
        Inputs:
            fig - Figure number
        """
                
        N   = self.NRange()
        Tb  = self.Tb(0, N)
        Pb  = self.Pb(0, N)
        Ib  = self.Ib(N)
        eta = self.Efficiency(N = N)
        
        Em = self.Effmax()
        NEm = self.N_Effmax()
        TbEm = self.Tb(0, NEm)
        PbEm = self.Pb(0, NEm)
        IEm = self.I_Effmax()

#        N   /= self.GearRatio
#        NEm /= self.GearRatio
        
        N  /= RPM
        Tb /= (IN*OZF)
        Pb /= self.PowerUnit
        Ib /= A
        
        NEm  /= RPM
        TbEm /= (IN*OZF)
        PbEm /= self.PowerUnit
        IEm  /= A
                 
        pyl.figure(fig)
        pyl.subplot(221)
        pyl.plot(NEm,TbEm,'or')
        pyl.plot(N,Tb)
        pyl.xlabel('RPM'); pyl.ylabel('Torque (in*ozf)')
        pyl.legend(['Max Eff.'], loc = 'best')
        ymin, ymax = pyl.ylim()
        pyl.ylim(ymin = 0, ymax = max(ymax, max(Tb)*1.1))
    
        pyl.subplot(222)
        pyl.plot(NEm,PbEm,'or')
        pyl.plot(N,Pb)
        pyl.xlabel('RPM'); pyl.ylabel('Power (' + self.PowerUnitName + ')') 
        pyl.legend(['Max Eff.'], loc = 'best')
        ymin, ymax = pyl.ylim()
        pyl.ylim(ymin = 0, ymax = max(ymax, max(Pb)*1.1))

        pyl.subplot(223)
        pyl.plot(NEm,IEm,'or')
        pyl.plot(N,Ib)
        pyl.xlabel('RPM'); pyl.ylabel('Current (A)') 
        pyl.legend(['Max Eff.'], loc = 'best')
        ymin, ymax = pyl.ylim()
        pyl.ylim(ymin = 0, ymax = max(ymax, max(Ib)*1.1))

        pyl.subplot(224)
        pyl.plot(NEm,Em,'or')
        pyl.plot(N,eta)
        pyl.xlabel('RPM'); pyl.ylabel('Efficiency') 
        pyl.legend(['Max Eff.'], loc = 'best')
        pyl.ylim( (0, 100) )
        
#        pyl.annotate("Motor Torque and Power Test Matching", xy=(.025, .975),
#                xycoords='figure fraction',
#                horizontalalignment='left', verticalalignment='top',
#                fontsize=20)

#===============================================================================
    def PlotTestData(self, fig = 1, TorqueLimit = None, PowerLimit = None, CurrentLimit = None):
        """
        Plots Torque and Power along with test data
        
        Inputs:
            fig - Figure number
        """
        
        if self.TestData is None:
            raise ACMotorrError("No test data has been given")
  
        TestData = npy.array(self.TestData)

        Nmax = max( max(TestData[:,0]) / RPM, self.Nmax / RPM )
        
        N    = npy.linspace(0,Nmax,20)*RPM
        Tb   = self.Tb(0, N)
        Pb   = self.Pb(0, N)
        Pin  = self.Pin(N)
        Ib   = self.Ib(N)
        eta  = self.Efficiency(N = N)
                
        Em = self.Effmax()
        NEm = self.N_Effmax()
        TbEm = self.Tb(0, NEm)
        PbEm = self.Pb(0, NEm)
        IEm = self.I_Effmax()
        
#        N   /= self.GearRatio
#        NEm /= self.GearRatio
        
        N   /= RPM
        Tb  /= (IN*OZF)
        Pb  /= self.PowerUnit
        Pin /= self.PowerUnit
        Ib  /= A
        
        NEm  /= RPM
        TbEm /= (IN*OZF)
        PbEm /= self.PowerUnit
        IEm  /= A
         
        Ntest = self.ToUnumList(TestData[:,0],RPM)
        Ttest = self.ToUnumList(TestData[:,1],IN*OZF)
        Itest = self.ToUnumList(TestData[:,2],A)
        Vtest = self.ToUnumList(TestData[:,3],V)
  
        Ptest = (Ntest*Ttest)/ self.PowerUnit
        Pintest = (Vtest*Itest)/ self.PowerUnit
        Ntest = Ntest / RPM
        Ttest = Ttest / (IN*OZF)
        Itest = Itest / A
        Vtest = Vtest / V
        etatest = Ptest/Pintest*100
        
        ylabel_offset = -0.1
        
        pyl.figure(fig)
        ax = pyl.subplot(221)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))

        junk, Mechanical, MaxEffPnt = pyl.plot( N,Tb, 'b'
                                              , Ntest,Ttest,'ob'
                                              , NEm,TbEm,'or')
        pyl.xlabel('RPM'); pyl.ylabel('Torque (in*ozf)')
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
        pyl.grid()
        #pyl.legend(['Model', 'Measured', 'Max Eff.'], loc = 'best')
        pyl.ylim(ymin = 0)
        if TorqueLimit is not None:
            pyl.ylim(ymax = TorqueLimit)
    
        ax = pyl.subplot(222)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
        MechModel, ElecModel, Electric, MaxLim, junk, junk \
            = pyl.plot( N,Pb, 'b'
                      , N,Pin,'g'
                      , Ntest,Pintest,'og'
                      , [min(N), max(N)],[self.Wmax/self.PowerUnit, self.Wmax/self.PowerUnit],'r'
                      , Ntest,Ptest,'ob'  
                      , NEm,PbEm,'or')
        
        pyl.figlegend([MechModel, ElecModel, Mechanical, Electric, MaxEffPnt, MaxLim], ['Mech. Model', 'Elec. Model', 'Mech. Test Data', 'Elec. Test Data', 'Max Eff.', 'Limit'], loc = 'upper center', numpoints = 1, ncol = 3, frameon = False, labelspacing = 0.1)

        pyl.xlabel('RPM'); pyl.ylabel('Power (' + self.PowerUnitName + ')')
        pyl.gca().yaxis.set_label_coords(-0.12,0.5)
        pyl.grid()
        #pyl.legend(['Break', 'Input'], loc = 'best')
        pyl.ylim(ymin = 0)
        if PowerLimit is None:
            PowerLimit = self.Wmax/self.PowerUnit*1.1
        pyl.ylim(ymax = PowerLimit)

        ax = pyl.subplot(223)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
        pyl.plot([min(N), max(N)],[self.Imax/A, self.Imax/A],'r')
        pyl.plot(N,Ib,'b')
        pyl.plot(Ntest,Itest,'og')
        pyl.plot(NEm,IEm,'or')
        pyl.xlabel('RPM'); pyl.ylabel('Current (A)') 
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
        pyl.grid()
        #pyl.legend(['Max'], loc = 'best')
        pyl.ylim(ymin = 0)
        if CurrentLimit is not None:
            pyl.ylim(ymax = CurrentLimit)

        ax = pyl.subplot(224)
        ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
        pyl.plot(N,eta)
        pyl.plot(Ntest,etatest,'og')
        pyl.plot(NEm,Em,'or')
        pyl.xlabel('RPM'); pyl.ylabel('Efficiency') 
        pyl.ylim( (0, 100) )
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
        pyl.grid()
                
#        pyl.annotate("Calibration of " + self.name + " Model with Experimental Data", xy=(.025, .975),
#                xycoords='figure fraction',
#                horizontalalignment='left', verticalalignment='top',
#                fontsize=20)

################################################################################
if __name__ == '__main__':
    Battery = ACBattery()
    Battery.Cells = 3
    Battery.Capacity = 450*mAh
    Battery.C_Rating = 65
    Battery.Weight = 41.5*GRAM*gacc
    
    Motor  = ACMotor()
    Motor.Battery = Battery
    Motor.Ri  = 0.35*OHM
    Motor.Io  = 0.46*A
    Motor.Kv  = 1900*RPM/V
    Motor.xRm = 1
    Motor.Pz = 0
    Nmax      = Motor.Nmax/RPM

    print AsUnit( Motor.Pmax(), 'W' )
    print AsUnit( Motor.I_Effmax(), 'A' )
    print AsUnit( Motor.N_Effmax(), 'rpm' )
    print 'Max Eff : ', Motor.Effmax()
    
    N  = npy.linspace(0, Nmax, 20)*RPM
    
    Tb    = Motor.Tb(0, N)
    Pb    = Motor.Pb(0, N)
    eta   = Motor.Efficiency(N = N)
    
    N     = N / RPM
    Tb    = Tb / (IN*OZF)
    Pb    = Pb / Motor.PowerUnit
    
    pyl.figure(1)
    pyl.subplot(131)
    pyl.plot(N,Tb)
    pyl.xlabel('RPM'); pyl.ylabel('Tb (IN*OZF)')
    
    pyl.subplot(132)
    pyl.plot(N,Pb)
    pyl.xlabel('RPM'); pyl.ylabel('Pb (' + Motor.PowerUnitName + ')')
    pyl.ylim( (0, max(Pb)*1.1 ))

    pyl.subplot(133)
    pyl.plot(N,eta)
    pyl.xlabel('RPM'); pyl.ylabel('Efficiency')
    pyl.ylim( (0, max(eta)*1.1 ))
        
    pyl.show()
    
