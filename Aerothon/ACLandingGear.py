"""
A landing-gear class
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACStructural, Length, Angle, g
from ACMass import ACMassBox, ACMassCyl
from scalar.units import IN, ARCDEG, RAD, LBF, LBM
from math import sin, cos
import numpy as npy
import pylab as pyl

################################################################################
class ACLandingGear(ACStructural):
    """
    A class for landing-gear
    
    TODO: Document attributes for the main gear
    """
#===============================================================================
    def __init__(self):
        ACStructural.__init__(self)
        UnitList = self.UnitList
        NoneList = self.NoneList

        # Strut parameters
        self.NoneList['StrutL']     = None       ; UnitList['StrutL']     = Length
        self.__dict__['StrutW']     = 1 * IN     ; UnitList['StrutW']     = Length
        self.__dict__['StrutH']     = 1 * IN     ; UnitList['StrutH']     = Length
        self.__dict__['Theta']      = 0 * ARCDEG ; UnitList['Theta']      = Angle
        self.NoneList['GearHeight'] = None       ; UnitList['GearHeight'] = Length

        # Wheel parameters
        self.__dict__['WheelDiam']      = 1 * IN ; UnitList['WheelDiam']      = Length
        self.__dict__['WheelThickness'] = 0.1 * IN   ; UnitList['WheelThickness'] = Length

        

        # Internal variables
        self.__dict__['Strut'] = ACMassBox()
        self.Strut.Weight      = 1e-5 * LBF
        self.Strut.name        = 'Strut'
        self.__dict__['Wheel'] = ACMassCyl()
        self.Wheel.Weight      = 1e-5 * LBF
        self.Wheel.name        = 'Wheel'

#===============================================================================
    def Refresh(self):
        """
        Update drawing array

        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        ACStructural.Refresh(self)
        self.param.refreshing = True
        
        self._CheckConsistent()

        Strut = self.Strut
        Wheel = self.Wheel

        Theta  = self.Theta / RAD

        #
        # Calculate strut parameteres
        #
        StrutAxis = npy.array([0, sin(Theta), -cos(Theta)])
        StrutVect = StrutAxis * self.StrutL

        #
        # Set the strut parameters
        # Position the height such that the wheel sits on the ground
        #
        Strut.Axis = StrutAxis
        Strut.LWH  = [self.StrutL, self.StrutW, self.StrutH]
        Strut.X    = self.X

        #
        # Position the wheel
        #
        WheelHub    = npy.array([0 * IN, self.StrutH, 0 * IN])
        Wheel.X     = npy.array(Strut.X) + StrutVect + WheelHub
        Wheel.Axis  = [0, 1, 0]
        Wheel.LenDi = [self.WheelThickness, self.WheelDiam]

        #
        # Make the wheel and strut symmetric if desired
        #
        Strut.Symmetric = self.Symmetric
        Wheel.Symmetric = self.Symmetric

        self.param.refreshing = False

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
        ACStructural.Draw(self, fig, top, side, front)
        self.Wheel.Draw(fig, top, side, front)
        self.Strut.Draw(fig, top, side, front)

#===============================================================================
    def _CalcGearHeight(self):
        " Computes the height of the landing gear including the wheel "
        Theta = self.Theta / RAD
        return cos(Theta) * self.StrutL + self.WheelDiam / 2

#===============================================================================
    def _CalcStrutL(self):
        " Computes the length of the landing gear strut "
        Theta = self.Theta / RAD
        return (self.GearHeight - self.WheelDiam / 2) / cos(Theta)

#===============================================================================
    def _CalcWeight(self):
        " Computes the weight of the landing gear "
        fac = 2. if self.Symmetric else 1
        return fac*(self.Strut.Weight + self.Wheel.Weight)

#===============================================================================
    def WheelBase(self):
        " Computes the location where the wheels rest on the ground "
        if self.dirty:
            self.Refresh()

        Wheel = self.Wheel

        Base = Wheel.GetX() - npy.array([0,0, 0.5 * self.WheelDiam / IN])*IN

        return Base
    
#===============================================================================
    def CG(self):
        " Computes the CG of the landing gear "
        if self.dirty:
            self.Refresh()

        Strut = self.Strut
        Wheel = self.Wheel

        CG = ((Strut.CG() * Strut.Weight + Wheel.CG() * Wheel.Weight) / (Strut.Weight + Wheel.Weight))

        if self.Symmetric:
            CG[1] = 0 * IN

        return CG
    
#===============================================================================
    def MOI(self):
        "Computes the MOI of the landing gear"
        if self.dirty:
            self.Refresh()
            
        Strut = self.Strut
        Wheel = self.Wheel
        
        CG = self.CG()
        
        TempMOI = npy.array([0,0,0])*LBM*IN**2
        
        def AddMOI(PCG, PM, PMOI):
            dCG = PCG - CG
            TempMOI[0] += PMOI[0] + PM*(dCG[1]**2 + dCG[2]**2) #Ixx
            TempMOI[1] += PMOI[1] + PM*(dCG[0]**2 + dCG[2]**2) #Iyy
            TempMOI[2] += PMOI[2] + PM*(dCG[1]**2 + dCG[0]**2) #Izz      
        
        # Add strut
        StrutCG = Strut.CG()
        StrutMOI = Strut.MOI()
        StrutM = Strut.Weight / g
        
        AddMOI(StrutCG, StrutM, StrutMOI)
     
        # Add wheel
        WheelCG = Wheel.CG()
        WheelMOI = Wheel.MOI()
        WheelM = Wheel.Weight / g
        
        AddMOI(WheelCG, WheelM, WheelMOI)
        
        return TempMOI
            
#===============================================================================
    def CD(self, Wing, Alpha):
        """
        Computes the drag of the landing gear
        
        Inputs:
            Alpha - 2-D angle of attack of the wing
        """
        if self.dirty:
            self.Refresh()
            
        WheelDiam = self.WheelDiam
        WheelThickness = self.WheelThickness
        StrutL = self.StrutL
        StrutT = self.StrutH
        WingLower = Wing.Lower(0 * IN)
        CLWing = Wing.CL(Alpha) 
        S = Wing.S
        b = Wing.b
        DqsTire = .3 # Does this value need read in from a table
        DqsStrut = 1.4 # Does this value need read in from a table
        
        STire = WheelDiam * WheelThickness
        cg = S / b      #Geometric Mean Chord
        luc = WingLower - (WheelDiam / 2)    #Distance from Tire Axle To Bottom of Wing
        Fuc = (1 - (0.04 * CLWing * cg)/luc)**2 # Thombeek pg550 G67
        
        # Thombeek pg550 G67
        Dqg1 = 1.5 * STire * Fuc
        Dqg2 = 2.7 * (DqsTire * WheelDiam * WheelThickness + DqsStrut * StrutL * StrutT)
        
        Dqg1 = self.MakeList(Dqg1)
        
        Dqg = []
        for dqg in Dqg1:
            if dqg > Dqg2:
                Dqg.append(dqg / (IN**2)) # Use the largest drag value
            else: 
                Dqg.append(Dqg2 / (IN**2))
        
        Dqg = npy.array(Dqg)*IN**2
        
        if self.Symmetric:
            Dqg *= 2
        
        return Dqg[0] if len(Dqg) == 1 else Dqg
            
#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the landing gear and its parts to the Weight Table
        
        Input:
            PartName    - The name give to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WeightTable[PartName] = self
        
        self.Strut.AddToWeightTable('Strut', WeightTable[PartName])
        self.Wheel.AddToWeightTable('Wheel', WeightTable[PartName])
        
#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.
        """
        
        self._CheckEquation(['StrutL','GearHeight'])
            
  
            
#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACLandingGear, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'GearHeight':
            return self._CalcGearHeight()
        elif key == 'StrutL':
            return self._CalcStrutL()

################################################################################
if __name__ == '__main__':

    MainGear = ACLandingGear()
    MainGear.Theta = 45*ARCDEG
    MainGear.GearHeight, MainGear.StrutW, MainGear.StrutH = 5 * IN, 0.2 * IN, 0.1 * IN
    MainGear.Symmetric = False
    MainGear.X = [5 * IN, 2 *IN, MainGear.GearHeight]
    print MainGear.MOI()
    MainGear.Symmetric = True
    print MainGear.MOI()

    NoseGear = ACLandingGear()
    NoseGear.Symmetric = False
    NoseGear.Theta = 0*ARCDEG
    NoseGear.GearHeight, NoseGear.StrutW, NoseGear.StrutH = 5 * IN, 0.2 * IN, 0.1 * IN
    NoseGear.X = [0 * IN, 0 * IN, MainGear.GearHeight]
#    print NoseGear.MOI()
    
    
    print 'Main Gear CG:' , MainGear.CG()

    MainGear.Draw()
    NoseGear.Draw()

    pyl.show()