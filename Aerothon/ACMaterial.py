"""
A class for defining material properties
"""
from __future__ import division # let 5/2 = 2.5 rather than 2

from ACBase import ACBase, UnderSpecifiedError, ForceDensity, Density, Length, AreaDensity, \
                   AreaForceDensity, LinearDensity, LinearForceDensity, Pressure, g
from scalar.scalar import _scalar                   
from scalar.units import LBF, FT, IN, M, SEC, KG, LBM, PSI, OZM
import copy

################################################################################
class MaterialError(Exception):
    """
    Manages errors for material properties
    """

#===============================================================================
    def __init__(self,count,params):

        self.msg = "Material properties over specified.\n"
        self.msg += "Must only specify " + str(count) + " of "
        for i in range(len(params)-1):
            self.msg += "\'" + params[i] + "\', "

        self.msg += "and \'" + params[-1] + "\'"

#===============================================================================
    def __str__(self):
        return self.msg


################################################################################
class ACMaterial(ACBase):
    """
    Handles material properties for a given material
    """

#===============================================================================
    def __init__(self):
        super(ACMaterial, self).__init__()
        UnitList = self.UnitList

        self.NoneList['Thickness']        = None    ;  UnitList['Thickness']        = Length
        self.NoneList['Width']            = None    ;  UnitList['Width']            = Length
        
        self.NoneList['Density']          = None    ;  UnitList['Density']          = Density
        self.NoneList['ForceDensity']     = None    ;  UnitList['ForceDensity']     = ForceDensity

        self.NoneList['AreaDensity']      = None    ;  UnitList['AreaDensity']      = AreaDensity
        self.NoneList['AreaForceDensity'] = None    ;  UnitList['AreaForceDensity'] = AreaForceDensity

        self.NoneList['LinearDensity']      = None  ;  UnitList['LinearDensity']      = LinearDensity
        self.NoneList['LinearForceDensity'] = None  ;  UnitList['LinearForceDensity'] = LinearForceDensity

        #
        # Modulus of elasticity
        #
        self.__dict__['E'] = 0*PSI ;  UnitList['E'] = Pressure
        
        #
        # A dictionary of the variables that can be used to calculations in this class
        # if this is false an desired equation is under specified
        #
        self.param.ValidEqVar = {}

#===============================================================================
    def copy(self):
        """
        Copies the instance of the ACMaterial
        """
        CopyMat = ACMaterial()
        
        CopyMat.NoneList = copy.deepcopy(self.NoneList)
        
        for attr, val in self.UnitList.iteritems():
            try:
                CopyMat.__dict__[attr] = copy.deepcopy(self.__dict__[attr])
            except KeyError:
                pass
        
        return CopyMat
        

#===============================================================================
    def _CalcDensity(self):
        """
        Calculates the density if a force density was specified
        """
        NoneList = self.NoneList
        if not NoneList.has_key('AreaDensity') and not NoneList.has_key('Thickness'):
            return self.AreaDensity / self.Thickness 
        elif not NoneList.has_key('AreaForceDensity') and not NoneList.has_key('Thickness'):
            return self.AreaForceDensity / self.Thickness / g
        else:
            return self.ForceDensity / g

#===============================================================================
    def _CalcAreaDensity(self):
        """
        Calculates the area density
        """
        NoneList = self.NoneList
        if not NoneList.has_key('Density') and not NoneList.has_key('Thickness'):
            return self.Density*self.Thickness
        elif not NoneList.has_key('ForceDensity') and not NoneList.has_key('Thickness'):
            return self.ForceDensity*self.Thickness/g
        else:
            return self.AreaForceDensity/g

#===============================================================================
    def _CalcLinearDensity(self):
        """
        Calculates the linear density
        """
        NoneList = self.NoneList
        if not NoneList.has_key('Thickness') and not NoneList.has_key('Width'):
            
            if not NoneList.has_key('Density'):
                return self.Density*self.Thickness*self.Width
            elif not NoneList.has_key('ForceDensity'):
                return self.ForceDensity*self.Thickness*self.Width/g
            
        elif not NoneList.has_key('Thickness'):
            
            if not NoneList.has_key('AreaDensity'):
                return self.AreaDensity*self.Thickness
            elif not NoneList.has_key('AreaForceDensity'):
                return self.AreaForceDensity*self.Thickness/g
            
        elif not NoneList.has_key('LinearForceDensity'):
            return self.LinearForceDensity/g

#===============================================================================
    def _CalcForceDensity(self):
        """
        Calculates the force-density if a density was specified
        """
        NoneList = self.NoneList
        if not NoneList.has_key('AreaDensity') and not NoneList.has_key('Thickness'):
            return self.AreaDensity / self.Thickness * g
        elif not NoneList.has_key('AreaForceDensity') and not NoneList.has_key('Thickness'):
            return self.AreaForceDensity / self.Thickness
        else:
            return self.Density*g

#===============================================================================
    def _CalcAreaForceDensity(self):
        """
        Calculates the area force density
        """
        NoneList = self.NoneList
        if not NoneList.has_key('Thickness') and not NoneList.has_key('ForceDensity'):
            return self.ForceDensity * self.Thickness
        elif not NoneList.has_key('Thickness') and not NoneList.has_key('Density'):
            return self.Density * self.Thickness * g
        else:
            return self.AreaDensity*g

#===============================================================================
    def _CalcLinearForceDensity(self):
        """
        Calculates the linear force density
        """
        NoneList = self.NoneList
        if not NoneList.has_key('Thickness') and not NoneList.has_key('Width'):
        
            if not NoneList.has_key('ForceDensity'):
                return self.ForceDensity * self.Thickness * self.Width
            elif not NoneList.has_key('Density'):
                return self.Density * self.Thickness * self.Width * g
        
        elif not NoneList.has_key('Thickness'):
            
            if not NoneList.has_key('AreaForceDensity'):
                return self.AreaForceDensity * self.Thickness
            elif not NoneList.has_key('AreaDensity'):
                return self.AreaDensity * self.Thickness * g
            
        elif  not NoneList.has_key('LinearDensity'):
            return self.LinearDensity*g

#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACMaterial, self).__getattr__(key)
        if ans is not None:
            return ans

        #
        # It is unusual to put Referesh here, but there is no other place to put it
        # for this class
        #
        if self.dirty:
            self.Refresh()
        
        if not self.param.ValidEqVar[key].Valid:
            message = '\nIn order to calculate ' + key
            raise UnderSpecifiedError(self.name,self.param.ValidEqVar[key].Need,self.param.ValidEqVar[key].EqVar, message = message)

        # Try to calculate the value instead
        if key == 'Density':
            return self._CalcDensity()
        elif key == 'AreaDensity':
            return self._CalcAreaDensity()
        elif key == 'ForceDensity':
            return self._CalcForceDensity()
        elif key == 'AreaForceDensity':
            return self._CalcAreaForceDensity()
        elif key == 'LinearDensity':
            return self._CalcLinearDensity()
        elif key == 'LinearForceDensity':
            return self._CalcLinearForceDensity()


#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACMaterial, self).Refresh()
        self.param.refreshing = True
        #
        # Perform all calculations
        #
        self._CheckConsistent()

        self.param.refreshing = False

#===============================================================================
    def _WeakCheckEquation(self, EqVar, Need = 1):
        """
        Does a weak equation check. Only checks the equation if at least nMin variables
        have been specified for the equation 
        
        Inputs:
            EqVar - Variables of the equation to check
            Need  - Minimum number of variables required to check the equation
        """
        #
        # Check the equation only if at least one variable has been specified
        #
        if sum([0 if self.NoneList.has_key(var) else 1 for var in EqVar]) >= Need:
            self._CheckEquation(EqVar, Need = Need)
            Valid = True
        else:
            Valid = False
            
        class Eq:
            def __init__(self):
                self.Valid = False
            pass
        
        ValidEqVar = self.param.ValidEqVar
        
        #
        # Save of if the variable can be calculated
        # If it cannot be calculated, save of the variables and number
        # of those variables need to calculate the variable
        #
        for var in EqVar:
            if not ValidEqVar.has_key(var): ValidEqVar[var] = Eq()
            ValidEqVar[var].Valid = ValidEqVar[var].Valid or Valid   
            if not Valid:
                ValidEqVar[var].EqVar = EqVar
                ValidEqVar[var].Need = Need

#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.

        For example, only Density or ForceDensity should be specified, not both.
        """

        self._WeakCheckEquation(['Density','ForceDensity'])
        self._WeakCheckEquation(['AreaDensity','AreaForceDensity'])       
        self._WeakCheckEquation(['Density','ForceDensity','Thickness','AreaDensity','AreaForceDensity'], 2)
        self._WeakCheckEquation(['LinearDensity','LinearForceDensity','Density','ForceDensity','Thickness','Width'], 3)
        self._WeakCheckEquation(['LinearDensity','LinearForceDensity','AreaDensity','AreaForceDensity','Thickness'], 2)


################################################################################
if __name__ == '__main__':

    from DefaultMaterialsLibrary import PinkFoam

    Heavy = ACMaterial()
#    Heavy.Density = 1000*LBM/IN**3
    Heavy.Thickness = 0.1*IN
#    Heavy.ForceDensity = 0.1 * LBF/IN**3
    Heavy.AreaForceDensity = 0.1 * LBF/IN**2

    Foam = PinkFoam.copy()
    Foam.ForceDensity = 10 * LBF/FT**3
    
#    print Heavy.AreaDensity
    print Heavy.ForceDensity
    print PinkFoam.ForceDensity
    print Foam.ForceDensity
#    print Heavy.Density

