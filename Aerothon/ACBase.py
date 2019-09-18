"""
 Aircraft base classes.

 This is a collection of classes used to create parts of an aircraft.
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
#
# Try to import psycho to improve performance
#
try:
    from psyco.classes import __metaclass__
    import psyco
    psyco.full()
except ImportError:
#    print "You do not have Psycho installed. Psycho can drastically improve performance"
    pass

import numpy as npy
import pylab as pyl
import inspect
from math import sqrt
from scalar.scalar import _scalar, _is_zero
from scalar.units import M, KG, LBF, IN, SEC, ARCDEG, N, HZ, LBM, V, A, OHM, RPM, Ah

################################################################################

# Useful Constants
g = 9.81 * M / SEC**2

################################################################################
class ACError(Exception):
    """ Error handler to setting attributes """
    def __init__(self):

        self.offender = ''
#        stack = inspect.stack()
#        for i in range(len(stack)):
#
#            if stack[i][4] is None:
#                continue

#            caller = stack[i][4][0]
#            if caller.find('raise') > -1:
#                break

#        i=i+1
#        if stack[i][4][0].find('__getattr__'): i=i+1

#        filename = stack[i][1]
#        line     = stack[i][2]
#        offender = stack[i][4][0].strip('\n')

#        self.offender = "\n\nError occurred on line " + str(line) + \
#                        " of\n" + filename + "\n\nErroneous statement: " + offender

#===============================================================================
class ParamKeyError(ACError):
    """ Error handler to setting attributes """
    def __init__(self, key, name):
        ACError.__init__(self)
        self.key = key
        self.name = name

    def __str__(self):
        return "Parameter '" + self.key + \
               "' does not exist in '" + self.name + "'" + self.offender

#===============================================================================
class ListLenError(Exception):
    """ Error handler for list not having the same length """
    def __init__(self, message, name):
        self.msg = message
        self.name = name

    def __str__(self):
        return repr(self.msg + "'" + self.name + "'")

#===============================================================================
class MessageError(ACError):
    """ Error handler for wrong units """
    def __init__(self, message):
        ACError.__init__(self)
        self.msg = message

    def __str__(self):
        return self.msg + self.offender

################################################################################
class OverSpecifiedError(Exception):
    """
    Manages errors for when equations are over specified
    """

#===============================================================================
    def __init__(self,type, count, params, message = ''):

        self.msg = type + " properties over specified.\n"
        self.msg += "Must only specify " + str(count) + " of "
        for i in range(len(params)-1):
            self.msg += "\'" + params[i] + "\', "

        self.msg += "and \'" + params[-1] + "\'"
        self.msg += message

#===============================================================================
    def __str__(self):
        return self.msg

################################################################################
class UnderSpecifiedError(Exception):
    """
    Manages errors for when equations are under specified
    """

#===============================================================================
    def __init__(self,type, count, params, message = ''):

        self.msg = type + " properties are under specified.\n"
        self.msg += "Must specify " + str(count) + " of "
        for i in range(len(params)-1):
            self.msg += "\'" + params[i] + "\', "

        self.msg += "and \'" + params[-1] + "\'"
        self.msg += message

#===============================================================================
    def __str__(self):
        return self.msg

################################################################################
class ACUnitCheck(object):
    """
    A class used to make sure that appropriate units are assigned to variables
    """
#===============================================================================
    def __init__(self, unit, name):
        """
        Inputs:
            unit - A _scalar unit (ex: M**2)
            name - The name associated with the unit (ex: 'Area')
        """
        self.unit = unit
        self.name = name

#===============================================================================
    def CheckUnit(self, key, value, owner):
        """
        Checks that the unit is consistent with self.unit

        Inputs
            key   - The name of the variable checked
            value - The variable assigned to key
            owner - The name of the object who own the variable key
        """
        #
        # A value of None will never have a unit and is a valid asignment
        #
        if value is None:
            return
        #
        # This is the case of a list of scalars
        #
        if hasattr(value, '__len__') and not isinstance(value, _scalar):
            for val in value:
                self._TryUnit(key, val, owner)
            return

        #
        # unit is a single scalar
        #
        self._TryUnit(key, value, owner)


#===============================================================================
    def _TryUnit(self, key, value, owner):
        """
        Test that the unit of unit is consistent with self.unit

        Inputs
            key   - The name of the variable checked
            value - The variable assigned to key
            owner - The name of the object who own the variable key
        """
        message = "'" + key + "' in '" + owner + "' must be in units of " + self.name

        #
        # The units do not matter if it is zero
        #
        if _is_zero(value):
            return
        #
        # Check that they both are or are not _scalars's
        #
        if isinstance(value, _scalar) != isinstance(self.unit, _scalar):
            raise MessageError(message)

        #
        # Try to do a simple operation. If this fails the units are inconsistent.
        #
        try:
            temp = self.unit + value
            return
        except _scalar.InconsistentUnits:
            raise MessageError(message)

#===============================================================================
#
# A collection of unit types
#
#===============================================================================
Unitless           = ACUnitCheck(1             , 'Unitless')
Time               = ACUnitCheck(SEC           , 'Time')
Length             = ACUnitCheck(M             , 'Length')
Area               = ACUnitCheck(M**2          , 'Area')
Volume             = ACUnitCheck(M**3          , 'Volume')
Mass               = ACUnitCheck(KG            , 'Mass')
Force              = ACUnitCheck(KG*M/SEC**2   , 'Force')
Power              = ACUnitCheck(KG*M**2/SEC**3, 'Power')
Torque             = ACUnitCheck(N*M           , 'Torque')
MassTorque         = ACUnitCheck(KG*M          , 'MassTorque')
Velocity           = ACUnitCheck(M/SEC         , 'Velocity')
Angle              = ACUnitCheck(ARCDEG        , 'Angle')
Density            = ACUnitCheck(KG/M**3       , 'Density')
ForceDensity       = ACUnitCheck(N/M**3        , 'ForceDensity')
AreaDensity        = ACUnitCheck(KG/M**2       , 'Density/Area')
AreaForceDensity   = ACUnitCheck(N/M**2        , 'ForceDensity/Area')
LinearDensity      = ACUnitCheck(KG/M          , 'Density/Length')
LinearForceDensity = ACUnitCheck(N/M           , 'ForceDensity/Length')
Pressure           = ACUnitCheck(N/M**2        , 'Pressure')
Frequency          = ACUnitCheck(HZ            , 'Frequency')
PerAngle           = ACUnitCheck(1/ARCDEG      , '1/Angle')
Thrust_SFC         = ACUnitCheck(SEC/M         , 'Thrust Specific Fuel Consumption (PSFC)')
Power_SFC          = ACUnitCheck(SEC**2/M**2   , 'Power Specific Fuel Consumption (PSFC)')
MomentOfInertia    = ACUnitCheck(KG*M**2       , 'Moment of Inertia (Mass*Length**2)')
AreaMomentOfInertia= ACUnitCheck(M**4          , 'Area Moment of Inertia (Length**4)')
Current            = ACUnitCheck(A             , 'Electrical Current')
Resistance         = ACUnitCheck(OHM           , 'Electrical Resistance')
RPMPerVolt         = ACUnitCheck(RPM/V         , 'RPM/Volt')
ElectricPotential  = ACUnitCheck(V             , 'Electric Potential (Voltage)')
ElectricCapacity   = ACUnitCheck(Ah            , 'Electric Capacity (Amp-Hour)')

################################################################################
class ACBase(object):
    """
    A base class for all classes with parameters

    Attributes:
        dirty - Determine if any attributes have been modified and
                a refresh is requires
        NoneList - A dictionary of attributes that have been set to None
        UnitList - A dictionary of units associated with attributes
        param    - A class for storing inAsUnition internal to the class
    """

#===============================================================================
    def __init__(self):

        self.__dict__['name'] = self._getname()
        self.__dict__['dirty'] = True

        #
        # A dictionary of values that have been set to None
        # This is used to implement a calculation instead
        # of a simple get for a variable
        #
        self.__dict__['NoneList'] = {}

        #
        # A dictionary of all the unit associated with a variable
        #
        self.__dict__['UnitList'] = {}

        #
        # A parent classes that will be 'dirtied' if self becomes dirty
        #
        self.__dict__['Parents'] = []
        
        #
        # A class to store internal inAsUnition
        # this will avoid making the class 'dirty'
        # when performing internal calculations
        #
        class ParamClass:
            pass
        self.__dict__['param'] = ParamClass()
        self.param.refreshing = False

#===============================================================================
    def SetDirty(self):
        """
        Used to set the present class and the parent class to be dirty
        """
        if self.param.refreshing:
            return
        
        self.__dict__['dirty'] = True
        for Parent in self.__dict__['Parents']:
            Parent.SetDirty()
            
#===============================================================================
    def __setattr__(self, key, value):
        """
        Overload the '.' on all objects such that new
        attributes cannot accidentally be created. Thus only
        attributes created in __init__ will ever be available
        in the class.
        """
        #
        # Make sure no new variables are created besides those set in __init__
        #
        if not self.__dict__.has_key(key) and not self.NoneList.has_key(key):
            raise ParamKeyError(key, self.__dict__['name'])

        #
        # Make sure that the proper units have been assigned
        #
        if self.UnitList.has_key(key):
            self.UnitList[key].CheckUnit(key, value, self.__class__.__name__)

        #
        # Assign the value if it is different from what is already there
        #
        self.__dict__[key] = value

        #
        # If the value is set to None, pull it out of __dict__ and put it in NoneList
        # so that __getattr__ gets called when the variable is asked for
        #
        if value is None:
            self.NoneList[key] = None
            del self.__dict__[key]
        elif self.NoneList.has_key(key):
            del self.NoneList[key]
            
        if key != 'dirty':
            self.SetDirty()

#===============================================================================
    def __getattr__(self, key):
        """
        Limit the access to only what exists in the __dict__ dictionary
        """

        #
        # If the value was set to None, it got taken out of the __dict__ and
        # put in NoneList to force the call to __getattr__
        #
        if self.NoneList.has_key(key):
            return None

        #
        # The key is not in the None list or __dict__
        #
        raise ParamKeyError(key, self.__dict__['name'])

#===============================================================================
    def __delattr__(self, key):
        #
        # Cannot delete anything from the dictionary
        #
        message = "Cannot del '" + key + "' from " + self.__class__.__name__
        raise MessageError(message)

#===============================================================================
    def _CheckSameLen(self, L1, name1, L2,  name2):
        "Determines if two lists are the same length."

        message = name1 + " and " + name2 + " must be of same length in "
        try:
            if len(L1) != len(L2):
                raise ListLenError(message, self.__class__.__name__)
        except TypeError:
            try:
                len(L1)
            except TypeError:
                L1 = [L1]

            try:
                len(L2)
            except TypeError:
                L2 = [L2]

            if len(L1) != len(L2):
                raise ListLenError(message, self.__class__.__name__)

        return L1, L2
    
#===============================================================================
    def _CheckEquation(self, EqVar, message = '', Need = None):
        """
        Verifies that an equation is not over or under specified
        
        Inputs:
            EqVar - A list of variables in the equation
            message - An optional additional error message
        """
        NoneList = self.NoneList
        if Need is None: Need = len(EqVar) - 1 #Assume equation needs one variable unspecified
        
        #
        # Check that the equation is not over sepcified. If non of the variables are
        # in the NoneList, they have have all been specified
        #
#        if all([not NoneList.has_key(var) for var in EqVar]):
#            raise OverSpecifiedError(self.name,Need,EqVar,message)

        #
        # Check that the equation is not under specified. Count the number of 
        # variables that are not in the NoneList. This should add up to the 
        # total number of variables needed
        #
        nVar = sum([0 if NoneList.has_key(var) else 1 for var in EqVar])
        if nVar < Need:
            raise UnderSpecifiedError(self.__class__.__name__, Need, EqVar, message)
        elif nVar > Need:
            raise OverSpecifiedError(self.__class__.__name__, Need, EqVar, message)
            

#===============================================================================
    def _getname(self):
        """
        Returns the name of the instance of the class when it was created

        WARNING: This should only be called from __init__ in ACBase.
                 Calling it any other time will not return a name
        """
        return self.__class__.__name__

        stack = inspect.stack()
        for i in range(len(stack)):

            if stack[i][4] is None:
                continue

            caller = stack[i][4][0]
            caller = caller.replace(' ','')
            caller = caller.split('=')
            name = None
            LR = []
            for j in range(len(caller)):
                LR.append(caller[j].split(','))

            if len(LR) == 1:
                continue

            for j in range(len(LR[0])):
                if LR[1][j].find(self.__class__.__name__) > -1:
                   name = LR[0][j]

            if name is not None:
               break

        if name is None:
            return self.__class__.__name__

        name = name.strip('self.')
        if name.find("'") != -1:
           name = name[name.find("'")+1:name.rfind("'")]

        return name

#===============================================================================
    def MakeList(self, Value):
        """
        Converts value into a list if it is not already
        
        Inputs:
            Value - The value to convert to a list
        """
        try:
            a = len(Value)
            return Value
        except:
            return [Value]

#===============================================================================
    def ToUnumList(self,List,Unit):
        """
        Converts a list of unums to a unum list
        
        Inputs:
            List - The list to convert
            Unit - The unit of the unums in the list
        """
        return npy.array([l / Unit for l in List])*Unit

#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
#        print "Refreshing ", self.name
        self.__dict__['dirty'] = False
        self.param.refreshing = True
        #
        # Perform all calculations
        #
        self.param.refreshing = False


################################################################################
class ACComponent(ACBase):
    """
    A base class for aircraft components.

    Components are sections of an aircraft part.

    For example, of a wing is an aircraft part, then an
    aileron will be a component of the wing.

    Note that everything should always be drawn in inches

    Attributes
       X        - Origin of the object
       Drawings - A dictionary of ComponentDrawing structures to be drawn
       ShowTop
       ShowSide
       ShowFront - Determines if a side is visible
    """

    class ComponentDrawing:
        """
        A class for drawing a part of a component.
        """
        def __init__(self):
            
            # Color of the component
            self.color = 'b'

            # Plot style
            self.style = '-'

            # Arrays for drawing the different views of the component
            # Front view
            self.yf = npy.array([])
            self.zf = npy.array([])

            # Side view
            self.xs = npy.array([])
            self.zs = npy.array([])

            # Top view
            self.xt = npy.array([])
            self.yt = npy.array([])

#===============================================================================
    def __init__(self):
        ACBase.__init__(self)

        # Initialize the global location array
        self.__dict__['X']=[0 * IN, 0 * IN, 0 * IN]

        class Drawings:
            pass

        # A dictionary of component drawings
        self.__dict__['Drawings'] = Drawings()

        # Determins if the component is symmetric about the y = 0 plane (which is an x-z plane)
        self.__dict__['Symmetric'] = False

        # Determinse if  side is shown
        self.__dict__['ShowTop'] = True
        self.__dict__['ShowSide'] = True
        self.__dict__['ShowFront'] = True

#===============================================================================
    def CopyVisibility(self, other):
        """
        Copy's the visibility of other to self
        """
        self.ShowFront = other.ShowFront
        self.ShowSide  = other.ShowSide
        self.ShowTop   = other.ShowTop

#===============================================================================
    def GetX(self):
        """
        Returns X as scalar with an array i.e. (x1, x2, x3) [in]
        """
        return npy.array([self.X[i] / IN for i in range(3)])*IN

#===============================================================================
    def _DrawFront(self, Drawing):
        """
        Draws a front view of the component.

        It is assumed the the desired plotting location has been
        chosen before the call
        """

        x = Drawing.yf
        y = Drawing.zf
        color = Drawing.color
        style = Drawing.style

        pyl.plot(x, y, color = color,  linestyle = style)
        if self.Symmetric:
            pyl.plot(-x, y, color = color,  linestyle = style)

#===============================================================================
    def _DrawSide(self, Drawing):
        """
        Draws a side view of the component.

        It is assumed the the desired plotting location has been
        chosen before the call
        """

        x = Drawing.xs
        y = Drawing.zs
        color = Drawing.color
        style = Drawing.style

        pyl.plot(x, y, color = color,  linestyle = style)

#===============================================================================
    def _DrawTop(self, Drawing):
        """
        Draws a top view of the component.

        It is assumed the the desired plotting location has been
        chosen before the call
        """

        y = Drawing.yt
        x = Drawing.xt
        color = Drawing.color
        style = Drawing.style

        pyl.plot(x, y, color = color,  linestyle = style)
        if self.Symmetric:
            pyl.plot(x, -y, color = color,  linestyle = style)

#===============================================================================
    def SWet(self):
        """
        Default value for wetted area of structural components.

        Returns 0. if no wetted area is defined.
        """
        # This function returns zero in case the structural component
        # does not have an effective wetted area. For example if it is
        # a fuel tank inside the aircraft, it does not contribute to the wetted area
        return 0.

#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """
        Draws this component

        Inputs:
            fig - Integer number of figure to draw into
            top - Subplot for top view
            side - Subplot for side view
            front - Subplot for front view
        """
        # Update the drawing arrays if needed
        if self.dirty:
            self.Refresh()
        
        pyl.figure(fig, figsize=(10,10))
        pyl.subplots_adjust(left=0.1)
            
        if self.ShowTop:
            pyl.subplot(top) #Top view
            pyl.axis('equal')
            pyl.title('Top View')
            pyl.xlabel('x (in)')
            pyl.ylabel('y (in)')

            for v in self.Drawings.__dict__.itervalues():
                self._DrawTop(v)

        if self.ShowSide:
            pyl.subplot(side) #Side view
            pyl.axis('equal')
            pyl.title('Side View')
            pyl.axhline(color = 'k')
            pyl.xlabel('x (in)')
            pyl.ylabel('z (in)')

            for v in self.Drawings.__dict__.itervalues():
                self._DrawSide(v)

        if self.ShowFront:
            pyl.subplot(front) #Front view
            pyl.axis('equal')
            pyl.title('Front View')
            pyl.axhline(color = 'k')
            pyl.xlabel('y (in)')
            pyl.ylabel('z (in)')
            #
            # The horizontal axis must be flipped here to such that the view looks
            # from the nose to the tail of the body
            #
            if not pyl.gca().xaxis_inverted():
                pyl.gca().invert_xaxis()

            for v in self.Drawings.__dict__.itervalues():
                self._DrawFront(v)

#===============================================================================
    def UpdateDrawing(self):
        """
        This method should be used to update all drawings of the object
        it defaults to calling Refresh if it is not implemented.
        """
        self.Refresh()

################################################################################
class ACStructural(ACComponent):
    """
    A class for aircraft structures.

    This includes bulkheads, servos, and other parts of the aircraft.
    """

#===============================================================================
    def __init__(self):
        super(ACStructural, self).__init__()

        UnitList = self.UnitList
        self.NoneList['Weight'] = None  ;   UnitList['Weight'] = Force
        self.__dict__['WeightGroup'] = 'None'

#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACStructural, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'Weight':
            return self._CalcWeight()
        #
        # If no calculation exist return None
        #
        return None
    
#===============================================================================
    def _CalcWeight(self):
        """
        A fall back CalcWeight that will just return zero
        """
        return 0*LBF

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the current part to the Weight Table
        
        Input:
            PartName    - The name give to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WeightTable[PartName] = self

#===============================================================================
    def _XYZOrient(self):
        """
        Calculates orientation of a vector and the X, Y and Z axes.
        0 = 0 deg, 1 = 90 deg
        """

        Axial = npy.array(self.Axis)
        #
        # Normalize the axis vector
        #
        Axial = Axial/sqrt(sum(Axial**2))

        #
        # Create some unit vectors to get the radial vector
        #
        x = npy.array([1,0,0])
        y = npy.array([0,1,0])
        z = npy.array([0,0,1])

        #
        # Cross the vector with an appropriate vector that is not parallel
        # to the axis
        #
        if abs(pyl.dot(Axial,z)) != 1:
            Radial1 = pyl.cross(z,Axial)
        elif abs(pyl.dot(Axial,x)) != 1:
            Radial1 = pyl.cross(x,Axial)
        elif abs(pyl.dot(Axial,y)) != 1:
            Radial1 = pyl.cross(y,Axial)
        

        Radial1 = Radial1/sqrt(sum(Radial1**2))

        #
        # Get a second radial orientation to draw in
        #
        Radial2 = pyl.cross(Axial,Radial1)
        Radial2 = Radial2/sqrt(sum(Radial2**2))

        #
        # Return the vectors
        #
        return Axial, Radial1, Radial2

#===============================================================================
    def _XFormMOI(self,Ix,Iy,Iz):
        """
        Transfer MOI into Global Coordinates
        """
        
        # Correct Component MOI for orientation...
        # LI, LJ, LK are the I,J,K vectors in the Components' local coordinates
        # Taken from:  http://www.kwon3d.com/theory/transform/transform.html
        LI,LJ,LK = self._XYZOrient()
        GI = npy.array([1,0,0])
        GJ = npy.array([0,1,0])
        GK = npy.array([0,0,1])
        
        dot = npy.dot
        
        XFormMatrix = npy.array([[dot(GI,LI),dot(GJ,LI),dot(GK,LI)],
                                 [dot(GI,LJ),dot(GJ,LJ),dot(GK,LJ)],
                                 [dot(GI,LK),dot(GJ,LK),dot(GK,LK)]])
        
        lambdaX = dot(XFormMatrix, GI)
        lambdaY = dot(XFormMatrix, GJ)
        lambdaZ = dot(XFormMatrix, GK)
        
        # Rotated MOI
        GIx = Ix*lambdaX[0]**2 + Iy*lambdaX[1]**2 + Iz*lambdaX[2]**2
        GIy = Ix*lambdaY[0]**2 + Iy*lambdaY[1]**2 + Iz*lambdaY[2]**2
        GIz = Ix*lambdaZ[0]**2 + Iy*lambdaZ[1]**2 + Iz*lambdaZ[2]**2
        
        return self.ToUnumList((GIx, GIy, GIz), LBM*IN**2)
       
################################################################################
if __name__ == '__main__':

    #Construct a part with a box and draw it

    box1 = ACComponent()
    box1.Drawings.box1 = ACComponent.ComponentDrawing()
    box1.Drawings.box1.yf = npy.array([0, 1, 1, 0, 0])
    box1.Drawings.box1.zf = npy.array([0, 0, 1, 1, 0])
    box1.Drawings.box1.xs = npy.array([0, 0.5, 0.5, 0, 0])
    box1.Drawings.box1.zs = npy.array([0, 0, 0.5, 0.5, 0])
    box1.Drawings.box1.xt = npy.array([0, 0.75, 0.75, 0, 0])
    box1.Drawings.box1.yt = npy.array([0, 0, 0.75, 0.75, 0])
    box1.Drawings.box1.color = 'r'
    box1.Symmetric = True

    box1.ShowFront = False
    
    box1.X = [IN, IN, IN]

    base = ACBase()
    
    box1.Draw(1)
    pyl.show()


