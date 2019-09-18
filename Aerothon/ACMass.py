"""
A collection of classes for fuselage parts
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACComponent, ACStructural
import numpy as npy
import pylab as pyl
from operator import mul
from matplotlib.patches import Ellipse
from scalar.units import IN, LBM
from math import sqrt, sin, cos, tan
from ACBase import Length, Mass, Unitless, g

################################################################################
class ACMassBox(ACStructural):
    """
    A class for fuselage parts that are box shaped

    Components locations are defined at the center of the front faces
    bottom edge

    i.e. servos, fuel tank, etc.

        LWH     - [Length, Width, Height] dimensions of component
        X       - [x,y,z] coordinates of component
        Weight  - Weight of component
        Axis    - Orientation vector [x,y,z] normalized to 1
        Xat     - (xi, eta, zeta) Fractional position of X
        Color   - The color
    """
#===============================================================================
    def __init__(self):
        ACStructural.__init__(self)
        UnitList = self.UnitList

        self.__dict__['LWH']  = [1 * IN, 1 * IN, 1 * IN]; UnitList['Dim']  = Length
        self.__dict__['Axis'] = [1, 0, 0]               ; UnitList['Axis'] = Unitless
        self.__dict__['Xat']  = [0, 0.5, 0]             ; UnitList['Xat'] = Unitless
        self.__dict__['Color'] = 'b'

        #
        # Setup the drawings for each face of the mass box
        #
        self.Drawings.Front  = ACComponent.ComponentDrawing()
        self.Drawings.Back   = ACComponent.ComponentDrawing()
        self.Drawings.Top    = ACComponent.ComponentDrawing()
        self.Drawings.Bottom = ACComponent.ComponentDrawing()
        self.Drawings.Left   = ACComponent.ComponentDrawing()
        self.Drawings.Right  = ACComponent.ComponentDrawing()

#===============================================================================
    def _GetOrient(self):
        """
        Creates a MassBox that is aligned with the 'Axis' vector
        """

        length, width, height = self._XYZOrient()

        length = self.LWH[0]*length
        width  = self.LWH[1]*width
        height = self.LWH[2]*height

        return length, width, height

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
        #
        # Perform all calculations
        #
        Drawings = self.Drawings.__dict__

        length, width, height = self._XYZOrient()
        length *= self.LWH[0] / IN
        width  *= self.LWH[1] / IN
        height *= self.LWH[2] / IN
        X = self.GetX() / IN
        xi, eta, zeta = self.Xat[0], self.Xat[1], self.Xat[2] 

        class Rectangle:
            def __init__(self, length, height, origin):
                self.NodesX = npy.array([0, length[0], length[0] + height[0], height[0], 0]) + origin[0]
                self.NodesY = npy.array([0, length[1], length[1] + height[1], height[1], 0]) + origin[1]
                self.NodesZ = npy.array([0, length[2], length[2] + height[2], height[2], 0]) + origin[2]
        
        origin = X - length*xi - width*eta - height*zeta
        Sides = {}
        Sides['Front']  = Rectangle(width, height, origin)
        Sides['Back']   = Rectangle(width, height, origin + length)
        Sides['Top']    = Rectangle(length, width, origin + height)
        Sides['Bottom'] = Rectangle(length, width, origin)
        Sides['Left']   = Rectangle(length, height, origin)
        Sides['Right']  = Rectangle(length, height, origin + width)

        for key, rec in Sides.iteritems():
            Drawings[key].yf = rec.NodesY
            Drawings[key].zf = rec.NodesZ
            Drawings[key].xs = rec.NodesX
            Drawings[key].zs = rec.NodesZ
            Drawings[key].xt = rec.NodesX
            Drawings[key].yt = rec.NodesY
            Drawings[key].color = self.Color

        self.param.refreshing = False

#===============================================================================
    def Center(self):
        """
        Returns the global center location of a MassBox as (x,y,z)
        """
        if self.dirty:
            self.Refresh()

        length, width, height = self._GetOrient()
        
        xi, eta, zeta = self.Xat[0], self.Xat[1], self.Xat[2] 
        X = self.GetX()
        
        return length*(0.5 - xi) + width*(0.5 - eta) + height*(0.5 - zeta) + X
             
#===============================================================================
    def CG(self):
        """
        Returns the global CG location of a MassBox as (x,y,z)
        """
        if self.dirty:
            self.Refresh()

        CG = self.Center()
        
        if self.Symmetric:
            CG[1] = 0*IN

        return CG
    
#===============================================================================
    def MOI(self, Local = False):
        """
        Returns the moments of inertia about its CG as (Ixx, Iyy, Izz)
        """
        if self.dirty:
            self.Refresh()

        CG = self.CG()
        [length, width, height] = self.LWH
        Mass = self.Weight / g
        
        if self.Symmetric:
            Mass *= 2
        
        Ix = 1/12 * Mass * (width**2 + height**2) 
        Iy = 1/12 * Mass * (length**2 + height**2) 
        Iz = 1/12 * Mass * (width**2 + length**2)
        
        MOI = self._XFormMOI(Ix, Iy, Iz)
        
        if self.Symmetric:
            dCG = self.CG() - self.Center()
            MOI[0] += Mass*(dCG[1]**2 + dCG[2]**2)#Ix
            MOI[1] += Mass*(dCG[0]**2 + dCG[2]**2)#Iy
            MOI[2] += Mass*(dCG[1]**2 + dCG[0]**2)#Iz

        if Local is False:
            return MOI
        else:
            return self.ToUnumList((Ix, Iy, Iz), LBM*IN**2)

#===============================================================================
    def Volume(self):
        """
        Computes the volume of the mass box
        """
        if self.dirty:
            self.Refresh()
            
        return self.LWH[0]*self.LWH[1]*self.LWH[2]
        
################################################################################
class ACMassCyl(ACStructural):
    """
    A class for fuselage parts that are cylindrically shaped

    i.e. wheels, tail boom, etc.

        LenDi   - [length, diameter]
        Xat     - (xi, eta, zeta) Fractional position of X
        Color   - The color
    """
#===============================================================================
    def __init__(self):
        ACStructural.__init__(self)
        UnitList = self.UnitList

        self.__dict__['LenDi']  = [1 * IN, 1 * IN]  ;   UnitList['LenDi']  = Length
        self.__dict__['Axis']   = [1, 0, 0]         ;   UnitList['Axis']   = Unitless
        self.__dict__['Xat']    = [0, 0.5, 0]             ; UnitList['Xat'] = Unitless
        self.__dict__['Color'] = 'b'

        #
        # This is new, it is a dictionary of drawings
        # You create what should be drawn and fill it's arrays instead
        #
        self.Drawings.Side1 = ACComponent.ComponentDrawing()
        self.Drawings.Side2 = ACComponent.ComponentDrawing()

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
        super(ACMassCyl, self).Refresh()
        self.param.refreshing = True
        #
        # Perform all calculations
        #
        Side1 = self.Drawings.Side1
        Side2 = self.Drawings.Side2
        
        Side1.color = self.Color
        Side2.color = self.Color
        
        xi, eta, zeta = self.Xat[0], self.Xat[1], self.Xat[2] 
        L = self.LenDi[0]
        R = self.LenDi[1] / IN / 2

        Axial, Radial1, Radial2 = self._XYZOrient()
                
        X = (self.GetX() - Axial*L*xi) / IN

        Axial *= L / IN
        Radial1 *= R
        Radial2 *= R

        #
        # Top
        #
        Side1.xt = npy.array([-Radial1[0], Axial[0] - Radial1[0], Axial[0] + Radial1[0], Radial1[0], -Radial1[0]]) + X[0]
        Side1.yt = npy.array([-Radial1[1], Axial[1] - Radial1[1], Axial[1] + Radial1[1], Radial1[1], -Radial1[1]]) + X[1]

        Side2.xt = npy.array([-Radial2[0], Axial[0] - Radial2[0], Axial[0] + Radial2[0], Radial2[0], -Radial2[0]]) + X[0]
        Side2.yt = npy.array([-Radial2[1], Axial[1] - Radial2[1], Axial[1] + Radial2[1], Radial2[1], -Radial2[1]]) + X[1]

        #
        # Front
        #
        Side1.yf = npy.array([-Radial1[1], Axial[1] - Radial1[1], Axial[1] + Radial1[1], Radial1[1], -Radial1[1]]) + X[1]
        Side1.zf = npy.array([-Radial1[2], Axial[2] - Radial1[2], Axial[2] + Radial1[2], Radial1[2], -Radial1[2]]) + X[2]

        Side2.yf = npy.array([-Radial2[1], Axial[1] - Radial2[1], Axial[1] + Radial2[1], Radial2[1], -Radial2[1]]) + X[1]
        Side2.zf = npy.array([-Radial2[2], Axial[2] - Radial2[2], Axial[2] + Radial2[2], Radial2[2], -Radial2[2]]) + X[2]

        #
        # Side
        #
        Side1.xs = npy.array([-Radial1[0], Axial[0] - Radial1[0], Axial[0] + Radial1[0], Radial1[0], -Radial1[0]]) + X[0]
        Side1.zs = npy.array([-Radial1[2], Axial[2] - Radial1[2], Axial[2] + Radial1[2], Radial1[2], -Radial1[2]]) + X[2]

        Side2.xs = npy.array([-Radial2[0], Axial[0] - Radial2[0], Axial[0] + Radial2[0], Radial2[0], -Radial2[0]]) + X[0]
        Side2.zs = npy.array([-Radial2[2], Axial[2] - Radial2[2], Axial[2] + Radial2[2], Radial2[2], -Radial2[2]]) + X[2]

        self.param.refreshing = False

#===============================================================================
    def _GetOrient(self):
        """
        Creates a MassBox that is aligned with the 'Axis' vector
        """

        L = self.LenDi[0]
        R = self.LenDi[1] / 2

        Axial, Radial1, Radial2 = self._XYZOrient()

        Axial   = L*Axial
        Radial1 = R*Radial1
        Radial2 = R*Radial2
        
        return Axial, Radial1, Radial2

#===============================================================================
    def Center(self):
        """
        Returns the global center location of a MassCyl as (x,y,z)
        """

        Axial, Radial1, Radial2 = self._GetOrient()

        X = self.GetX()
        xi = self.Xat[0]
                
        return Axial *(0.5 - xi) + X

#===============================================================================
    def CG(self):
        """
        Returns the global CG location of a MassCyl as (x,y,z)
        """

        Axial, Radial1, Radial2 = self._GetOrient()

        CG = self.Center()
        
        if self.Symmetric:
            CG[1] = 0*IN
        
        return CG

#===============================================================================
    def MOI(self, Local = False):
        """
        Returns the moments of inertia about its CG as (Ixx, Iyy, Izz)
        """
        if self.dirty:
            self.Refresh()

        CG = self.CG()
        length = self.LenDi[0]
        radius = self.LenDi[1] / 2
        Mass = self.Weight / g
        
        if self.Symmetric:
            Mass *= 2
        
        Ix = 1/2 * Mass * (radius**2) 
        Iy = Iz = 1/12 * Mass * (3 * radius**2 + length**2) 

        MOI = self._XFormMOI(Ix, Iy, Iz)
        
        if self.Symmetric:
            dCG = self.CG() - self.Center()
            MOI[0] += Mass*(dCG[1]**2 + dCG[2]**2)#Ix
            MOI[1] += Mass*(dCG[0]**2 + dCG[2]**2)#Iy
            MOI[2] += Mass*(dCG[1]**2 + dCG[0]**2)#Iz
            
        if Local is False:
            return MOI
        else:
            return self.ToUnumList((Ix, Iy, Iz), LBM*IN**2)

#===============================================================================
    def _DrawFront(self, Drawing):
        """
        Draws a front view of the component.
        """
        super(ACMassCyl, self)._DrawFront(Drawing)

        #
        # Do any additional drawing for a front view here
        # This is called from ACComponent
        #

        Axial, Radial1, Radial2 = self._XYZOrient()

        if abs(Axial[0]) > 0.:
            xi, eta, zeta = self.Xat[0], self.Xat[1], self.Xat[2] 
            #
            # Remember X is [Unum, Unum, Unum]
            #
            L = self.LenDi[0]
            D = self.LenDi[1] / IN
            R = D / 2
            X = (self.GetX() - Axial*L*xi) / IN

            Axial *= L / IN
            Radial1 *= R 
            Radial2 *= R 

            #
            # Get the current axis so the ellipse patch can be added to it
            #
            FrontAxis = pyl.gca()

            #
            # Create the ellipse
            #
            y = X[1]
            z = X[2]
            W = D # Fix later
            H = D # Fix later
            el1 = Ellipse( (y, z), W, H, 0.0, fill=False, edgecolor=Drawing.color)
            el2 = Ellipse( (y + Axial[1], z + Axial[2] ), W, H, 0.0, fill=False, edgecolor=Drawing.color)
            #
            # Add it to the axis
            #
            FrontAxis.add_patch(el1)
            FrontAxis.add_patch(el2)
            if self.Symmetric:
                el1 = Ellipse( (-y, z), W, H, 0.0, fill=False, edgecolor=Drawing.color)
                el2 = Ellipse( (-(y + Axial[1]), z + Axial[2] ), W, H, 0.0, fill=False, edgecolor=Drawing.color)
                FrontAxis.add_patch(el1)
                FrontAxis.add_patch(el2)
                    

#===============================================================================
    def _DrawSide(self, Drawing):
        """
        Draws a side view of the component.
        """
        ACStructural._DrawSide(self, Drawing)

        #
        # Do any additional drawing for a side view here
        # This is called from ACComponent
        #
        Axial, Radial1, Radial2 = self._XYZOrient()

        if abs(Axial[1]) > 0.:
            xi, eta, zeta = self.Xat[0], self.Xat[1], self.Xat[2] 
            #
            # Remember X is [Unum, Unum, Unum]
            #
            L = self.LenDi[0]
            D = self.LenDi[1] / IN
            R = D / 2
            X = (self.GetX() - Axial*L*xi) / IN

            Axial *= L / IN
            Radial1 *= R
            Radial2 *= R 
            
            #
            # Get the current axis so the ellipse patch can be added to it
            #
            FrontAxis = pyl.gca()
            
            #
            # Create the ellipse
            #
            x = X[0]
            z = X[2]
            W = D # Fix later
            H = D # Fix later
            el1 = Ellipse( (x, z), W, H, 0.0, fill=False, edgecolor=Drawing.color)
            el2 = Ellipse( (x + Axial[0], z + Axial[2] ), W, H, 0.0, fill=False, edgecolor=Drawing.color)
            #
            # Add it to the axis
            #
            FrontAxis.add_patch(el1)
            FrontAxis.add_patch(el2)              

#===============================================================================
    def _DrawTop(self, Drawing):
        """
        Draws a top view of the component.
        """
        ACStructural._DrawTop(self, Drawing)

        #
        # Do any additional drawing for a top view here
        # This is called from ACComponent
        #
        Axial, Radial1, Radial2 = self._XYZOrient()

        if abs(Axial[2]) > 0.:
            xi, eta, zeta = self.Xat[0], self.Xat[1], self.Xat[2]
            #
            # Remember X is [Unum, Unum, Unum]
            #
            L = self.LenDi[0]
            D = self.LenDi[1] / IN
            R = D / 2
            X = (self.GetX() - Axial*L*xi) / IN

            Axial *= L / IN
            Radial1 *= R 
            Radial2 *= R

            #
            # Get the current axis so the ellipse patch can be added to it
            #
            FrontAxis = pyl.gca()

            #
            # Create the ellipse
            #
            x = X[0]
            y = X[1]
            W = D # Fix later
            H = D # Fix later
            el1 = Ellipse( (x, y), W, H, 0.0, fill=False, edgecolor=Drawing.color)
            el2 = Ellipse( (x + Axial[0], y + Axial[1] ), W, H, 0.0, fill=False, edgecolor=Drawing.color)
            #
            # Add it to the axis
            #
            FrontAxis.add_patch(el1)
            FrontAxis.add_patch(el2)
            if self.Symmetric:
                el1 = Ellipse( (x, -y), W, H, 0.0, fill=False, edgecolor=Drawing.color)
                el2 = Ellipse( (x + Axial[0], -(y + Axial[1]) ), W, H, 0.0, fill=False, edgecolor=Drawing.color)
                FrontAxis.add_patch(el1)
                FrontAxis.add_patch(el2)
 

################################################################################
if __name__ == '__main__':

##    BHF = ACMassBox()
##    BHF.X = [1 * IN, 0 * IN, 6 * IN]
##    BHF.Dim = [1/8 * IN, 2 * IN, 2 * IN]

    BHM = ACMassBox()
    BHM.X = [1 * IN, 20 * IN, 1 * IN]
    BHM.LWH = [1 * IN, 2 * IN, 3 * IN]
    BHM.Axis = [1, 0, 0]
    BHM.Xat = [0.5, 0.5, 0.5]
    BHM.Weight = 3*LBM*g
    BHM.Symmetric = True
    print BHM.CG()
    print "Note:  These two MOI should match"
    print BHM.MOI()
    print BHM.MOI(Local = True)

##    BHA = ACMassBox()
##    BHA.X = [10 * IN, 0 * IN, 3 * IN]
##    BHA.Dim = [1/8 * IN, 5 * IN, 5 * IN]

    Wheel1 = ACMassCyl()
    Wheel1.X = [0 * IN, 0 * IN, 0 * IN]
    Wheel1.LenDi = [5 * IN, 1 * IN]
    Wheel1.Axis = [0, 1, 0]
    Wheel1.Weight = 1*LBM*g
    print Wheel1.CG()
    print "Note:  These two should be off by 90 degrees (Swap Ixx and Iyy)"
    print Wheel1.MOI()
    print Wheel1.MOI(Local = True)

##    BHF.Draw()
    BHM.Draw()
##    BHA.Draw()
##    Wheel1.Draw()
    pyl.show()














