"""
A collection of classes for fuselages. 

NOTE: The Micro2018 version required trapezoidal bulkheads, so I've patched the code 
to allow this. To switch to a trapezoidal bulkhead,
set 'Width' to -1 and then define 'WidthBot' and 'WidthTop' to create the new parameters.

                                                -Caleb Bisig
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACBase, ACComponent, ACStructural, Length, Unitless, g, Force
from ACMass import ACMassBox
from ACMaterial import ACMaterial
from ACPropulsion import ACPropulsion
import AeroUtil
from scalar.units import IN, LBF, SEC, LBM, SLUG, FT, RAD
from scalar.units import AsUnit
import cmath as math
import pylab as pyl
import numpy as npy
################################################################################
"""
NOTE: I had to edit a number of things to add in the trapezoid, made a list of sections to keep 
track of where I needed to make changes.

(Also for anyone in the future who wants to figure out this code for other shapes
or edit a new custom team version of ACFuselage by copying this code and renaming it:)

        --TABLE OF CONTENTS--
lINE---TOPIC---DESCRIPTION (** = "Change required to allow for Trapezoidal Bulkheads"
XX99        ---ACFuselage ---Creates the ACFuselage class, that's all.
    X111        ---Bulkhead   ---Creates BULKHEAD class 
        X120**      ---__init__   ---Creates list of attributes available for Bulkhead (**added attributes WidthBot and WidthTop)
        X141**      ---Corners    ---Calculate CORNER POINTS from Width & Height (** or WidthBot & WidthTop if Width=-1)
        X162        ---Top          ---Defines the TOP HEIGHT of the given bulkhead
        X170        ---Bottom       ---Defines the BOTTOM HEIGHT of the given bulkhead
        X178**      ---Area         ---Bulkhead Face AREA (**Added Trapezoidal Area Calculation)
        CX193**      ---CG           ---Bulkhead CG (**Added Trapezoidal CG Calculation)
        X208**        ---MOI          ---Bulkhead MOMENT OF INERTIA (**Assumed a formula for thin walled trapezoid, see note in MOI section.)
        X255        ---GetWeight    ---Gets Weight found in _CalcWeight
        X267**        ---_CalcWeight  ---Determine Weight from area*density. (**Added pure Area*density line for trapezoids)
        X282**        ---Refresh    ---Something about refreshing variables? (*In any event, fixed Trapezoid coordinates too.)
    X335        ---Section  ---Creates SECTION class between bulkheads.
        X351        ---Face    ---Creates FACE class to calculate and draw faces.
            X352**      ---__init__   ---Creates a normal vector to the bulkhead face. (Added an ASCII image to help) 
            X395        ---PointLocation  ---Creates interpolating points across top and bottom, with a sweeping vector.
        X431        ---Component  ---Lists attributes of placed objects, like servos.
            X440**      ---__init__   ---For defining section component is found in. 
            X457        ---Refresh  ---For defining section face where component attaches.
        X484**      ---__init__   ---Defines "Front" and "Back" Bulkheads for future section creation
            X529        ---Faces    ---Faces of Components (not sections just yet)
        X544        ---AddComponent    ---Defines position information of the component
        X588        ---PointLocation    ---Defines vectors between FrontBulkhead and BackBulkhead (Added some ASCII art to explain)
        X538        ---SWet
        XXXX        ---Top
        XXXX        ---Bottom
        XXXX**      ---Height
        XXXX**      ---CG
        X670**      ---MOI
        XXXX        ---_CalcComponentWeight
        X700        ---_CalcSectionWeight
        XXXX        ---_CalcWeight
        XXXX        ---_CalcLength
        XXXX        ---__getattr__
        XXXX        ---__getattribute__
        XXXX        ---Draw
        XXXX        ---AddToWeightTable
        XXXX        ---Refresh
    X932        ---ACPayload
        XXXX        ---__init__   ---Creates list of attributes available for 
        XXXX        ---Refresh
    1012        ---__init__   ---Creates list of attributes available for 
    XXXX        ---AddSection
    XXXX        ---AddToWeightTable
    XXXX        ---Swet
    XXXX        ---_CalcWeight
    XXXX        ---CG
    XXXX        ---MOI
    1222        ---AircraftCG
    XXXX        ---Top
    XXXX        ---Bottom
    XXXX        ---Height
    XXXX        ---GetNoseBulk
    XXXX        ---GetNoseSection
    XXXX        ---_FrontalArea
    1323        ---_Length
    XXXX        ---CD
    XXXX        ---Draw
    1415        ---Refresh
    1463       ---FinalOutputs





"""
################################################################################
class ACFuselage(ACStructural):
    """
    A class for describing a fuselage
    
    Attributes:
        Xcg        - The desired CG location of the aircraft
        XcgSection - The section which should contain the CG (first section by default)
        XcgSecFrac - Fraction of the XcgSection where the CG is positioned (0.5 by default)
        Sections - A class of sections
        
    """            
#===============================================================================
    class BulkHead(ACStructural):
        """
        A bulk head class
        
        Attributes:
            Width    - The width of the bulkhead
            Height   - The height of the bulkhead
            Material - The material of the bulkhead
        """
        def __init__(self):
            super(ACFuselage.BulkHead, self).__init__()
            UnitList = self.UnitList
            
            self.__dict__['Width']    = 1*IN  ; UnitList['Width']  = Length
            self.__dict__['WidthBot']    = 1*IN  ; UnitList['WidthBot']  = Length
            self.__dict__['WidthTop']    = 1*IN  ; UnitList['WidthTop']  = Length
            self.__dict__['Height']   = 1*IN  ; UnitList['Height'] = Length
            self.__dict__['Material'] = ACMaterial()
            
            #
            # Set it to have no weight by default
            #
            self.Material.AreaForceDensity = 0*LBF/IN**2
            
            #
            # Create a component drawing so the bulk head can be drawn
            #
            self.Drawings.BulkHead = ACComponent.ComponentDrawing()
            
#==============================================================================
        def Corners(self):
            """
            Returns list of coordinates of the bulkhead
            """
            Width = self.Width #Standard Rectangular bulkhead
            WidthBot = self.WidthBot #Trapezoidal bulkhead only
            WidthTop = self.WidthTop #Trapezoidal bulkhead only
            Height = self.Height
            X = self.X #Center point (origin) of Bulkhead (see ACBase section 'ACComponent' line 521)
            if Width==-1*IN: #Trapezoidal bulkhead only
                pt1 = npy.array([X[0],X[1]-WidthBot/2,X[2]-Height/2])
                pt2 = npy.array([X[0],X[1]+WidthBot/2,X[2]-Height/2])
                pt3 = npy.array([X[0],X[1]+WidthTop/2,X[2]+Height/2])
                pt4 = npy.array([X[0],X[1]-WidthTop/2,X[2]+Height/2])
            else: #Standard Rectangular bulkhead
                pt1 = npy.array([X[0],X[1]-Width/2,X[2]-Height/2])
                pt2 = npy.array([X[0],X[1]+Width/2,X[2]-Height/2])
                pt3 = npy.array([X[0],X[1]+Width/2,X[2]+Height/2])
                pt4 = npy.array([X[0],X[1]-Width/2,X[2]+Height/2])             
            return (pt1,pt2,pt3,pt4)
#===============================================================================
        def Top(self):
            """
            Returns the top of the bulkhead
            """
            
            return self.X[2] + self.Height/2 #(Center point + 1/2 Height)

#===============================================================================
        def Bottom(self):
            """
            Returns the bottom of the bulkhead
            """
            
            return self.X[2] - self.Height/2 #(Center point - 1/2 Height)
        
#===============================================================================
        def Area(self):
            """
            Returns the area of the bulkhead
            """
            Width = self.Width #Standard Rectangular bulkhead
            WidthBot = self.WidthBot #Trapezoidal bulkhead only
            WidthTop = self.WidthTop #Trapezoidal bulkhead only
            Height = self.Height
            X = self.X #Center point (origin) of Bulkhead (see ACBase section 'ACComponent' line 521)
            if Width==-1*IN: #Trapezoidal bulkhead only
                return self.Height * ((self.WidthBot+self.WidthTop)/2) #Trapezoidal Area Formula
            else: #default square bulkhead
                return self.Width * self.Height #Basic Rectangle Area Formula

#===============================================================================
        def CG(self):
            """
            Returns the center of gravity for the bulkhead
            """
            X = self.GetX()
            Width = self.Width #Standard Rectangular bulkhead
            WidthBot = self.WidthBot #Trapezoidal bulkhead only
            WidthTop = self.WidthTop #Trapezoidal bulkhead only
            Height = self.Height
            
            if Width==-1*IN: #Trapezoidal bulkhead only
                return (X - Height/2)+((WidthBot+2*WidthTop)/(3*(WidthTop+WidthBot))*Height) #(Center point - 1/2 Height, a.k.a. "Bottom")+equation for y_bar (see http://mathworld.wolfram.com/Trapezoid.html)
            else: #default square bulkhead                
                return X

#===============================================================================
        def MOI(self):
            """
            Returns the moments of inertia about its CG as (Ixx, Iyy, Izz)
            Assumes X always nose to tail of aircraft and bulkhead is thin rectangular plate
            
            (**Or a thin-walled Trapezoidal plate! -Caleb)
            
            """
            if self.dirty:
                self.Refresh()

            CG = self.CG()
            Area=self.Area
            width = self.Width
            height = self.Height
            WidthBot = self.WidthBot #Trapezoidal bulkhead only
            WidthTop = self.WidthTop #Trapezoidal bulkhead only
            mass = self.Weight / g
            if width==-1*IN: #Trapezoidal bulkhead only
                """
                Note from Caleb: 
                See this discussion for derivation of flat rectangular plate:
                https://www.physicsforums.com/threads/moment-of-inertia-for-a-rectangular-plate.12903/
                It appears that the equation for the thin plate assumes mass=height*width=area
                (that is, removing density from mass=area*density for a theoretical flat plate)
                to modify the general equations for a rectangle, as seen here:
                http://www.efunda.com/math/areas/rectangle.cfm
                
                I have taken the same approach to modifying the equations for trapezoids from this site:
                http://www.efunda.com/math/areas/IsosTrapezoid.cfm
                (where for the given X=nose to tail, Z=vertical axis system, Ix -> Iy, Iy -> Iz, and Jz -> Ix))
                The equations are divided by the bulkhead area, and multiplied by mass.
                """
                
                Ix = 1/48 * mass * height*(4*(height**2)*WidthBot+12*(height**2)*WidthTop+7*(WidthBot**3)+7*WidthTop*(WidthBot**2)+WidthBot*(WidthTop**2)+(WidthTop**3)) / (self.Height * ((self.WidthBot+self.WidthTop)/2))
                Iy = 1/12 * mass * (height**3)*(3*WidthTop+WidthBot) / (self.Height * ((self.WidthBot+self.WidthTop)/2))
                Iz = 1/48 * mass * height*(WidthTop+WidthBot)*((WidthTop**2)+7*(WidthBot**2)) / (self.Height * ((self.WidthBot+self.WidthTop)/2))
            else: #default square bulkhead  
                #Beer-Johnston
                Ix = 1/12 * mass * (width**2 + height**2)
                Iy = 1/12 * mass * height**2
                Iz = 1/12 * mass * width**2

            return self.ToUnumList([Ix, Iy, Iz], LBM*IN**2)

#===============================================================================
        def GetWeight(self):
            """
            Retrieves Weight
            """
            if self.dirty:
                self.Refresh()
                
            Weight = self._CalcWeight
                
            return Weight
        
#===============================================================================
        def _CalcWeight(self):
            """
            Calculates the weight of the bulk head
            """
            Width = self.Width #Standard Rectangular bulkhead            
            if Width==-1*IN: #Trapezoidal bulkhead only
                return self.Height * ((self.WidthBot+self.WidthTop)/2) * self.Material.AreaForceDensity
            else: #default square bulkhead  
                #
                # Simply Area * AreaForceDensity
                #
                return self.Width * self.Height * self.Material.AreaForceDensity
        
               
#===============================================================================
        def Refresh(self):
            """
            Virtual function used to recompute attributes
    
            All derived instances of this methods must set
            self.param.refreshing = True at the beginning of the method and
            self.param.refreshing = False when completed.
    
            Otherwise an infinite loop will be created!!!
            """
            self.param.refreshing = True
            super(ACFuselage.BulkHead, self).Refresh()
            #
            # Perform all calculations
            #
            BulkHead = self.Drawings.BulkHead
            Width = self.Width / IN
            WidthTop = self.WidthTop / IN
            WidthBot = self.WidthBot / IN
            Height = self.Height / IN
            X = self.GetX() / IN
            Weight = self._CalcWeight()
            
            
            #
            # Top
            #
            if Width==-1: #Trapezoidal bulkhead only
                BulkHead.xt = npy.array([0       , 0      ]) + X[0]
                BulkHead.yt = npy.array([-WidthTop/2, WidthTop/2]) + X[1]    
            else: #default square bulkhead 
                BulkHead.xt = npy.array([0       , 0      ]) + X[0]
                BulkHead.yt = npy.array([-Width/2, Width/2]) + X[1]
        
            #
            # Front
            #
            if Width==-1: #Trapezoidal bulkhead only
                BulkHead.yf = npy.array([-WidthBot/2, WidthBot/2, WidthTop/2, -WidthTop/2, -WidthBot/2]) + X[1]
                BulkHead.zf = npy.array([-Height/2, -Height/2, Height/2,  Height/2, -Height/2]) + X[2]  
            else: #default square bulkhead 
                BulkHead.yf = npy.array([-Width/2, Width/2, Width/2, -Width/2, -Width/2]) + X[1]
                BulkHead.zf = npy.array([-Height/2, -Height/2, Height/2,  Height/2, -Height/2]) + X[2]
        
            #
            # Side
            #        Note: Works for Rectangle or Trapezoid Bulkhead.
            BulkHead.xs = npy.array([0, 0     ]) + X[0]
            BulkHead.zs = npy.array([-Height/2 , Height/2]) + X[2]
    
            self.param.refreshing = False
        
#===============================================================================
    class Section(ACStructural):
        """
        A section of the fuselage. The section always exists between two bulkheads.
        
        Attributes:
            Length      - The length of the section
            Align       - The relative alignment of the back bulkhead to the front bulkhead
            FrontBulk   - The front bulk head
            BackBulk    - The back bulk head
            SkinMat     - The 'skin' material of the section (i.e the area density)
            StringerMat - The 'stringer' material of the section (i.e. linear density)
                          The total weight is the sum of the weights due to
                          linear and area density
        """
        
#==============================================================================
        class Face():
            def __init__(self,pt1,pt2,pt3,pt4,sign,coord):
                """
                Inputs:
                    All points are counterclockwise around the box
                    pt1 - bottom left 
                    pt2 - bottom right
                    pt3 - top right
                    pt4 - top left
                    
                    The sign assures that normal vector always faces inward on the section
                    
                        4--3
                        |  |
                        1--2
                        
                    (Note: This bit confused me at first, so I've added a few more diagrams below here. -Caleb)
                """
                self.pt1 = pt1
                self.pt2 = pt2
                self.pt3 = pt3
                self.pt4 = pt4
                
                v1 = pt3 - pt1
                v2 = pt4 - pt2
                
                Normal = pyl.cross(v1,v2) * sign
                """
                Result: (Ascii art is weird in italics, just pretend 1 through 4 are in a square)
                          4     3
                           \  /
                            X      
                          /  \     
                        1     2
                                    The "X" where v1 and v2 cross is the "Normal" pointing out of the screen.
                """
                Normal = npy.array([Normal[i] / IN**2 for i in range(3)])
                norm = sum(Normal**2)**0.5
                
                self.coord = coord
                self.Area = norm / 2 * IN**2
                self.Normal = Normal / norm

#==============================================================================                
            def PointLocation(self,FractPos):
                """
                3D Bilinear Interpolation of point on a Face
                Input:
                    FractPos - tuple of (xi, eta, zeta) signifying percent location in axis
                    
                    TODO: Make ASCII picture explaining what the hell is going on here!
                    
                """
                coord = self.coord
                
                xi  = FractPos[coord[0]]
                eta = FractPos[coord[1]]
                
                v12 = self.pt2 - self.pt1
                v43 = self.pt3 - self.pt4

                
                p12 = xi * v12 + self.pt1
                p43 = xi * v43 + self.pt4
                
                veta = p43 - p12
         
                return eta * veta + p12
                """
                Result:
                            4....P43.->.3
                           |     ^     |    
                          |     |     |   
                         |     |     | 
                        1....P12.->.2    
                        
                        Where Veta is the center vertical vector moving from P12 to P43, 
                        as P12 and P43 move from left to right along vector 12 and vector 43.
                """       
#===============================================================================
        class Component(ACMassBox):
            """
            A Component to be placed within the fuselage.
            
            Attributes:
                Position    - Relative placement of component within the section
                Face        - The Face the component is attached to (None if positioned in 3D space)
            """
#===============================================================================            
            def __init__(self, Section):
                """
                Initializes the class
                
                Inputs:
                    Section    - The Section in which the Component is located.
                """
                
                super(ACFuselage.Section.Component, self).__init__()
                UnitList = self.UnitList
                
                self.__dict__['Position'] = (0,0,0);  UnitList['Position'] = Unitless
                self.__dict__['Face']     = None
                
                self.param.Section = Section
                
#===============================================================================                
            def Refresh(self):
                """
                Virtual function used to recompute attributes
        
                All derived instances of this methods must set
                self.param.refreshing = True at the beginning of the method and
                self.param.refreshing = False when completed.
        
                Otherwise an infinite loop will be created!!!
                """
                self.param.refreshing = True
                
                # Find spatial position
                if self.Face is not None:
                    MyFace    = self.param.Section.param.Faces.__dict__[self.Face]
                    self.X    = MyFace.PointLocation(self.Position)
                    self.Axis = MyFace.Normal
                else:
                    self.X = self.param.Section.PointLocation(self.Position)
                
                self.param.refreshing = False
                #
                # Call this after X and Axis have been updated
                # 
                super(ACFuselage.Section.Component,self).Refresh()
                
#==============================================================================                
        def __init__(self, FrontBulk, BackBulk, index, Fuselage):
            """
            Initializes the class.
            
            Inputs:
                FrontBulk - The front bulk head
                BackBulk  - The back bulk head
                index     - An index for sorting sections
                Fuselage  - The ACFuselage instance that owns this section
                L         - The length of the section
                A         - The relative alignment of the back bulkhead to the front bulkhead
            """
            super(ACFuselage.Section, self).__init__()
            UnitList = self.UnitList
                        
            self.__dict__['FrontBulk']   = FrontBulk
            self.__dict__['BackBulk']    = BackBulk
            self.NoneList['Length']      = None         ; UnitList['Length'] = Length
            self.NoneList['Align']       = None         ; UnitList['Align'] = Unitless
            self.__dict__['SkinMat']     = ACMaterial()
            self.__dict__['StringerMat'] = ACMaterial()
            
            #
            # Just set the AreaForceDensity to zero so the section weighs nothing
            #
            self.SkinMat.AreaForceDensity       = 0*LBF/IN**2
            self.StringerMat.LinearForceDensity = 0*LBF/IN
            
            #
            # This should never be changed. It is used to sort the sections.
            #
            self.param.index = index
            self.param.Fuselage = Fuselage
            
            #
            # Set the fuselage as the parent of the section
            #
            self.Parents.append( Fuselage )
            
            #
            # Make sure the section is dirtied if a bulkhead is dirtied
            #
            self.FrontBulk.Parents.append( self )
            self.BackBulk.Parents.append( self )
            
            class Faces:
                pass
            self.param.Faces = Faces() 
            
            self.Drawings.Edge1 = ACComponent.ComponentDrawing()
            self.Drawings.Edge2 = ACComponent.ComponentDrawing()
            
            self.param.Components = []
            
            #
            # A section is defined by the X locations of it's bulkheads
            #
            del self.__dict__['X']
            
#===============================================================================
        def AddComponent(self, ComponentName, Weight, LWH, Face = None, Position = (0.5,0.5,0.5)):
            """
            Adds a component to the fuselage
            
            Inputs:
                ComponentName - The name of the component
                Weight        - Weight of the component
                Position      - Tuple containing (xi, eta, zeta) percent locations for placement
                Face          - String Variable where None implies placed within section, otherwise 
                                Top, Bottom, Left, Right, Front, Back 
                                if attached to face of the section.
            """
            Components = self.param.Components
            
            #
            # Create the component in the sections class so it can be referenced later
            #
            Components.append(ACFuselage.Section.Component(self))
            Components[-1].Face     = Face
            Components[-1].Position = Position
            Components[-1].Weight   = Weight
            Components[-1].LWH      = LWH
            Components[-1].name     = ComponentName
            Components[-1].Parents.append( self )
            
            #
            # Initial guess for ease of placement
            #   Places location point in center of box if floating
            #   In center on face if attached to another face.
            #
            if Face is not None:
                Components[-1].Xat = [0, 0.5, 0.5]
            else:
                Components[-1].Xat = [0.5, 0.5, 0.5]
            
            #
            # Save of a reference to the component to make it so the user can type
            # Fus = ACFuselage()
            # Fus.AddComponent('Battery',)
            # Fus.Battery
            #
            self.__dict__[ComponentName] = Components[-1]

#===============================================================================
        def PointLocation(self,FractPos):
            """
            3D Trilinear Interpolation of point within a Section
            
            Inputs:
                FractPos - Tuple of (Xi,Eta,Zeta) containing fractional
                           position (0 to 1) within section
            """
            xi = FractPos[0]
            eta = FractPos[1]
            zeta = FractPos[2]
            
            fp = self.FrontBulk.Corners()
            bp = self.BackBulk.Corners()
            
            # Cube coordinates
            """
                p1 thru p4 are front bulkhead pts 1 thru 4
                p5 thru p8 are back bulkhead pts 1 thru 4
            """
            p1 = fp[0]
            p2 = fp[1]
            p3 = fp[2]
            p4 = fp[3]
            p5 = bp[0]
            p6 = bp[1]
            p7 = bp[2]
            p8 = bp[3]
            
            # Interpolate Xi plane corners in Section
            """
                p15 is Xi point along vector between p1 and p5
                p26 is Xi point along vector between p2 and p6
                p37 is Xi point along vector between p3 and p7
                p48 is Xi point along vector between p4 and p8
            """
            p15 = (p5 - p1) * xi + p1
            p26 = (p6 - p2) * xi + p2
            p37 = (p7 - p3) * xi + p3
            p48 = (p8 - p4) * xi + p4
            
            """                                8-----------7
                            4-----------3     |           |   
                           |           | ->  |  BackBulk |  
                          | FrontBulk |     |           |   
                         |           |     5-----------6 
                        1-----------2     
                        lines of interpolated points connect respective corners between front and back faces.
            """
            
            # Interpolate Eta line endpoints on Xi plane
            """
                p1526 is Eta point along vector between p15 and p26
                p4837 is Eta point along vector between p48 and p37
            """
            p1526 = (p26 - p15) * eta + p15
            p4837 = (p37 - p48) * eta + p48
            """
            Result: Interpolated points created from side to side on top and bottom faces of sections
            """
            # Interpolate Zeta point along Eta Line and return
            return (p4837 - p1526) * zeta + p1526
            """
            Result: Interpolated points created bottom to top
            """
#===============================================================================
        def SWet(self):
            """
            Returns the wetted area for the section. This is primarily used for drag
            and weight calculations. All four sides are used even though all four 
            sides of the section may not be exposed. This simply gives a 
            more conservative estimate.
            """
            if self.dirty or self.param.Fuselage.dirty: self.param.Fuselage.Refresh()
            
            SWet = 0 * IN**2
            
            Faces = self.param.Faces
            
            SWet += Faces.Top.Area
            SWet += Faces.Bottom.Area
            SWet += Faces.Left.Area
            SWet += Faces.Right.Area
                
            return SWet

#===============================================================================
        def Top(self, x):
            """
            Calculates the top of the section for the given x location
            
            Inputs:
                x - An global x location between the front and back bulkhead
            """
            if self.dirty or self.param.Fuselage.dirty: self.param.Fuselage.Refresh()
            
            FX = self.FrontBulk.X[0]
            BX = self.BackBulk.X[0]
                    
            if FX > x or BX < x:
                return None #A flag to indicate this is not the correct section for this X

            FT = self.FrontBulk.Top()
            BT = self.BackBulk.Top()
            
            #
            # Linearly interpolate the top location between the two bulkheads
            #
            return FT + (BT-FT)/(BX-FX)*(x-FX)

#===============================================================================
        def Bottom(self, x):
            """
            Calculates the bottom of the section for the given x location
            
            Inputs:
                x - An global x location between the front and back bulkhead
            """
            if self.dirty or self.param.Fuselage.dirty: self.param.Fuselage.Refresh()
            
            FX = self.FrontBulk.X[0]
            BX = self.BackBulk.X[0]
                    
            if FX > x or BX < x:
                return None #A flag to indicate this is not the correct section for this X

            FB = self.FrontBulk.Bottom()
            BB = self.BackBulk.Bottom()
            
            #
            # Linearly interpolate the top location between the two bulkheads
            #
            return FB + (BB-FB)/(BX-FX)*(x-FX)

#===============================================================================
        def Height(self, x):
            """
            Calculates the height of the section for the given x location
            
            Inputs:
                x - A global x location between the front and back bulkhead
            """
            if self.dirty or self.param.Fuselage.dirty: self.param.Fuselage.Refresh()

            return self.Top(x) - self.Bottom(x)

#===============================================================================
        def CG(self):
            """
            CG is simply half way between the bulk heads
            """
            if self.dirty or self.param.Fuselage.dirty: self.param.Fuselage.Refresh()
            
            Components = self.param.Components
            
            TempCG = npy.array([0,0,0])*IN*LBF
            TempWeight = 0*LBF
            
            #Loop over Components within section 
            for Component in Components:
                CompWeight = Component.Weight
                CompCG = Component.CG()

                TempCG += CompCG*CompWeight
                TempWeight += CompWeight
            
            #Get Section structural weight
            if self.NoneList.has_key("Weight"):
                SWeight = self._CalcSectionWeight()
            else:
                SWeight = self.__dict__["Weight"]
            
            #
            #Calculate Section CG
            #
            BBCG = self.BackBulk.CG()
            FBCG = self.FrontBulk.CG()
            FArea = self.param.Faces.Front.Area
            BArea = self.param.Faces.Back.Area
            SectCG = (self.BackBulk.CG()*FArea + self.FrontBulk.CG()*BArea)/(2*(FArea + BArea))
            
            TempCG += SectCG*SWeight
            TempWeight += SWeight
            
            #Avoid devide by zero case in case of no components or undefined section weight 
            if (TempWeight)==0*LBF:  
                CG = SectCG #use rough CG between bulkheads
            else:          
                CG = TempCG/TempWeight
                
            return CG

#===============================================================================
        def MOI(self):
            """
            Sums Moments of Inertia of all components within section about CG
            Parallel Axis Theorem used for components, and skin and bulkheads not taken into account.
            """
            if self.dirty or self.param.Fuselage.dirty: self.param.Fuselage.Refresh()
            
            Components = self.param.Components
            
            CG = self.CG()
            TempMOI = npy.array([0,0,0])*LBM*IN**2
            
            for Component in Components:
                CompCG = Component.CG()
                CompMOI = Component.MOI()            
                
                CompM = Component.Weight / g
                dCG = CompCG - CG
                
                #Add current Component MOI to the TempMOI calc
                TempMOI[0] += CompMOI[0] + CompM*(dCG[1]**2 + dCG[2]**2)#Ixx
                TempMOI[1] += CompMOI[1] + CompM*(dCG[0]**2 + dCG[2]**2)#Iyy
                TempMOI[2] += CompMOI[2] + CompM*(dCG[1]**2 + dCG[0]**2)#Izz
            
            #TODO: Find Moments of inertia component for skin...
           
            return TempMOI

#===============================================================================
        def _CalcComponentWeight(self):
            """
            Calculates the weight of the components within section. 
            """
            
            Weight = 0 * LBF
            for Component in self.param.Components:
                Weight += Component.Weight
            
            return Weight
        
#===============================================================================
        def _CalcSectionWeight(self):
            """
            Calculates the weight of the section. 
            
            This is conservative as it assumes the
            skin wraps around all 4 sides of the section.
            """
            
            Weight  = self.SWet() * self.SkinMat.AreaForceDensity
            Weight += self.Length * self.StringerMat.LinearForceDensity

            return Weight

#===============================================================================
        def _CalcWeight(self):
            """
            Sums Weight of section
            """
            
            Weight = self._CalcSectionWeight()
            Weight += self._CalcComponentWeight()
            
            #
            # Simply Wetted Area * AreaForceDensity
            #
            return Weight

#===============================================================================
        def _CalcLength(self):
            """
            Calculates the length of the section.
            """
            #
            # Simply the distance between the bulkheads
            #
            return self.BackBulk.X[0] - self.FrontBulk.X[0]

#===============================================================================
        def __getattr__(self, key):
            """
            Attempts to calculate any variables set to None
            """
            ans = super(ACFuselage.Section, self).__getattr__(key)
            if ans is not None:
                return ans
    
            # Try to calculate the value instead
            if key == 'Length':
                return self._CalcLength()
            #
            # If no calculation exist return None
            #
            return None

#===============================================================================
        def __getattribute__(self, key):
            """
            Used to hand the case where a section weight is specified.
            """
            ans = super(ACFuselage.Section, self).__getattribute__(key)
            if key == 'Weight' and ans is not None:
                return ans + ACFuselage.Section._CalcComponentWeight(self)
            
            return super(ACFuselage.Section, self).__getattribute__(key)
            
#===============================================================================
        def Draw(self, fig = 1, top = 221, side = 222, front = 223):
            """
            Draws the section and its components
    
            Inputs:
                fig   - Integer number of figure to draw into
                top   - Subplot for top view
                side  - Subplot for side view
                front - Subplot for front view
            """
            if self.dirty or self.param.Fuselage.dirty: self.param.Fuselage.Refresh()
            
            super(ACFuselage.Section, self).Draw(fig, top, side, front)
            
            Components = self.param.Components
    
            for Component in Components:
                Component.Draw(fig, top, side, front)

#===============================================================================
        def AddToWeightTable(self, PartName, WeightTable):
            """
            Adds the components and its parts to the Weight Table
            
            Input:
                PartName    - The name give to this part
                WeightTable - The table to add self to 
            """
            if self.dirty: self.Refresh()
            
            WeightTable[PartName] = self
            
            Components  = self.param.Components
            
            for Component in Components:
                Component.AddToWeightTable(Component.name, WeightTable[PartName])

#===============================================================================
        def Refresh(self):
            """
            Virtual function used to recompute attributes
    
            All derived instances of this methods must set
            self.param.refreshing = True at the beginning of the method and
            self.param.refreshing = False when completed.
    
            Otherwise an infinite loop will be created!!!
            """
            super(ACFuselage.Section, self).Refresh()
            self.param.refreshing = True
            #
            # Perform all calculations
            #
            FrontBulk = self.FrontBulk
            BackBulk  = self.BackBulk
            
            #
            # Position the Back bulkhead Length behind the Front Bulkhead
            #
            if self.Length is not None:
                BackBulk.X[0] = FrontBulk.X[0] + self.Length
                
            if self.Align is not None:
                FtBlkAlignPt = FrontBulk.CG()[2] + FrontBulk.Height/2*self.Align
                BkBlkAlignPt = BackBulk.CG()[2] + BackBulk.Height/2*self.Align
                if (BkBlkAlignPt >= BackBulk.CG()[2] + BackBulk.Height/2):
                    BkBlkAlignPt = BackBulk.CG()[2] + BackBulk.Height/2
                elif (BkBlkAlignPt <= BackBulk.CG()[2] - BackBulk.Height/2):
                    BkBlkAlignPt = BackBulk.CG()[2] - BackBulk.Height/2
                BackBulk.X[2] += FtBlkAlignPt - BkBlkAlignPt

                
            BackBulk.Refresh()
            
            #
            # Create Faces for internal placement within section
            #
            frntpts = FrontBulk.Corners()
            backpts = BackBulk.Corners()
            Faces = self.param.Faces
            Face = ACFuselage.Section.Face
            
            x=0;y=1;z=2
            
            Faces.Front  = Face(frntpts[0],frntpts[1],frntpts[2],frntpts[3], 1,(y,z))
            Faces.Back   = Face(backpts[0],backpts[1],backpts[2],backpts[3],-1,(y,z))
            Faces.Left   = Face(frntpts[0],backpts[0],backpts[3],frntpts[3],-1,(x,z))
            Faces.Right  = Face(frntpts[1],backpts[1],backpts[2],frntpts[2], 1,(x,z))
            Faces.Top    = Face(frntpts[3],backpts[3],backpts[2],frntpts[2],-1,(x,y))
            Faces.Bottom = Face(frntpts[0],backpts[0],backpts[1],frntpts[1], 1,(x,y))
 
 
            #
            # Refresh the components
            #
            Components = self.param.Components
    
            for Component in Components:
                Component.Refresh()

            #
            # Draw the lines to distinguish the section
            #
            Edge1 = self.Drawings.Edge1
            Edge2 = self.Drawings.Edge2
            
            
            FW = FrontBulk.Width / IN
            FH = FrontBulk.Height / IN
            
            BW = BackBulk.Width / IN
            BH = BackBulk.Height / IN
            
            FX = FrontBulk.X
            BX = BackBulk.X
            
            FX = [FX[0] / IN, FX[1] / IN, FX[2] / IN]
            BX = [BX[0] / IN, BX[1] / IN, BX[2] / IN]
            
            #
            # Top
            #
            Edge1.xt = npy.array([FX[0]       , BX[0]      ]) 
            Edge1.yt = npy.array([FX[1] + FW/2, BX[1] + BW/2])

            Edge2.xt = npy.array([FX[0]       , BX[0]      ]) 
            Edge2.yt = npy.array([FX[1] - FW/2, BX[1] - BW/2])
        
            #
            # Sections cannot bee seen from the Front
            #
        
            #
            # Side
            #
            Edge1.xs = npy.array([FX[0]       , BX[0]       ]) 
            Edge1.zs = npy.array([FX[2] + FH/2, BX[2] + BH/2])

            Edge2.xs = npy.array([FX[0]       , BX[0]       ]) 
            Edge2.zs = npy.array([FX[2] - FH/2, BX[2] - BH/2])
    
            self.param.refreshing = False

#===============================================================================
    class ACPayload(Section.Component):
        """
        A payload class
        
        Attributes:
            Length   - The length of the payload
            Width    - The width of the payload
            Weight   - The Weight of the payload
            Material - The material of the payload
        """
        def __init__(self):
            super(ACFuselage.ACPayload, self).__init__(None)
            UnitList = self.UnitList
            
            self.__dict__['Length']     = 1*IN  ; UnitList['Length'] = Length
            self.__dict__['Width']      = 1*IN  ; UnitList['Width']  = Length
            self.__dict__['Weight']     = 0*LBF ; UnitList['Weight'] = Force
            self.__dict__['XcgSecFrac'] = 0.5   ; UnitList['XcgSecFrac'] = Unitless
            self.__dict__['Material'] = ACMaterial()
            
            #Set the CG in the middle to start
            self.Position = (self.XcgSecFrac,0.5,0)
            
            #
            # Set it to have no weight by default
            #
            self.Material.ForceDensity = 1*LBF/IN**3
            
#==============================================================================
        def Refresh(self):
            """
            Update drawing array
    
            Virtual function used to recompute attributes
    
            All derived instances of this methods must set
            self.param.refreshing = True at the beginning of the method and
            self.param.refreshing = False when completed.
    
            Otherwise an infinite loop will be created!!!
            """
            self.param.refreshing = True
            #
            # Perform all calculations
            #
            
            #Length = self.Length
            #Width  = self.Width
            Weight = self.Weight
            
            #Move the Payload to the desired location
            self.Position = (self.XcgSecFrac,0.5,0)
            
            Length = 0.95*self.param.Section.Length*min( abs(1 - self.XcgSecFrac), abs( self.XcgSecFrac ) )*2
            Width  = 0.95*min( self.param.Section.FrontBulk.Width, self.param.Section.BackBulk.Width )

            Volume = Weight / self.Material.ForceDensity
            
            self.LWH = [Volume / (Width*Length), Width, Length ]
            
            self.param.refreshing = False
            #
            # Need to refresh on at the end as we are modifying LWH
            #
            super(ACFuselage.ACPayload, self).Refresh()
            
#===============================================================================
# Define comparison operators for the Section so sections can be sorted 
#===============================================================================
#        def __gt__(self, otherSection):
#            return self.param.index > otherSection.param.index
#        def __lt__(self, otherSection):
#            return self.param.index < otherSection.param.index

#===============================================================================
# ACFuselage
#===============================================================================
    def __init__(self):
        super(ACFuselage, self).__init__()
        
        UnitList = self.UnitList
        
        self.name = "Fuselage"
        
        self.__dict__['XOffset']    = 0*IN  ; UnitList['XOffset']    = Length
        self.NoneList['XcgSection'] = None
        self.NoneList['TailBulk']   = None
        self.__dict__['XcgSecFrac'] = 0.5   ; UnitList['XcgSecFrac'] = Unitless
        self.__dict__['C_dc']       = 1.7   ; UnitList['C_dc']       = Unitless
        
        self.__dict__['Payload'] = ACFuselage.ACPayload()
        self.Payload.Xat = [0, 0.5, 0.5]
        self.Payload.Color = 'magenta'
        
        #
        # The sections of the fuselage
        #
        self.__dict__['Sections'] = []
                
#===============================================================================
    def AddSection(self, SectionName, Length = None, Align = None):
        """
        Adds a section to the fuselage
        
        Inputs:
            SectionName - The name of the section
            Length      - The length of the section
            Align       - The Alignment of the section
                          Note: This value can be any positive or negative decimal, 
                                but behaviour is somewhat arbitrary.  
                                Following values hold true at all times...
                                  0.0  = align by centerline
                                  1.0  = align by top of bulkheads
                                  -1.0 = align by bottom of bulkheads
                                  None = centers bulkhead on CG Z location
        """
        Sections = self.Sections
        
        if len(Sections) == 0:
            #
            # This is the first section created, so the front bulk head must be constructed
            #
            FrontBulk = ACFuselage.BulkHead()
        else:
            #
            # Appending a section, so the last sections BackBulk is the new sections FrontBulk
            #
            FrontBulk = Sections[-1].BackBulk
            
        BackBulk = ACFuselage.BulkHead()
        index = len(Sections)
        
        #
        # Create the section in the sections class so it can be referenced later
        #
        Sections.append(ACFuselage.Section(FrontBulk, BackBulk, index, self))

        Sections[-1].Length = Length
        Sections[-1].Align  = Align
        Sections[-1].name   = SectionName
        
        #
        # Save of a reference to the section to make so the user can type
        # Fus = ACFuselage()
        # Fus.AddSection('Nose', 10*IN)
        # Fus.Nose
        #
        self.__dict__[SectionName] = Sections[-1]
        
        if len(Sections) == 1:
            self.XcgSection = Sections[-1]

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the fuselage and its sections to the Weight Table
        
        Input:
            PartName    - The name give to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WeightTable[PartName] = self
        
        Sections = self.Sections
        
        for Section in Sections:
            WeightTable[PartName][Section.name + 'FrontBulk'] = Section.FrontBulk
            Section.AddToWeightTable(Section.name, WeightTable[PartName])

        WeightTable[PartName][Section.name + 'BackBulk'] = Section.BackBulk
                            
#===============================================================================
    def SWet(self):
        """
        Sums up the wetted area of all sections
        """
        if self.dirty: self.Refresh()
        
        SWet = 0*IN**2
        for Section in self.Sections:
            SWet += Section.SWet()
        return SWet

#===============================================================================
    def _CalcWeight(self):
        """
        Calculates the weight of the fuselage as the sum of the weight of
        the bulk heads and sections.
        """
        if self.dirty: self.Refresh()
        
        Sections = self.Sections
        Weight = Sections[0].FrontBulk.Weight        
        for Section in Sections:
            Weight += Section.Weight
            Weight += Section.BackBulk.Weight #Only have to sum back bulk head because bulkheads are shared between sections
        
        return Weight

#===============================================================================
    def CG(self):
        """
        Calculates the CG for the Fuselage by summing over all sections
        """
        if self.dirty: self.Refresh()
        
        Sections = self.Sections
        
        TempCG = npy.array([0,0,0])*IN*LBF
        TempWeight = 0*LBF
        
        BulkCG      = Sections[0].FrontBulk.CG()
        BulkWeight  = Sections[0].FrontBulk.Weight
        TempCG     += BulkCG*BulkWeight
        TempWeight += BulkWeight
        
        for Section in Sections:
            SectCG      = Section.CG()
            SectWeight  = Section.Weight
            TempCG     += SectCG*SectWeight
            TempWeight += SectWeight
            
            BulkCG      = Section.BackBulk.CG()
            BulkWeight  = Section.BackBulk.Weight
            TempCG     += BulkCG*BulkWeight
            TempWeight += BulkWeight
            
        CG = TempCG/TempWeight
        
        return CG

#===============================================================================
    def MOI(self):
        """
        Calculates the Moments of Inertia for the Fuselage by summing over all sections
        via the parallel axis theorem
        """
        if self.dirty: self.Refresh()
        
        CG = self.CG()
        Sections = self.Sections
        
        TempMOI = npy.array([0,0,0])*LBM*IN**2
        
        #First bulkhead        
        BulkCG      = Sections[0].FrontBulk.CG()
        BulkM       = Sections[0].FrontBulk.Weight / g
        BulkMOI     = Sections[0].FrontBulk.MOI()
        dCG = BulkCG - CG
        
        TempMOI[0] += BulkMOI[0] + BulkM*(dCG[1]**2 + dCG[2]**2)#Ix
        TempMOI[1] += BulkMOI[1] + BulkM*(dCG[0]**2 + dCG[2]**2)#Iy
        TempMOI[2] += BulkMOI[2] + BulkM*(dCG[1]**2 + dCG[0]**2)#Iz
        
        for Section in Sections:
            SectCG      = Section.CG()
            SectM       = Section.Weight / g
            SectMOI     = Section.MOI()
            dCG         = SectCG - CG
            
            TempMOI[0] += SectMOI[0] + SectM*(dCG[1]**2 + dCG[2]**2)#Ix
            TempMOI[1] += SectMOI[1] + SectM*(dCG[0]**2 + dCG[2]**2)#Iy
            TempMOI[2] += SectMOI[2] + SectM*(dCG[1]**2 + dCG[0]**2)#Iz
                        
            BulkCG      = Section.BackBulk.CG()
            BulkM       = Section.BackBulk.Weight / g
            BulkMOI     = Section.BackBulk.MOI()
            dCG         = BulkCG - CG
                        
            TempMOI[0] += BulkMOI[0] + BulkM*(dCG[1]**2 + dCG[2]**2)#Ix
            TempMOI[1] += BulkMOI[1] + BulkM*(dCG[0]**2 + dCG[2]**2)#Iy
            TempMOI[2] += BulkMOI[2] + BulkM*(dCG[1]**2 + dCG[0]**2)#Iz
            
        return TempMOI

#===============================================================================
    def AircraftCG(self):
        """
        Returns the desired location of the CG for the aircraft
        """
        if self.dirty: self.Refresh()
        
        CG     =  self.XcgSection.FrontBulk.GetX()
        CG[0] += (self.XcgSection.Length*self.XcgSecFrac)
        
        return CG

#===============================================================================
    def Top(self, x):
        """
        Returns the top of the fuselage at the given x location
        
        Inputs:
            x - The global x location where the top of the fuselage is desired
        """
        if self.dirty: self.Refresh()
        
        Sections = self.Sections
        
        for Section in Sections:
            Top = Section.Top(x)
            if Top is not None:
                return Top
            
        return Sections[-1].BackBulk.Top()

#===============================================================================
    def Bottom(self, x):
        """
        Returns the bottom of the fuselage at the given x location
        
        Inputs:
            x - The global x location where the bottom of the fuselage is desired
        """
        if self.dirty: self.Refresh()
        
        Sections = self.Sections
        
        for Section in Sections:
            Bottom = Section.Bottom(x)
            if Bottom is not None:
                return Bottom

        return Sections[-1].BackBulk.Bottom()

#===============================================================================
    def Height(self, x):
        """
        Returns the height of the fuselage at the given x location
        
        Inputs:
            x - The global x location where the bottom of the fuselage is desired
        """
        if self.dirty: self.Refresh()
        
        return self.Top(x) - self.Bottom(x)

#===============================================================================
    def GetNoseBulk(self):
        """
        Returns the foremost bulkhead
        """
        if self.dirty: self.Refresh()
        
        return self.Sections[0].FrontBulk

#===============================================================================
    def GetNoseSection(self):
        """
        Returns the foremost section
        """
        if self.dirty: self.Refresh()
        
        return self.Sections[0]
    
#===============================================================================
    def _FrontalArea(self):
        """
        Returns the maximum cross-sectional area of the fuselage
        """
        if self.dirty: self.Refresh()
        
        Sections = self.Sections
        Area = Sections[0].FrontBulk.Area()
        
        for Section in Sections:
            Area = max(Area,Section.BackBulk.Area())
        
        return Area
    
#===============================================================================
    def _Length(self):
        """
        Returns the length of the fuselage
        """
        Sections = self.Sections
        return Sections[-1].BackBulk.X[0] - Sections[0].FrontBulk.X[0]

#===============================================================================
    def CD(self, alpha2d, V, h, Wing):
        """
        Returns the drag of the fuselage
        
        Inputs:
            alpha2d - 2D angle of attack of the main wing
            V       - Velocity
            h       - Altitude
            Wing    - A wing used to calculate fuselage angle of attack
        """
        if self.dirty: self.Refresh()
        
        Vs       = self.MakeList(V)
        AlphaFus = self.MakeList( Wing.AlphaFus(alpha2d) )

        Length   = self._Length()
        Swet     = self.SWet()
        Sections = self.Sections
        Area     = self._FrontalArea()
        Deff     = (4 * Area / math.pi)**0.5 #Effective diameter of fuselage


        #Calculates effective fuselage slenderness ratio
        Slender1 = Length / Deff
        Slender2 = (Sections[0].Length + Sections[-1].Length) / Deff + 2
        if Slender1 < Slender2:
            Slender = Slender1
        else: 
            Slender = Slender2

        #TODO: Check Slender2 equation in Torenbeek
   
        ShapeFactor = 2.2 / (Slender ** 1.5) + 3.8 / Slender ** 3 #Calculates fuselage shape factor

        if len(Vs) > 1 and len(AlphaFus) > 1:
            raise "Can only specify a range of angles of attack or velocities"
        
        DqFuse   = []
        for i in xrange(len(Vs)):
            #
            # The log will error if the velocity is zero
            #
            V = Vs[i] if Vs[i] > 0.1*FT/SEC else 0.1*FT/SEC
                
            ReFuse   = AeroUtil.Re(h, V, Length)
                        
            CfFuse = 0.455 / (abs(math.log(ReFuse))**2.58) #Fuselage Skin Coefficient
            
            DqFuse0 = CfFuse * Swet * (1 + ShapeFactor)
            
            for a in AlphaFus:
                sina = abs(math.sin(a / RAD)**3)
            
                DqFuse.append( (DqFuse0 + self.C_dc * Swet * 0.25 * sina)/FT**2 )
        
        return DqFuse[0]*FT**2 if len(DqFuse) == 1 else npy.array( DqFuse )*FT**2

#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """
        Draws this fuselage and its sections

        Inputs:
            fig   - Integer number of figure to draw into
            top   - Subplot for top view
            side  - Subplot for side view
            front - Subplot for front view
        """
        super(ACFuselage, self).Draw(fig, top, side, front)
        
        Sections = self.Sections

        Sections[0].FrontBulk.Draw(fig, top, side, front)
        for Section in Sections:
            Section.Draw(fig, top, side, front)
            Section.BackBulk.Draw(fig, top, side, front)

        if self.Payload:
            self.Payload.Draw(fig, top, side, front)
        
#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACFuselage, self).Refresh()
        
        self.param.refreshing = True
        #
        # Perform all calculations
        #

        #
        # Position the first bulkhead and update it.
        # The first section is always assumed the nose section.
        #
        Sections = self.Sections
        
        Nose = Sections[0]
        
        
        Nose.FrontBulk.X[0] = self.XOffset
        Nose.FrontBulk.Refresh()
                       
        #
        # Refresh all the sections. This will position all remaining bulkheads.
        #
        for Section in Sections:
            Section.Refresh()
        
 
        #
        # Update the payload
        #
        if self.Payload:
            #self.Payload.Position = (self.XcgSecFrac,0.5,0)
            self.Payload.param.Section = self.XcgSection
            self.Payload.Refresh()
                
       
        self.param.refreshing = False
    
################################################################################
if __name__ == '__main__':
#    
#    from TestPropulsion import Propulsion
    
    Fuselage = ACFuselage()
    
    Fuselage.AddSection('Nose'     , 7*IN,  1)
    Fuselage.AddSection('PyldBay'  , 10*IN, 1)
    Fuselage.AddSection('TailTaper', 15*IN)
    Fuselage.AddSection('Tail')

    #Fuselage.Nose.Length = 10*IN
    Fuselage.Nose.FrontBulk.Width = 3*IN
    Fuselage.Nose.FrontBulk.Height = 3*IN
    Fuselage.Nose.FrontBulk.Weight = 0.5*LBF
    Fuselage.Nose.Weight = 0.2*LBF
    
    Fuselage.PyldBay.FrontBulk.Width = 5*IN
    Fuselage.PyldBay.FrontBulk.Height = 5*IN
    Fuselage.PyldBay.SkinMat.AreaForceDensity = 0.01*LBF/IN**2

    Fuselage.PyldBay.BackBulk.Width = 5*IN
    Fuselage.PyldBay.BackBulk.Height = 5*IN
    
    Fuselage.TailTaper.Align = 1
    
    Fuselage.Tail.BackBulk.X[0] = 50*IN
    Fuselage.Tail.BackBulk.X[2] = 0*IN
    Fuselage.Tail.Weight = 0.7*LBF
    
    Fuselage.Nose.AddComponent    (     "Battery"   , 0.5*LBF , (0.25*IN,1.5*IN,1*IN)    , "Right"  , (0.2,0.5,0.5)   )
#    Fuselage.PyldBay.AddComponent (     "Payload"   , 15.*LBF , (5*IN,3.5*IN,3*IN)       ,  None    , (0.5, 0.5, 0.5)   )
    Fuselage.Nose.AddComponent    (     "FuelTank"  , 0.5*LBF , (2.5*IN,2*IN,1.25*IN)    , "Back"   , (0.75, 0.5, 0.7) )
    Fuselage.Nose.AddComponent    ("NoseWheelServo" , 0.1*LBF , (.5*IN,1*IN,1*IN)        , "Bottom" , (0.6, 0.2, 0) )
    Fuselage.Nose.AddComponent    ("Engine"         , 2.0*LBF , (-3.5*IN,2*IN,2.5*IN)      , "Front"  )
    
#    Fuselage.Nose.SetPropulsion(Propulsion, "Front")
    
    Fuselage.XcgSection = Fuselage.PyldBay
    
    Fuselage.XOffset = 3*IN

    print 'Payload Bay SWet:', Fuselage.PyldBay.SWet()
    print 'Fuselage weight :', Fuselage.Weight
    print 'Tail length     :', Fuselage.Tail.Length
    print 'Nose CG         :', Fuselage.Nose.CG()
    print 'Nose Weight     :', Fuselage.Nose.Weight
    print 'Payload Bay CG  :', Fuselage.PyldBay.CG()
    print 'Tail Taper CG   :', Fuselage.TailTaper.CG()
    print 'Tail CG         :', Fuselage.Tail.CG()
    print 'Fuselage CG     :', Fuselage.CG()
    print 'Aircraft CG     :', Fuselage.AircraftCG()
    print 'Fuselage MOI    :', AsUnit(Fuselage.MOI(), "slug*ft**2", "%3.5f")

    Fuselage.Draw()
    
    pyl.show()
