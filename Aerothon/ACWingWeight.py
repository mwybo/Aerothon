"""
A collection of classes for calculating the weight of a lifting surface
"""
from __future__ import division # let 5/2 = 2.5 rather than 2 
from ACMaterial import ACMaterial
from ACBase import ACBase, ACComponent, ACStructural
from ACBase import ForceDensity, Length, Force, Area, AreaMomentOfInertia, Pressure, Unitless, g
from ACMass import ACMassBox, ACMassCyl
from scalar.units import LBF, IN, FT, PSI, SEC, ARCDEG, LBM, SLUG
from scalar.units import AsUnit
import numpy as npy
import pylab as pyl
from scipy import integrate
import math

################################################################################
class ACWingWeight(ACStructural):
    """
    A base class for calculating the weight of a wing

    Attributes:
        SkinMat           - ACMaterial class defining the material of the skin
    """
#===============================================================================
    class Spar(ACStructural):
        """
            A Spar to be placed within the Wing.  This spar can be either a tubespar or a square spar
            
            Attributes:
                I        - Area moment of inertia for spar (Optional, will be calculated otherwise)
                Area     - Cross-sectional Area of Spar (square assumed)
                
                Position - Relative placement of component within the Wing (see Note 1)
                   Note 1:  Position is a nondimensional two value tuple
                            the first value is position along chord (LE = 0, TE = 1)
                            the second value is position within the thickness (-1 = bottom surface, 1 = top surface, 0 = camberline
                PctSpan  - Percent of span that this spar occurs within, assumes centered and defaults to 1.0 (full span)
                TubeSpar - Logical value (T/F) if spar is a tube
                Dim1     - Primary Dimension, typically Height   (see Note below)
                Dim2     - Secondary Dimension, typically Width  (see Note below)
                   Note 2:  If TubeSpar = True
                              Dim1 = Outer Diameter,  Dim2 = Inner Diameter
                            If TubeSpar = False
                              Dim1 = Height,          Dim2 = Width
                
                MaxBendStress - Maximum allowable bending stress of the spar 
                ScaleToWing - 
        """
#===============================================================================
        def __init__(self, Wing):
            """
            Initializes the class
            
            Inputs:
                Wing    - The Wing in which the Spar is located.
            """
            super(ACWingWeight.Spar, self).__init__()
            UnitList = self.UnitList

            self.__dict__['Position'] = (0,0);  UnitList['Position'] = Unitless
            self.NoneList['I']        = None ;  UnitList['I']        = AreaMomentOfInertia
            self.NoneList['Area']     = None ;  UnitList['Area']     = Area
            self.NoneList['Dim1']     = None ;  UnitList['Dim1']     = Length
            self.NoneList['Dim2']     = None ;  UnitList['Dim2']     = Length
            self.__dict__['PctSpan']  = 1.0  ;  UnitList['PctSpan']  = Unitless
            self.__dict__['SparMat']  = ACMaterial()
            self.__dict__['TubeSpar'] = False
            self.__dict__['DrawDetail']=False
            self.__dict__['Structural']=True
            self.__dict__['ScaleToWing']= [True,True]
            
            self.NoneList['MaxBendStress']  = None  ;   UnitList['MaxBendStress'] = Pressure
            
            self.SparMat.ForceDensity = 0 * LBF/IN**3  # Default weight = 0
            
            self.param.Wing = Wing
            
            #
            # Create the drawing for the spar
            #
            self.Drawings.SparDrawing = ACComponent.ComponentDrawing()
                        
#===============================================================================
        def CG(self):
            """
            Calculates the Center of Gravity of the spar assuming thin, long bar which means CG is at center
            """
            if self.dirty: self.Refresh()
                        
            Wing = self.param.Wing
            
            b = Wing.b
            b = b/2. if Wing.FullWing else b

            spartip = b*self.PctSpan
            
            y = 0*IN if Wing.FullWing else spartip/2
            
            CG = self._CalcCentroid(y)
            
            return CG
            
#            
#            #The width is Dim2 for a rectangle, and Dim1 for a tube
#            Width = self.Dim2 if self.TubeSpar else self.Dim1
#                 
#            
#            CG = npy.array([0, 0, 0]) * IN   
#            
#            CG[0] = Wing.LE(spartip) + Wing.Chord(spartip)*self.Position[0] + Width/2
#            CG[1] = (self.X[1] + spartip) / 2
#            CG[2] = Wing.Lower(spartip, x = CG[0]) + (Wing.Thickness(spartip, x = CG[0])*self.Position[1]) / 2
#            
#            if Wing.FullWing:
#                CG[1] = self.X[1]
#                           
#            if abs(Wing.Axis[1]) == 1:  #Vertical Wing
#                temp = CG[0]
#                CG[0] = CG[2]
#                CG[2] = temp
#            
#            return CG
            
#===============================================================================
        def MOI(self):
            """
            Calculates the moment of inertia of the spar
            """
            if self.dirty: self.Refresh()
            
            Wing = self.param.Wing
        
            MOI = npy.array([0,0,0])*LBM*IN**2
            
            #
            # Calculate the MOI of the spar
            #
            Weight = self.Weight
            
            #NOTE: Do correct thing for TubeSpar
            
            #
            # MOI of spar in Local coordinates (y is out wingtip)
            #
            MOI[0] = 1/12 * Weight * (Wing.b**2) / g
            MOI[1] = 0.0*LBM*IN**2
            MOI[2] = 1/12 * Weight * (Wing.b**2) / g

            if abs(Wing.Axis[1]) == 1:  #Vertical Wing
                temp = MOI[1]
                MOI[1] = MOI[2]
                MOI[2] = temp
            
            return MOI

#===============================================================================
        def BendingStress(self, y, V, Moments = None):
            """
            Calculates the bending stress of the spar
            
            Input:
                y - Spanwise location(s)
                V - Velocity
            """
            if self.dirty: self.Refresh()
            
            Wing = self.param.Wing
            
            if Moments is not None:
                M = Moments
            else:
                M = self.MomentLoading(y, V)
            
            h = []
            for ay in y:
                h.append(Wing.Thickness(ay, Wing.LE(ay) + 0.25*Wing.Chord(ay)) / IN)
            
            h = h[0]*IN if len(h) == 1 else npy.array(h)*IN
            I = self._CalcI(y)
            
            return M*h/(2.*I)

#===============================================================================
        def MaxBendingStress(self, V):
            """
            Calculates maximum the bending stress of the spar
            
            Input:
                V - Velocity
            """        
            return self.BendingStress(0*IN, V)

#===============================================================================
        def PlotBendingStress(self, V, fig = 1):
            """
            Plots the bending stress of the spar
    
            Inputs:
                V   - Velocity to compute loads
                fig - Figure Number
            """
            if self.dirty: self.Refresh()
            
            Wing = self.param.Wing
            b = Wing.b/2. if Wing.FullWing else Wing.b
            
            y      = npy.linspace(0.0, b / IN, 30)*IN
            Load   = Wing.LoadDistSchreiner(y, V)
            Shear  = Wing.ShearLoading(y, V, Load)
            Moment = Wing.MomentLoading(y, V, ShearLoad = Shear)
            
            BStress   = self.BendingStress(y, V, Moments = Moment) / PSI
            
            MaxBStress = max(BStress)
            
            if self.MaxBendStress is not None:
                MaxAllowStress = self.MaxBendStress / PSI
            else:
                MaxAllowStress = None
            
            y = y / IN
            
            pyl.figure(fig)
            pyl.title("Bending Stress for " + Wing.name)
            pyl.plot(y, BStress)
            pyl.plot([y[0], y[-1]], [MaxAllowStress, MaxAllowStress])
            legend = [self.name + ' Max Stress %1.2f (psi)' % MaxBStress]
            if MaxAllowStress is not None:
                pyl.ylim(ymin = 0., ymax = 1.1*MaxAllowStress)
                legend.append(self.name + ' Max Allowed')
            
            pyl.xlabel("Semi Span (in)"); pyl.ylabel("Bending Stess (psi)")
            pyl.legend(legend, loc = 'best')
        
#===============================================================================
        def _CalcCentroid(self, y = 0*IN):
            """
            Assumes Centroid remains at relative X position all the way out span and returns this X at span location y
            """
            Centroid = npy.array([0, 0, 0]) * IN 
            
            Wing = self.param.Wing

            Centroid[0] = Wing.LE(y) + Wing.Chord(y)*self.Position[0]
            if abs(Wing.Axis[1]) == 1:
                Centroid[2] = Wing.X[2] + y #span direction
                Centroid[1] = Wing.Mid(y,Centroid[0]) + Wing.Thickness(y,Centroid[0])/2 * self.Position[1]
            else:
                Centroid[1] = Wing.X[1] + y #span direction
                Centroid[2] = Wing.Mid(y,Centroid[0]) + Wing.Thickness(y,Centroid[0])/2 * self.Position[1]
                
            return Centroid            


#===============================================================================
        def _CalcWeight(self):
            """
            Calculates the Weight of the spar using rough integration every inch
            """
            if self.dirty: self.Refresh()
            
            if not self.SparMat.NoneList.has_key('LinearForceDensity'):
                b = self.param.Wing.b
                return self.SparMat.LinearForceDensity*b* self.PctSpan
            
            A = self.Area
            
            sym = 1
            
            b = self.param.Wing.b
            if self.param.Wing.FullWing:
                b /= 2
                sym = 2
                
            sym = 2*sym if self.param.Wing.Symmetric else sym
            
            L = b * self.PctSpan
            Weight = 0*LBF
            
            for ay in range(int(math.floor(L / IN))):
                A = self.GetCrossSection(ay*IN)[2]
                dL = 1*IN
                Weight += dL*A*self.SparMat.ForceDensity * sym
            dL = ((L / IN) % 1)*IN
            A = self.GetCrossSection(L)[2]
            Weight += dL*A*self.SparMat.ForceDensity * sym
            
            return Weight
                        
#===============================================================================
        def _CalcI(self, y = 0*IN):
            """
            Calculates the area moment of inertia of the spar assuming square spar
            """
            if self.dirty: self.Refresh()
            
            y = self.MakeList(y)
            
            I = []
            
            for yi in y:
                Dim1,Dim2,Area = self.GetCrossSection(yi)
                
                Height = OD = Dim1 
                Width  = ID = Dim2
                                    
                if self.TubeSpar:
                    II = (math.pi / 64) * (OD**4 - ID**4)
                else:
                    II = Width * Height**3 / 12
                    
                I.append( II / IN**4 )
        
            return I[0] * IN**4 if len(I) == 1 else npy.array(I) * IN**4
            
#===============================================================================
        def _CalcDim1(self):
            """
            Calculates the height of the spar
            """
            if self.dirty: self.Refresh()
            
            
            Dim1 = self.Area / self.Dim2

            if self.TubeSpar is True:
                Dim1 = ((self.Area / math.pi) + self.Dim2**2)**0.5
                    
            return Dim1
            
#===============================================================================
        def _CalcDim2(self):
            """
            Calculates the width of the spar (Dim2)
            """
            if self.dirty: self.Refresh()
            
            Dim2 = self.Area / self.Dim1
    
            if self.TubeSpar is True:
                Dim2 = (self.Dim1**2 - (self.Area / math.pi))**0.5
    
            return Dim2
            
#===============================================================================
        def GetCrossSection(self, y = 0*IN):
            """
            Gets the crossection dimensions at a location along span
            """
            
            if self.dirty: self.Refresh()
            
            Centroid = self._CalcCentroid(y)
            d1 = self.Dim1
            d2 = self.Dim2
            
            # Check to see if Spar exceeds thickness of wing
            if self.param.Wing.Thickness(y,Centroid[0]) < self.Dim1 and self.ScaleToWing[0]:
                d1 = self.param.Wing.Thickness(y,Centroid[0])
                if self.ScaleToWing[1]:
                    d2 *= (self.param.Wing.Thickness(y,Centroid[0]) / self.Dim1)  #scale d2
            
            A = d2 * d1
                
            if self.TubeSpar is True:
                A = math.pi * (d1**2 - d2**2)
            
            return d1,d2,A
        
#===============================================================================
        def _CalcArea(self):
            """
            Calculates the crossectional area of the spar
            """
            if self.dirty: self.Refresh()
            
            d1 = self.Dim1
            d2 = self.Dim2
              
            Area = d2 * d1
                
            if self.TubeSpar is True:
                Area = math.pi * (d1**2 - d2**2)
            
            return Area

#===============================================================================
        def __getattr__(self, key):
            """
            Attempts to calculate any variables set to None
            """
            ans = super(ACWingWeight.Spar, self).__getattr__(key)
            if ans is not None:
                return ans
#    
            # Try to calculate the value instead
            if key == 'Area':
                return self._CalcArea()
            elif key == 'Dim1':
                return self._CalcDim1()
            elif key == 'Dim2':
                return self._CalcDim2()
            elif key == 'I':
                return self._CalcI()
    
            #
            # If no calculation exist return None
            #
            return None
        
#===============================================================================
        def _CheckConsistent(self):
            """
            Makes sure that a consistent set of parameters has been specified.
            """
            self._CheckEquation(['Area','Dim1','Dim2'], Need = 2)
        
#===============================================================================
        def Refresh(self):
            """
            Virtual function used to recompute attributes

            All derived instances of this methods must set
            self.param.refreshing = True at the beginning of the method and
            self.param.refreshing = False when completed.

            Otherwise an infinite loop will be created!!!
            """
            super(ACWingWeight.Spar, self).Refresh()
            self.param.refreshing = True
            
            Wing = self.param.Wing

            self._CheckConsistent()
            
            self.CopyVisibility(Wing)
            
            self.X = self._CalcCentroid(0*IN)
            
            self.UpdateDrawing()

            self.param.refreshing = False
            
#===============================================================================
        def UpdateDrawing(self):
            """
            Updates all drawings of the spar
            """

            SparDrawing = self.Drawings.SparDrawing
            Wing = self.param.Wing
        
            xpos = self.X[0] / IN
            if abs(Wing.Axis[1]) == 1:  #Vertical Wing
                #Swap z and y since vertical wings have local z and y swapped from global system...
                ypos = self.X[2] / IN
                zpos = self.X[1] / IN
            else:
                ypos = self.X[1] / IN
                zpos = self.X[2] / IN
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            def GenCircle(xy, r, n):
                coordx = npy.empty(n+1)
                coordy = npy.empty(n+1)
                r      = r / IN
                
                for npt in range(n+1):
                    theta = npt/n * 2*math.pi
                    coordx[npt] = math.cos(theta) / r + xy[0]
                    coordy[npt] = math.sin(theta) / r + xy[1]
                
                return (coordx,coordy)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            
            
            #
            # DRAW Crossection View of Spar (plots local x vs local z)
            #
            if self.TubeSpar is True:
                npts = 31 # Number of points to draw the ellipse
                xID,yID = GenCircle((xpos,zpos), self.Dim2, npts)
                xOD,yOD = GenCircle((xpos,zpos), self.Dim1, npts)
                x = npy.append(xOD,xID)
                y = npy.append(yOD,yID)
                
        
            else: #rectangle
                x = npy.empty(5)
                y = npy.empty(5)
                
                H = self.Dim1 / IN
                W = self.Dim2 / IN
                
                x[0] = x[4] = xpos - W/2    #upper left
                y[0] = y[4] = zpos + H/2
                x[1] = xpos - W/2           #lower left
                y[1] = zpos - H/2
                x[2] = xpos + W/2           #lower right
                y[2] = zpos - H/2
                x[3] = xpos + W/2           #upper right
                y[3] = zpos + H/2
                
            XSectionX = x
            XSectionY = y
             
            x = y = None #clear the variables
               
            #
            # DRAW Planform View of Spar (plots local x vs local y)
            #
            x = npy.empty(5)
            y = npy.empty(5)
            
            Wr = self.Dim2 / IN #rectangle spar width is Dim2
            if self.TubeSpar is True:
                Wr = self.Dim1 / IN #tubespar width is Dim1
            
            b = Wing.b
            b = b/2. if Wing.FullWing else b

            spartip = (b*self.PctSpan)
            
            Wt = self.GetCrossSection(spartip)[1] / IN
            if self.TubeSpar is True:
                Wt = self.GetCrossSection(spartip)[0] / IN #tubespar width is Dim1
                
            xtip = (Wing.LE(spartip) + (Wing.Chord(spartip)*self.Position[0])) / IN
            ytip = ypos + spartip / IN
            
            x[0] = x[4] = xpos - Wr/2   #front root
            y[0] = y[4] = ypos
            x[1] = xpos + Wr/2          #back root
            y[1] = ypos
            x[2] = xtip + Wt/2          #back tip
            y[2] = ytip
            x[3] = xtip - Wt/2          #front tip
            y[3] = ytip
            
            PlanformX = x
            PlanformY = y
            
            if abs(Wing.Axis[1]) == 1:  #Vertical Wing
                #Swap z and y since vertical wings have local z and y swapped from global system...
                SparDrawing.xt = XSectionX
                SparDrawing.yt = XSectionY
                SparDrawing.xs = PlanformX
                SparDrawing.zs = PlanformY
            else:
                SparDrawing.xs = XSectionX
                SparDrawing.zs = XSectionY
                SparDrawing.xt = PlanformX
                SparDrawing.yt = PlanformY
                
            pass
        
#===============================================================================
#
#   ACWingWeight Begins here!
#
#===============================================================================
    def __init__(self, Wing):
        super(ACWingWeight, self).__init__()
        UnitList = self.UnitList

        #
        # Save of the wing
        #
        self.param.Wing = Wing

        self.__dict__['SkinMat'] = ACMaterial()
        self.__dict__['DrawDetail'] = False
        self.__dict__['DSparRootSpar'] = None
        
        #
        # Set the material to have a zero weight by default
        #
        self.SkinMat.AreaForceDensity   = 0 * LBF/IN**2
        
        # The Spars in the wing
        self.__dict__['Spars'] = []

#===============================================================================
    def AddSpar(self, SparName, Height, Width, Position = (0.25,0), PctSpan = 1.0, DSpar = False, Mirror = False, Structural = True):
        """
        Adds a new SquareSpar to the wing
        Inputs: 
            SparName
            Height
            Width
            Position    (Defaults to quarterchord, camber line)
            PctSpan     (Defaults to full span)
            DSpar       Add a D-spar to wing that wraps around LE to this spar
            Mirror      Mirrors the spar dimensions across camber line (spar caps)  
                    NOTE: This creates two spars at once, appending _1 and _2 to the SparName. 
                         _2 is the mirror and modifications made to one spar WILL NOT TRANSFER
                         They become two separate spars.
        """
        
        Spars = self.Spars
        Wing = self.param.Wing
        
        #
        # Create the spar in the wing class so it can be referenced later
        #
        Spars.append(ACWingWeight.Spar(Wing))
        Spars[-1].Dim1     = Height
        Spars[-1].Dim2     = Width
        Spars[-1].Position = Position
        Spars[-1].PctSpan  = PctSpan
        Spars[-1].name     = SparName
        Spars[-1].DrawDetail = self.DrawDetail
        Spars[-1].Structural = Structural
        
        if Mirror is True:
            #modify name of previous spar
            Spars[-1].name = SparName + "_1"
            # Save of a reference to the first spar to make so the user can type
            self.__dict__[Spars[-1].name] = Spars[-1]
            
            Spars.append(ACWingWeight.Spar(Wing))
            Spars[-1].Dim1     = Height
            Spars[-1].Dim2     = Width
            Spars[-1].Position = (Position[0],-Position[1])
            Spars[-1].PctSpan  = PctSpan
            Spars[-1].name     = SparName + "_2"
            Spars[-1].DrawDetail = self.DrawDetail
            Spars[-1].Structural = Structural
        
        # Save of a reference to the spar to make so the user can type
        self.__dict__[Spars[-1].name] = Spars[-1]
        
        if DSpar is True:
            self.DSparRootSpar = self.Spars[-1].name
        

#===============================================================================
    def AddTubeSpar(self, SparName, OD, ID, Position = (0.25,0), PctSpan = 1.0):
        """
        Adds a new TubeSpar to the wing
        Inputs: 
            SparName
            Outer Diameter
            Inner Diameter
            Position    (Defaults to quarterchord, camber line)
            PctSpan     (Defaults to full span)
        """
        
        self.AddSpar(SparName,OD,ID,Position = Position,PctSpan = PctSpan)
        
        Spars = self.Spars
        Spars[-1].TubeSpar = True

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the wing weight parts to the Weight Table
        
        Input:
            PartName    - The name given to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        WingWeight = self
        class CSkin:
            def CG(self):      
                return WingWeight.CG()
            def MOI(self):
                return npy.array([0,0,0])*LBM*IN**2          
                
        Skin = CSkin()
        Skin.name   = 'Skin'
        Skin.Weight = self.SkinWeight()
        Skin.WeightGroup = WingWeight.WeightGroup
                 
        WeightTable[PartName] = Skin
                
        for Spar in self.Spars:
            Spar.AddToWeightTable(Spar.name, WeightTable)
        

#===============================================================================
    def PlotBendingStress(self, V, fig = 1):
        """
        Plots the bending stress of all the spars

        Inputs:
            V   - Velocity to compute loads
            fig - Figure Number
        """
        if self.dirty: self.Refresh()
        
        for Spar in self.Spars:
            Spar.PlotBendingStress(V, fig)

#===============================================================================
    def SkinWeight(self):
        """
        Calculates the Skin weight
        """
        Wing = self.param.Wing
        
        return self.SkinMat.AreaForceDensity * Wing.SWet()
    
#===============================================================================
    def SparWeight(self):
        """
        Calculates the weight all the spars
        """
        Spars = self.Spars
        
        Weight = 0*LBF
        
        #
        # Calculate the weight of the spar
        #
        for Spar in Spars:
            Weight += Spar.Weight
            
        return Weight
    
#===============================================================================
    def _CalcWeight(self):
        """
        Calculates the weight of the Spar and Skin
        """
        #
        # Calculate the weight of the spar
        #
        Weight = self.SparWeight()

        #
        # Add the weight of the skin
        #
        Weight += self.SkinWeight()
        
        if self.param.Wing.Symmetric:
            Weight *= 2

        return Weight 

#===============================================================================
    def _CalcI(self, y = 0 * IN):
        """
        Calculates the sum area moment of inertia for the wing about the camber line at 25% chord
        """
        Spars = self.Spars
        Wing = self.Wing
        
        WingI = 0 * IN**4
        
        for Spar in Spars:
            if Spar.Structural is False:
                continue
            
            SC = Spar._CalcCentroid(y)
            SH,SW,SA = Spar.GetCrossSection(y)
            SI = Spar._CalcI(y)
            
            if abs(Wing.Axis[1]) == 1:  #Vertical Wing
                WingI += SI + SA * (SC[1] - Wing.Mid(y,0.25*Wing.Chord(y)))**2
            else:
                WingI += SI + SA * (SC[2] - Wing.Mid(y,0.25*Wing.Chord(y)))**2
            
        return WingI
    
#===============================================================================
    def MOI(self):
        """
        Sums the MOI of the spars of the wing
        """
        
        Spars = self.Spars
        
        CG = self.CG()
        
        MOI = npy.array([0,0,0])*LBM*IN**2
        
        def AddMOI(PCG, PM, PMOI):
            dCG = PCG - CG
            MOI[0] += PMOI[0] + PM*(dCG[1]**2 + dCG[2]**2) #Ixx
            MOI[1] += PMOI[1] + PM*(dCG[0]**2 + dCG[2]**2) #Iyy
            MOI[2] += PMOI[2] + PM*(dCG[1]**2 + dCG[0]**2) #Izz
        
        for Spar in Spars:
            AddMOI(Spar.CG(), Spar.Weight / g, Spar.MOI())
            
        return MOI
    
#===============================================================================
#    def SparBendingStress(self, y, V, Moments = None):
#        """
#        Calculates the bending stress of the spar
#        
#        Input:
#            y - Spanwise location(s)
#            V - Velocity
#        """
#        if self.dirty: self.Refresh()
#        
#        Wing = self.param.Wing
#        
#        if Moments is not None:
#            M = Moments
#        else:
#            M = self.MomentLoading(y, V)
#        
#        h = []
#        for ay in y:
#            h.append(Wing.Thickness(ay, Wing.LE(ay) + self.SparFc*Wing.Chord(ay)) / (IN))
#        
#        h = h[0]*IN if len(h) == 1 else npy.array(h)*IN
#        I = self.SparI(y)
#        
#        return M*h/(2.*I)

#===============================================================================
#    def SparMaxBendingStress(self, V):
#        """
#        Calculates maximum the bending stress of the spar
#        
#        Input:
#            V - Velocity
#        """        
#        return self.SparBendingStress(0*IN, V)

#===============================================================================
#    def PlotBendingStress(self, V, fig = 1):
#        """
#        Plots the bending stress of the spar
#
#        Inputs:
#            V   - Velocity to compute loads
#            fig - Figure Number
#        """
#        if self.dirty: self.Refresh()
#        
#        Wing = self.param.Wing
#        b = Wing.b/2. if Wing.FullWing else Wing.b
#        
#        y      = npy.linspace(0.0, b / (IN), 30)*IN
#        Load   = Wing.LoadDistSchreiner(y, V)
#        Shear  = Wing.ShearLoading(y, V, Load)
#        Moment = Wing.MomentLoading(y, V, ShearLoad = Shear)
#        
#        BStress   = self.SparBendingStress(y, V, Moments = Moment) / (PSI)
#        
#        MaxBStress = max(BStress)
#        
#        if self.SparMaxBendStress is not None:
#            MaxAllowStress = self.SparMaxBendStress / (PSI)
#        else:
#            MaxAllowStress = None
#        
#        y = y / (IN)
#        
#        pyl.figure(fig)
#        pyl.title("Bending Stress for " + Wing.name)
#        pyl.plot(y, BStress)
#        pyl.plot([y[0], y[-1]], [MaxAllowStress, MaxAllowStress])
#        legend = ['Bending Stress %1.2f (psi)' % MaxBStress]
#        if MaxAllowStress is not None:
#            pyl.ylim(ymin = 0., ymax = 1.1*MaxAllowStress)
#            legend.append('Maximum Allowed Stress')
#        
#        pyl.xlabel("Semi Span (in)"); pyl.ylabel("Bending Stess (psi)")
#        pyl.legend(legend, loc = 'best')
        
        
#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACWingWeight, self).Refresh()
        self.param.refreshing = True
        #
        # Perform all calculations
        #
        
        for Spar in self.Spars:
            Spar.Refresh()
            
        self.param.refreshing = False
        
#===============================================================================
    def Draw(self, fig = 1, top = 221, side = 222, front = 223):
        """
        Draws the wing weight

        Inputs:
            fig - Integer number of figure to draw into
            top - Subplot for top view
            side - Subplot for side view
            front - Subplot for front view
        """
        super(ACWingWeight, self).Draw(fig, top, side, front)
        
        if self.DrawDetail is True:
            for Spar in self.Spars:
                Spar.Draw(fig, top, side, front)
            
################################################################################
class ACSolidWing(ACWingWeight):
    """
    Computes the weight of a solid foam wing

    Attributes:
        WingMat  - ACMaterial class defining the material of the solid wing
    """
#===============================================================================
    def __init__(self, Wing):
        super(ACSolidWing, self).__init__(Wing)
        UnitList = self.UnitList

        #
        # Save of the wing
        #
        self.param.Wing = Wing
        self.__dict__['WingMat'] = ACMaterial()
        
        #
        # Set the material to have a zero density
        #
        self.WingMat.ForceDensity = 0*LBF/IN**3
        
        self.param.Weight = 0*LBF
        
#===============================================================================
#    def _IntegrateMOI(self, iter = 5, tol = 10.0):
#        
#        #DO NOT USE UNTIL FIXED: Currently too computationally complex and does not correctly apply both wings.
#        
#        
#        Wing = self.param.Wing
#        rho = self.WingMat.Density 
#        
#        MOI = npy.array([0,0,0])*LBM*IN**2
#
#        #Redefine Equations to deal with UNUM issues...
#        def Upper(y,x):
#            ans = wing.Upper(y*IN,x*IN)
#            return ans / (IN)
#        def Lower(y,x):
#            ans = wing.Lower(y*IN,x*IN)
#            return ans / (IN)
#        def LE(y):
#            ans = wing.LE(y*IN)
#            return ans / (IN)
#        def TE(y):
#            ans = wing.TE(y*IN)
#            return ans / (IN)
#        
#        b = Wing.b / (IN)
#            
#        def Ixx(x,y,z):
#            return y**2 + z**2
#        def Iyy(x,y,z):
#            return x**2 + z**2
#        def Izz(x,y,z):
#            return y**2 + x**2
#
#        def IntChord(x,y,I):
#            def Thick(z):
#                return I(x,y,z)  
#            int,err = integrate.quadrature(Thick, Lower(y,x), Upper(y,x), tol=tol, maxiter=iter, vec_func = False)  #Integrates function I across thickness
#            return int
#
#        def IntSpan(y, I):
#            def Chord(x):
#                return IntChord(x,y,I)
#            int,err = integrate.quadrature(Chord,LE(y), TE(y), tol=tol, maxiter=iter, vec_func = False) #Integrates integral of function I across thickness across chord
#            return int
#
#        def IntWing(I):
#            def Span(y):
#                return IntSpan(y,I)
#            int,err = integrate.quadrature(Span,0, b/2, tol=tol, maxiter=iter, vec_func = False) #Integrates integral of function I across thickness and chord across span
#            return int*IN**5
#
#        # Units are IN**5 because triple integral of IN**2 values
##        def IntWing(I):
##            return integrate.tplquad(I, 0, b/2, LE, TE, Lower, Upper, epsabs = 10000., epsrel = 10000.)*IN**5
#        
#        MOI[0] = IntWing(Ixx)*rho
#        MOI[1] = IntWing(Iyy)*rho
#        MOI[2] = IntWing(Izz)*rho
#        
#        return MOI

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
        super(ACSolidWing, self).Draw(fig, top, side, front)
        #
        # This just draws the equivalent mass box representation of the wing.
        # Mostly for debugging purposes
        #
        
        Panels = self._BuildMassBoxEquivalent() #get Massbox simplified wing
        for panel in Panels:
            panel.Draw(fig, top, side, front)

#===============================================================================
    def _BuildMassBoxEquivalent(self):
        """
        Creates an equivalent wing out of massboxes
        """
        if self.dirty:
            self.param.Wing.Refresh()
        
        MBWingPanels = []
        Wing    = self.param.Wing
        rho     = self.WingMat.ForceDensity
        AFArea  = Wing.af.Area
        WingX   = Wing.GetX()
        WWCF    = 1.05  # Wing Weight Correction Factor
        
        npanel = len(Wing.param.span) - 1
        
        crd = Wing.param.chord
        spn = Wing.param.span
        dih = Wing.param.dihed
        swp = Wing.param.sweep 
        
        for panel in range(npanel):
            
            pspn = (spn[panel+1] - spn[panel]) * IN         #span
            pcrd = (crd[panel+1] + crd[panel]) / 2 * IN     #average chord
            px   = (swp[panel+1] + swp[panel]) / 2 * IN     #x offset from WingX
            py   = (spn[panel+1] + spn[panel]) / 2 * IN     #y offset from WingX
            pz   = (dih[panel+1] + dih[panel]) / 2 * IN     #z offset from WingX 
            
            #Use MassBox class for rough calculations
            #
            # MassBox is constructed such that its CG is roughly at same location on wing.
            #
            for i in range(2):
                if Wing.FullWing is True and i==1:
                    py = -py
                elif Wing.FullWing is False and i==1:
                    break
                
                MBWingPanels.append(ACMassBox())
                if Wing.Axis == [0,1]:
                    MBWingPanels[-1].Axis = (0, 0, 1)
                    MBWingPanels[-1].X    = [px + WingX[0], pz + WingX[1], py + WingX[2]]     # Place X relative to wing Local CSys.
                    MBWingPanels[-1].LWH  = [pspn, pcrd * Wing.af.Thickness()*WWCF, pcrd/2]  # Uses average chord with the length corrected to put
                    MBWingPanels[-1].Xat  = [0.5, 1.0, 0.5]
                else:
                    MBWingPanels[-1].Axis = (0, 1, 0)
                    MBWingPanels[-1].X    = [px + WingX[0], py + WingX[1], pz + WingX[2]]     # Place X relative to wing Local CSys.
                    MBWingPanels[-1].LWH  = [pspn, pcrd/2, pcrd * Wing.af.Thickness()*WWCF]  # Uses average chord with the length corrected to put
                    MBWingPanels[-1].Xat  = [0.5, 0.5, 0.0] 
                MBWingPanels[-1].Weight = rho * Wing.af.Area * pcrd**2 * pspn
                MBWingPanels[-1].Symmetric = Wing.Symmetric
            
        return MBWingPanels
        
#===============================================================================
    def CG(self):
        """
        Calculates the center of gravity of the foam wing relative to wing
        """
        
        Panels = self._BuildMassBoxEquivalent() #get Massbox simplified wing
        
        Wing  = self.param.Wing
        Weight = 0*LBF
        LocCG = npy.array([0,0,0])*IN*LBF
        CG    = npy.array([0,0,0])*IN
        X     = Wing.GetX()
        
        for Panel in Panels:
            Weight += Panel.Weight
            LocCG += Panel.CG()*Panel.Weight
        
        CG = LocCG / Weight
        
        return CG

#===============================================================================
    def MOI(self):
        """
        Calculates an approximate moments of inertia for the wing using massboxes
        """
        
        MOI = super(ACSolidWing, self).MOI()  #returns MOI due to spar in local coordinates
        
        Panels = self._BuildMassBoxEquivalent()  #get Massbox simplified wing
        
        Wing   = self.param.Wing
        LocCG  = self.CG()
        LocMOI = npy.array([0,0,0])*LBM*IN**2
        
        for Panel in Panels:
            sym  = 2 if Panel.Symmetric else 1
            pCG  = Panel.CG()
            pMOI = Panel.MOI()
            pM   = Panel.Weight * sym / g
            dCG  = pCG - LocCG
            
            LocMOI[0] += pMOI[0] + pM*(dCG[1]**2 + dCG[2]**2)#Ix
            LocMOI[1] += pMOI[1] + pM*(dCG[0]**2 + dCG[2]**2)#Iy
            LocMOI[2] += pMOI[2] + pM*(dCG[1]**2 + dCG[0]**2)#Iz
        
        
        MOI += LocMOI
        
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
        super(ACSolidWing, self).Refresh()
        self.param.refreshing = True
        #
        # Perform all calculations
        #
        
        self.CopyVisibility(self.param.Wing)
 
        self.param.refreshing = False

#===============================================================================
    def PannelWeight(self):
        """
        Calculates the weight of the mass box pannels
        """
        Weight = super(ACSolidWing, self).Weight()  #get spar and skin weight
        Wing = self.param.Wing
        
        Panels = self._BuildMassBoxEquivalent()  #get Massbox simplified wing

        sym  = 2 if Wing.Symmetric else 1
        
        for Panel in Panels:
            Weight += sym * Panel.Weight

        return Weight 

#===============================================================================
    def _CalcWeight(self):
        """
        Adds Weight of structure to spar and skin weight
        """
        Weight = super(ACSolidWing, self)._CalcWeight()  #get spar and skin weight
        Wing = self.param.Wing
        Vol = Wing.Volume()

        sym  = 2 if Wing.Symmetric else 1
        
        #
        # Add the weight of a solid wing to the skin and spar
        #
        Weight += Vol * self.WingMat.ForceDensity * sym

        return Weight 

################################################################################
class ACRibWing(ACWingWeight):
    """
    Computes the weight of a wing made from ribs

    Attributes:
        RibMat   - ACMaterial class defining the material of the ribs
        RibSpace - The maximum rib spacing
    """

#===============================================================================
    class Rib(ACComponent):

        def __init__(self, Wing):
            super(ACRibWing.Rib, self).__init__()
            UnitList = self.UnitList

            self.__dict__['Weight'] = 0*LBF               ; UnitList['Weight'] = Force
            self.__dict__['C']      = 0*IN                ; UnitList['C']      = Length
            self.__dict__['CG']     = [0*IN, 0*IN, 0*IN]  ; UnitList['CG']     = Length
            self.__dict__['LE']     = 0*IN                ; UnitList['LE']     = Length
            self.__dict__['Axis']   = (1, 0)              ;   UnitList['Axis'] = Unitless

#            self.ShowSide  = False
#            self.ShowFront = False
            
            self.param.Wing = Wing
            self.Drawings.Rib = ACComponent.ComponentDrawing()

#===============================================================================
        def Refresh(self):
            """
            Virtual function used to recompute attributes

            All derived instances of this methods must set
            self.param.refreshing = True at the beginning of the method and
            self.param.refreshing = False when completed.

            Otherwise an infinite loop will be created!!!
            """
            super(ACRibWing.Rib, self).Refresh()
            self.param.refreshing = True
            #
            # The rib is only drawn for the top view
            #
            Rib = self.Drawings.Rib

            LE = self.LE / (IN)
            C  = self.C / (IN)

            if abs(self.Axis[0]) == 1: #Horizontal surface
                y  = self.CG[1] / (IN)
                Rib.xt = npy.array([LE,LE+C])
                Rib.yt = npy.array([y,y])
            else: #Vertical surface
                z  = self.CG[2] / (IN)
                Rib.xs = npy.array([LE,LE+C])
                Rib.zs = npy.array([z,z])
                
            self.CopyVisibility(self.param.Wing)

            self.param.refreshing = False

#===============================================================================
        def _DrawTop(self, Drawing):
            """
            Draw the y location of the ribs
            """
            super(ACRibWing.Rib, self)._DrawTop(Drawing)

            if self.Axis[0] == 0: #Not Horizontal surface
                return
            #
            # Get the current axis
            #
            ax = pyl.gca()

            y = self.CG[1]
            TE = (self.LE + self.C*1.5) / (IN)

            lbl = AsUnit( y, "in", '%3.1f ')
            
            ax.annotate(lbl, xy=(TE, y / IN), xycoords='data')

#===============================================================================
        def _DrawSide(self, Drawing):
            """
            Draw the y location of the ribs
            """
            super(ACRibWing.Rib, self)._DrawSide(Drawing)

            if self.Axis[1] == 0: #Not Vertical surface
                return
            #
            # Get the current axis
            #
            ax = pyl.gca()

            z = self.CG[2]
            TE = (self.LE + self.C*1.25) / IN

            lbl = AsUnit( z, "in", '%3.1f ')
            
            ax.annotate(lbl, xy=(TE, z / IN), xycoords='data')

#===============================================================================

#===============================================================================
#   ACRibWing continues here
#===============================================================================
    def __init__(self, Wing):
        super(ACRibWing, self).__init__(Wing)
        UnitList = self.UnitList


        self.__dict__['RibMat'] = ACMaterial()
        self.__dict__['RibSpace'] = 5*IN  ;  UnitList['RibSpace'] = Length
        self.__dict__['DrawRibs'] = False

        # Give it a better name than just the name of the class
        self.name = 'Ribs'
        
        #
        # Storage list for the ribs
        #
        self.param.Ribs = []

        self.ShowSide  = False
        self.ShowFront = False

#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        super(ACRibWing, self).Refresh()
        self.param.refreshing = True
        #
        # Perform all calculations
        #
        self._CalcRibs()
        
        self.CopyVisibility(self.param.Wing)

        self.param.refreshing = False

#===============================================================================
    def _CalcRibs(self):
        """
        Calculates rib positions based on the maximum rib spacing
        """
        Wing   = self.param.Wing
        
        if Wing.dirty:
            Wing.Refresh()

        span   = Wing.param.span * IN
        AFArea = Wing.af.Area
        self.param.Ribs = []
        Ribs   = self.param.Ribs


        for i in range(len(span)-1):
            #
            # Compute the number of ribs that would fit this spacing
            #
            nRib = int(math.ceil(((span[i+1]-span[i])/self.RibSpace)))
            #
            # Once the number of ribs is determined, find the true spacing
            #
            dRib = (span[i+1]-span[i])/nRib
            for j in range(nRib):
                Ribs.append(self.Rib(Wing))
                y = float(j)*dRib + span[i]
                C = Wing.Chord(y)
                Area = AFArea*C**2

                Ribs[-1].C      = C
                Ribs[-1].Axis   = Wing.Axis
                Ribs[-1].Weight = Area * self.RibMat.AreaForceDensity
                Ribs[-1].LE     = Wing.LE(y)
                if abs(Wing.Axis[0]) == 1:
                    Ribs[-1].CG     = [Wing.LE(y) + 0.25*C, y + Wing.X[1], Wing.Mid(y)]
                else:
                    Ribs[-1].CG     = [Wing.LE(y) + 0.25*C, Wing.Mid(y), y + Wing.X[2]]

        #
        # Finish of the wing tip rib
        #
        Ribs.append(self.Rib(Wing))
        y = span[-1]
        C = Wing.Chord(y)
        Area = AFArea*C**2

        Ribs[-1].C      = C
        Ribs[-1].Axis   = Wing.Axis
        Ribs[-1].Weight = Area * self.RibMat.AreaForceDensity
        Ribs[-1].LE     = Wing.LE(y)
        if abs(Wing.Axis[0]) == 1:
            Ribs[-1].CG     = [Wing.LE(y) + 0.25*C, y + Wing.X[1], Wing.Mid(y)]
        else:
            Ribs[-1].CG     = [Wing.LE(y) + 0.25*C, Wing.Mid(y), y + Wing.X[2]]


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
        super(ACRibWing, self).Draw(fig, top, side, front)
        Ribs = self.param.Ribs
        
        if self.DrawRibs is False:
                return

        for i in range(len(Ribs)):
            Ribs[i].Draw(fig, top, side, front)

#===============================================================================
    def CG(self):
        """
        Calculates the center of gravity of the rib wing
        """
        
        Wing = self.param.Wing
        TempCGW = npy.array([0,0,0])*IN*LBF
        CG = npy.array([0,0,0])*IN
        TempWeight = 0.0 * LBF
        
        if self.dirty:
            self.Refresh()
            
        Ribs = self.param.Ribs
        for i in range(len(Ribs)):
            Ribs_i_CG   = Ribs[i].ToUnumList(Ribs[i].CG, IN)
            TempCGW    += Ribs_i_CG * Ribs[i].Weight
            TempWeight += Ribs[i].Weight

        CG = TempCGW / TempWeight

        if Wing.FullWing:
            ax = 1 if abs(Wing.Axis[0]) == 1 else 2
            CG[ax] = Wing.X[ax]
        
        if Wing.Symmetric:
            CG[1] = 0*IN

        return CG

#===============================================================================
    def RibMOI(self):
        """
        Calculates the moments of inertia of the rib wing
        """
        MOI = npy.array([0,0,0])*LBM*IN**2
        
        Wing = self.param.Wing
        CG   = self.CG()
        
        if self.dirty:
            self.Refresh()
                   
        Ribs = self.param.Ribs
        for i in range(len(Ribs)):
            Rx = Ribs[i].CG[0] - CG[0]
            Ry = Ribs[i].CG[1] - CG[1]
            Rz = Ribs[i].CG[2] - CG[2]
                
            sym  = 2 if Wing.Symmetric else 1
            full = 2 if Wing.FullWing else 1
                
            MOI[0] += full * sym * Ribs[i].Weight * (Rz**2 + Ry**2) / g
            MOI[1] += full * sym * Ribs[i].Weight * (Rz**2 + Rx**2) / g
            MOI[2] += full * sym * Ribs[i].Weight * (Rx**2 + Ry**2) / g
        
        return MOI


#===============================================================================
    def MOI(self):
        """
        Calculates the moments of inertia of the rib wing
        """
        MOI = super(ACRibWing, self).MOI()
        
        MOI += self.RibMOI()
        
        return MOI

#===============================================================================
    def RibWeight(self):
        """
        Calculates the weight of the ribs
        """
        Wing   = self.param.Wing
        
        Ribs = self.param.Ribs
        RibWeight = 0*LBF
        for i in range(len(Ribs)):
            RibWeight += Ribs[i].Weight

        sym  = 2 if Wing.Symmetric else 1
        full = 2 if Wing.FullWing else 1

        RibWeight *= full * sym
            
        return RibWeight        
        
#===============================================================================
    def _CalcWeight(self):
        """
        Adds Weight of structure to spar and skin weight
        """
        
        Weight = super(ACRibWing, self)._CalcWeight() #get spar and skin weight    
        
        return Weight + self.RibWeight()

#===============================================================================
    def AddToWeightTable(self, PartName, WeightTable):
        """
        Adds the wing weight parts to the Weight Table
        
        Input:
            PartName    - The name given to this part
            WeightTable - The table to add self to 
        """
        if self.dirty: self.Refresh()
        
        super(ACRibWing, self).AddToWeightTable(PartName, WeightTable)
        
        RibWing = self

        class CRibs:
            def CG(self):      
                return RibWing.CG()
            def MOI(self):
                return RibWing.RibMOI()          
                
        
        Ribs = CRibs()
        Ribs.name   = self.name
        Ribs.Weight = self.RibWeight()
        Ribs.WeightGroup = self.WeightGroup
                         
        WeightTable[PartName] = Ribs
        

################################################################################
if __name__ == '__main__':

    from ACWing import ACMainWing
    from DefaultMaterialsLibrary import PinkFoam

    RibBalsa = ACMaterial()
    RibBalsa.ForceDensity = 0.00577 * LBF/FT**3
    RibBalsa.Thickness = 1/8*IN
    RibBalsa.E = 500000 * PSI

    wing = ACMainWing(1)
    wing.Lift_LO       = 26.799 * LBF
    wing.V_max_climb   = 65 * FT/SEC
    wing.V_Stall       = 39.2 * FT/SEC
    wing.Alt_LO        = 600 * FT
    
    wing.X = [30*IN, 10*IN, 10*IN]
    wing.Axis = [0,1]
#    wing.Symmetric = True

#    wing.S = 900 * IN**2
#    wing.AR = 8
    wing.b = 10. * FT
    wing.TR = [1, 0.5, 0.4]
    wing.Fb = [0.2, 0.8, 1]
    wing.Gam = [0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG]
#    wing.Lam = [3.0 * ARCDEG, 10.0 * ARCDEG, 4.0 * ARCDEG]
    wing.CEdge = 'LE'
    wing.ConstUpper = True
    wing.Airfoil = 'S1223'
    wing.FullWing = True
    wing.Symmetric = True
#    wing.color = 'r'
    wing.LLFile = 'LLTPolar.txt'


#    RibWing = wing.SetWeightCalc(ACRibWing)
#    RibWing.RibMat = RibBalsa.copy()
#    RibWing.RibSpace = 5 * IN
#    MOI = RibWing.MOI()
#    wing.ShowSide = False
#    wing.ShowFront = False
#    wing.Draw(fig=1,top=111)
#    RibWing.Draw(fig=1,top=111)

    SolidWing = wing.SetWeightCalc(ACSolidWing)
    SolidWing.WingMat = PinkFoam.copy()
    MOI = SolidWing.MOI() * g
    CG = SolidWing.CG()
#    wing.ShowSide = False
#    wing.ShowFront = False
    wing.Draw(fig=1) #,top=111)
    SolidWing.Draw(fig=1)

    #SolidWing.SparMaxBendStress = 30e3 *LBF/IN**2
    #SolidWing.PlotBendingStress(40*FT/SEC, fig = 5)

    print "MOI :", MOI
    print "CG  :", CG
    print "Weight :", SolidWing.Weight()
#    print "Weight :", RibWing.Weight()

    pyl.show()



