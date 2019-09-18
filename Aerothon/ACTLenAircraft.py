from __future__ import division # let 5/2 = 2.5 rather than 2
from ACAircraft import ACTailAircraft
from ACBase import Length, Unitless, ARCDEG
from scalar.units import IN
from ACTail import ACTailSurf

#
# Create an aircraft class to impose the total length restriction
#
class ACTLenAircraft(ACTailAircraft):
    """
    An aircraft class restricted by the sum of its length, width and height.
    
    Attributes:
        TotalLengths - The sum of the length width and height of the aircraft
        
    """
#===============================================================================
    def __init__(self):
        
        super(ACTLenAircraft, self).__init__()
        self.__dict__['TotalLengths']  = 173*IN   ; self.UnitList['TotalLengths'] = Length
        self.__dict__['VTail2Pos']         = 0    ; self.UnitList['VTail2Pos']    = Unitless
        self.__dict__['VTail2Root']        = None ; self.UnitList['VTail2Pos']    = Unitless

        self.__dict__['VTail2'] = ACTailSurf(self.Wing, 4)
        self.VTail2.Parents.append( self )
        self.VTail2.name = "Vertical Tail 2"
        
        rud = self.VTail2.AddControl('Rudder')

        #
        # Rudder typically extends the entire vertical tail
        #
        rud.Fb = 1.
        rud.Fc = 0.25
        rud.Ft = 0.
        
#===============================================================================
    def SetWing(self, Wing):
        """
        Assigns a new wing to the aircraft
        
        Inputs:
            Wing - The new wing to assign to the aircraft
        """
        super(ACTLenAircraft, self).SetWing(Wing)
        #
        # Update the references in the second vertical tail surface
        #
        if self.VTail2 is not None:
            self.VTail2.SetWing( Wing )
        

#===============================================================================
    def _PositionParts(self):
        """
        Impose the total length restriction
        """        
        self._PositionLandingGear()
        self._PositionWing()
        self._PositionTail()
        
        Xcg = self.Xcg()
        
        H = 1*IN
        Hold = 2*IN
        Crht = 1*IN
        Crhtold = 2*IN
        #
        # Because the height of the vertical tail may change, iterate on the height
        #
        while(abs(H - Hold) + abs(Crht - Crhtold) > 0.01*IN):
            W  = max(self.Wing.b,self.HTail.b)
            Hold    = H
            Crhtold = Crht
            
            H  = max(self.Wing.Upper(0*IN), self.VTail.Tip()[2])
 
            #TODO: Need to fix this so it is correct in the future
            #Crht = max( self.HTail.Chord(0*IN), self.VTail.Chord(0*IN) )
            Crht = self.HTail.Chord(0*IN)
            
            #Need to set VTail.L here to get an estimate for the VTail chord
            self.HTail.L = self.VTail.L = max( self.TotalLengths - Xcg - W - H - 3/4*Crht, self.Wing.TE(0*IN) - Xcg )
                        
            #Now move the VTail so it does not exted beyond the HTail
            self.VTail.L = self.HTail.TE(0*IN) - self.VTail.Chord(0*IN)*0.75 - Xcg
            
            if self.VTail2 is not None:
                #if self.VTail2.NoneList.has_key('S'): 
                self.VTail2.L = self.VTail.L
                self.VTail2.L = self.HTail.TE(0*IN) - self.VTail2.Chord(0*IN)*0.75 - Xcg
                #else:
                #    self.VTail2.L = self.VTail.L

            self._PositionTail()
                    
        #
        # The tail must be positioned before the wing
        #
        self._PositionWingX()
        self._PositionPropulsion()
    

#===============================================================================
    def _PositionTail(self):
        """
        Positions the fuselage on the horizontal tail
        """
        super(ACTLenAircraft, self)._PositionTail()        
        
        if self.VTail2 is None:
            return
        
        # Update CG
        self.VTail2.SetAircraftXcg(self.Xcg())

        self.VTail2.X[1] = self.HTail.Tip()[1]*self.VTail2Pos
        
        if self.VTail2Root is None:
            self.VTail2.X[2] = self.VTail.X[2]
        else:
            self.VTail2.X[2] = self.HTail.X[2] - self.VTail2.b*self.VTail2Root
                
        self.VTail2.SetDirty()

#===============================================================================
    def UpdateDrawing(self):
        """
        Updates all drawings of the tail aircraft
        """
        super(ACTLenAircraft, self).UpdateDrawing()
                        
        if self.VTail2 is not None:
            self.VTail2.UpdateDrawing()


#===============================================================================
    def _CDComponents(self, a2dw, del_e, V = None):
        """
        Calculates the total area weighted CD 
        of the aircraft with the given elevator deflection

        TODO: Check if the there should be sin(DHW) here or not
        
        Inputs:
            a2dw  - 2D angle of attack of the main wing
            del_e - Elevator deflection
            V     - Velocity 
        """
        CD, CDLegend = super(ACTLenAircraft, self)._CDComponents(a2dw, del_e, V)
        
        if self.VTail2 is None:
            return CD, CDLegend
        
        VTail2 = self.VTail2

        if V is None:
            V = self.GetV_LO()
               
        if len( self.MakeList(V) ) > 1 and len( self.MakeList(a2dw) ) > 1:
            raise ACAircraftError("Cannot compute drag from a list of angles of attack and a list of velocities")

        Sw = self.Wing.S
        Sv  = VTail2.S
        CDv = VTail2.CD(a2dw*0, 0*ARCDEG) #No rudder deflection for now
        CD.append((CDv*Sv/Sw) )
        CDLegend.append('Vertical Tail 2')

        return CD, CDLegend        
