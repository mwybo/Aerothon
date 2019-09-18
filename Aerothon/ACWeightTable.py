from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import ACBase, Force, Length, MomentOfInertia, g, Unitless
from scalar.units import LBF, SLUG, IN, FT
from scalar.units import AsUnit
from ACTableBase import ACTableBase
import pylab as pyl


import numpy as npy
from matplotlib.colors import colorConverter

################################################################################
def pastel(colour, weight=2.4):
    """ Convert colour into a nice pastel shade"""
    rgb = npy.asarray(colorConverter.to_rgb(colour))
    # scale colour
    maxc = max(rgb)
    if maxc < 1.0 and maxc > 0:
        # scale colour
        scale = 1.0 / maxc
        rgb = rgb * scale
    # now decrease saturation
    total = rgb.sum()
    slack = 0
    for x in rgb:
        slack += 1.0 - x

    # want to increase weight from total to weight
    # pick x s.t.  slack * x == weight - total
    # x = (weight - total) / slack
    x = (weight - total) / slack

    rgb = [c + (x * (1.0-c)) for c in rgb]

    return rgb

################################################################################
def sumstr(aList):
    """
    Used to sum up a list of strings
    """
    sums = ''
    for str in aList:
        sums += str
    return sums

################################################################################
class ACWeightTable(ACBase, ACTableBase):
    """
    This class will summarize all the weights specified in the Aircraft.     
    """
    
    WeightUnit = LBF
    WeightUnitName = 'lbf'

    MOIUnit = SLUG*FT**2
    MOIUnitName = 'slug ft^2'
    
################################################################################
    class CPart(ACBase):
        """
        A class for managing parts of the aircraft
        
        If the weight and CG are not specified they will be computed
        from the sub-parts.
        
        Attributes:
            Order    - An integer value representing the sorting order in decreasing order
            Collapse - Logical if the sub parts should be collapsed and not displayed
            Group    - A group category that the part belongs to
            SortBy   - 'Weight', 'CG', 'Group', or 'Order': A string indicating what to sort the weight table by
            Weight   - The weight of the sub part
            CG       - The center of gravity of the sub part
            MPI      - The moment of inertia of the sub part relative to the aircraft CG
        """
#===============================================================================
        def __init__(self, PartName, Part, Level, AicraftCG):
            super(ACWeightTable.CPart, self).__init__()
            
            self.name = Part.name
            
            self.__dict__['Order']    = 1000000   ; self.UnitList['Order'] = Unitless
            self.__dict__['Collapse'] = False
            self.__dict__['PartName'] = PartName
            self.__dict__['Group']    = Part.WeightGroup
            
            self.__dict__['SortBy'] = 'Order'
            
            self.__dict__['Weight'] = Part.Weight ; self.UnitList['Weight'] = Force
            self.__dict__['CG']     = Part.CG()   ; self.UnitList['CG']     = Length
            
            #
            # Calculate the moment of inertia relative to the aircraft CG
            #
            MOI = Part.MOI()
            dCG = AicraftCG - self.CG
            PM  = self.Weight / g
            MOI[0] += PM*(dCG[1]**2 + dCG[2]**2) #Ixx
            MOI[1] += PM*(dCG[0]**2 + dCG[2]**2) #Iyy
            MOI[2] += PM*(dCG[1]**2 + dCG[0]**2) #Izz
            self.__dict__['MOI'] = MOI  ; self.UnitList['MOI'] = MomentOfInertia

            #
            # A list of sub parts
            #
            self.param.Parts = []
            
            #
            # The recursion level of the parts, determins the number of tabs in the table
            #
            self.param.Level = Level
            self.param.AicraftCG = AicraftCG            
       
#===============================================================================
        def AddPart(self, PartName, Part):
            """
            Adds a part to this part
            
            Input:
                PartName - The name given to the part (not used in the table)
                Part     - A part that is added to this part
            """
            self.param.Parts.append(ACWeightTable.CPart(PartName, Part, self.param.Level + 1, self.param.AicraftCG))
            self.__dict__[PartName] = self.param.Parts[-1]

#===============================================================================
        def _GetParts(self):
            """
            Returns the parts sorted appropriately
            """
            Parts = self.param.Parts
            
            for Part in Parts:
                Part.SortBy = self.SortBy
            Parts.sort()
            
            return Parts

#===============================================================================
        def GetGroups(self, Groups):
            """
            Returns a list of all the group names
            """
            Groups.add(self.Group)

            for Part in self.param.Parts:
                Part.GetGroups(Groups)
    
#===============================================================================
        def GetGroupWeights(self, GroupWeights):
            """
            Computes the total weight of each group
            """
            if self.param.Level != 0:
                GroupWeights[self.Group] += self.Weight
            
            for Part in self.param.Parts:
                Part.GetGroupWeights(GroupWeights)
    
#===============================================================================
        def PrintParts(self):
            """
            Prints all the weights to the screen
            """
            print sumstr(['  ']*self.param.Level) + self.PartName
            
            Parts = self._GetParts()
            for Part in Parts:
                Part.PrintParts()

#===============================================================================
        def AddTableCells(self, rowLables, cellText, cellColours, VectSep):
            """
            Adds the current part to the cells of the weight table
            
            Input:
                rowLables - A list row labels of the table
                cellText  - A list text of the cells 
                VectSep   - String vector value separator
            """
            wAsUnit   = '%1.3f'
            cgAsUnit  = '%1.1f'
            moiAsUnit = '%1.3f'
                        
            rowLables.append(sumstr(['   ']*self.param.Level) + self.name)
            
            Weight = self.Weight / ACWeightTable.WeightUnit
            CG     = self.CG / (IN)
            MOI    = self.MOI / ACWeightTable.MOIUnit
            
            cellText.append([ (wAsUnit % Weight)                      
                            , (cgAsUnit % CG[0]).rjust(7) + VectSep + 
                              (cgAsUnit % CG[1]).rjust(7) + VectSep + 
                              (cgAsUnit % CG[2]).rjust(7)
                            , (moiAsUnit % MOI[0]).rjust(7) + VectSep + 
                              (moiAsUnit % MOI[1]).rjust(7) + VectSep + 
                              (moiAsUnit % MOI[2]).rjust(7)            
                            , self.Group ])

            if self.name == 'Servo':
                pass
            
            if self.param.Level == 0:
                colour = pastel([1., 1., 0.], weight=2.7)
                cellColours.append([ colour  for i in range(len(cellText[0])) ])
            else:
                cellColours.append([ [1,1,1] for i in range(len(cellText[0])) ])
            
            #
            # Add all the parts of this part if it is not collapsed
            #
            if not self.Collapse:
                Parts = self._GetParts()
                for Part in Parts:
                    Part.AddTableCells(rowLables, cellText, cellColours, VectSep)
        
#===============================================================================
        def __getitem__(self, PartName):
            """
            Overloads the [] operator to retrieve a part.
            """
            return self.__dict__[PartName]

#===============================================================================
        def __setitem__(self, PartName, Part):
            """
            Overloads the [] operator to set a part.
            """
            return self.AddPart(PartName, Part)

#===============================================================================
        def _cmp(self):
            """
            Returns the numerical value used for comparison in sorting
            """
            if   self.SortBy.upper()[0] == 'W': #Weight
                return -self.Weight
            elif self.SortBy.upper()[0] == 'C': #CG
                return self.CG[0]
            elif self.SortBy.upper()[0] == 'G': #Group
                return self.Group
            else:                               #Order
                return self.Order

#===============================================================================
# Overloaded comparison operators
        def __le__(self, other):
            return self._cmp() <= other._cmp()
        def __ge__(self, other):
            return self._cmp() >= other._cmp()
        def __lt__(self, other):
            return self._cmp() < other._cmp()        
        def __gt__(self, other):
            return self._cmp() > other._cmp()

################################################################################
# ACWeightTable
################################################################################
    def __init__(self, AircraftCG):
        super(ACWeightTable, self).__init__()
        
        self.param.AircraftCG = AircraftCG
        self.__dict__['SortBy'] = 'Order'
       
        #
        # A list of parts
        #
        self.param.Parts = []
        self.param.AicraftCG = AircraftCG
        
#===============================================================================
    def AddPart(self, PartName, Part):
        """
        Adds a part to the weight table
        
        Input:
            PartName - The name given to the part (not used in the table)
            Part     - A part added to the weight table
        """
        self.param.Parts.append(ACWeightTable.CPart(PartName, Part, 0, self.param.AicraftCG))
        self.__dict__[PartName] = self.param.Parts[-1]

#===============================================================================
    def __getitem__(self, PartName):
        """
        Overloads the [] operator to retrieve a part.
        """
        return self.__dict__[PartName]

#===============================================================================
    def __setitem__(self, PartName, Part):
        """
        Overloads the [] operator to set a part.
        """
        return self.AddPart(PartName, Part)

#===============================================================================
    def _GetParts(self):
        """
        Returns the parts sorted appropriately
        """
        Parts = self.param.Parts
        
        for Part in Parts:
            Part.SortBy = self.SortBy
        Parts.sort()
        
        return Parts

#===============================================================================
    def _GetTableCells(self, VectSep):
        """
        Returns labels and cells for the weight table
        
        Inputs:
            VectSep - String vector value separator
        """
        Parts = self._GetParts()
                
        TotalWeight = 0*LBF
        TotalCG     = npy.array([0,0,0])*IN*LBF
        TotalMOI    = npy.array([0,0,0])*SLUG*FT**2
                               
        colLabels = ['Weight (' + ACWeightTable.WeightUnitName + ')', 'CG [x,y,z] (in)', 'MOI [Ixx, Iyy, Izz] (' + ACWeightTable.MOIUnitName + ')', 'Group']
        rowLabels = []
        cellText  = []
        cellColours = []
        
        for Part in Parts:
            TotalWeight += Part.Weight
            TotalCG     += Part.Weight*Part.CG
            TotalMOI    += Part.MOI
            
            Part.AddTableCells(rowLabels, cellText, cellColours, VectSep)
            
        TotalCG /= TotalWeight
        
        rowLabels.append('Totals:')
        
        wAsUnit   = '%1.2f'
        cgAsUnit  = '%1.1f'
        moiAsUnit = '%1.3f'
                
        TotalWeight = TotalWeight / ACWeightTable.WeightUnit
        TotalCG     = TotalCG / (IN)
        TotalMOI    = TotalMOI / ACWeightTable.MOIUnit
        
        cellText.append([(wAsUnit % TotalWeight), 
                         (cgAsUnit % TotalCG[0]).rjust(7) + VectSep + 
                         (cgAsUnit % TotalCG[1]).rjust(7) + VectSep + 
                         (cgAsUnit % TotalCG[2]).rjust(7),
                         (moiAsUnit % TotalMOI[0]).rjust(7) + VectSep + 
                         (moiAsUnit % TotalMOI[1]).rjust(7) + VectSep + 
                         (moiAsUnit % TotalMOI[2]).rjust(7),
                         ''            ])
        
        cellColours.append([ [1,1,1] for i in range(len(cellText[0])) ])
        
        return colLabels, rowLabels, cellText, cellColours

#===============================================================================
    def GetGroups(self):
        """
        Returns a list of all the group names
        """
        Groups = set()
        for Part in self.param.Parts:
            Part.GetGroups(Groups)

        return list(Groups)

#===============================================================================
    def GetGroupWeights(self):
        """
        Computes the total weight of each group
        """
        Groups = self.GetGroups()
        
        GroupWeights = {}
        for Group in Groups:
            GroupWeights[Group] = 0
            
        for Part in self.param.Parts:
            Part.GetGroupWeights(GroupWeights)
            
        return GroupWeights
        
#===============================================================================
    def PrintParts(self):
        """
        Prints all the weights to the screen
        """
        Parts = self._GetParts()
        for Part in Parts:
            Part.PrintParts()

#===============================================================================
    def WriteWeightTable(self, FileName = 'WeightTable.txt'):
        """
        Writes the weight table to a file
        """
        fTable = open(FileName, 'w')
        
        #
        # Get the labels and cell text
        #
        colLabels, rowLabels, cellText, cellColours = self._GetTableCells(VectSep = ', ')
        
        #
        # Determine the max row label so the rest can be padded
        #
        maxRowWidth = max([len(rowLabels[i]) for i in range(len(rowLabels))])
        
        #
        # Write the column labels
        #
        fTable.write(sumstr([' ']*maxRowWidth) + '\t')
        fTable.write(sumstr([colLabels[i] + '\t' for i in range(len(colLabels))]) + '\n' )
        
        #
        # Write the rest of the table
        #
        for i in range(len(rowLabels)):
            pad = maxRowWidth - len(rowLabels[i])
            fTable.write(rowLabels[i] + sumstr([' ']*pad) + '\t')
            fTable.write(sumstr([cellText[i][j] + '\t' for j in range(len(cellText[i]))]) + '\n')

        fTable.close()

#===============================================================================
    def PlotWeightTable(self):
        """
        Plots the table of weights
        """
        #
        # Get the labels and cell text
        #
        colLabels, rowLabels, cellText, cellColours = self._GetTableCells(VectSep = ', ')
        
        rowColours = []
        for row in cellColours:
            rowColours.append(row[0])
        
        colWidths = self._GetTableColWidths(rowLabels, cellText, colLabels)
        
        CurrentAxis = pyl.axes([0.07, 0.1, 0.9, 0.9])
        CurrentAxis.autoscale_view(tight=False, scalex=True, scaley=False)
        pyl.axis('off')
        table = pyl.table(cellText    = cellText , cellLoc = 'center',
                          rowLabels   = rowLabels, rowLoc  = 'left'  , rowColours = rowColours, 
                          colLabels   = colLabels, colLoc  = 'center',
                          colWidths   = colWidths,
                          cellColours = cellColours,
                          loc = 'upper right')
        
#        table.auto_set_font_size(False)
#        table.set_fontsize(11)
        

################################################################################
if __name__ == '__main__':
    
    from TestAircraft import aircraft
    
#    WeightTable = ACWeightTable([0*IN, 0*IN, 0*IN])
#    
#    WeightTable.AddPart('Wing', aircraft.Wing)
#    WeightTable.Wing.AddPart('Wing', aircraft.Wing)
#    
#    WeightTable['Wing2'] = aircraft.Wing
#    
#    WeightTable.AddPart('HTail', aircraft.HTail)
#    WeightTable.Wing.Wing.AddPart('HTail', aircraft.HTail)
#    
#    aircraft.Fuselage.AddToWeightTable('Fuselage', WeightTable)
    
    WeightTable = aircraft.GetWeightTable()
    
#    WeightTable.Wing.Collapse = True
#    WeightTable.Wing.Order = 0

    WeightTable.PlotWeightTable()
    WeightTable.WriteWeightTable()
    WeightTable.PrintParts()
   
    pyl.show()
    print 'Weight Table Test Completed'