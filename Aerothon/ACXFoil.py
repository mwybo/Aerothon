from ACBatchRun import ACBatchRun
from scalar.units import ARCDEG, IN
import numpy as npy
import os
import sys


################################################################################
class ACXFoil(ACBatchRun):
    """
    A class for writing, running and reading xfoil files.
    """

#===============================================================================        
    def __init__(self):
        """
        
        """
        ACBatchRun.__init__(self)
        #
        # The software directory is one folder up 
        #
        path = __file__
        path = path[:path.rfind(os.sep)] #Remove this file name
        path = path[:path.rfind(os.sep)] #Remove the Aerothon folder
        
        SoftwareDir = path + os.sep + 'Software' + os.sep
        
        if sys.platform.startswith('win') or sys.platform.startswith('cygwin'):
            self.__dict__['Exec'] = SoftwareDir + 'xfoil.exe'
        elif sys.platform.startswith('darwin'):
            self.__dict__['Exec'] = SoftwareDir + 'xfoil.mac'
        else:
            raise "An executable for XFoil is not available for this Operating system. Add it to the 'Software' folder and modify 'ACXFoil.py'"
                 
        self.__dict__['CommandFile'] = 'XFoilCommands.dat'
        
#===============================================================================
    class XFoilRun(ACBatchRun.BatchRun):
        """
        
        """

#===============================================================================
        def __init__(self, AirfoilFile, ScreenOutputFile, RunDir, CommandFile):
            """
            
            """
            ACBatchRun.BatchRun.__init__(self, AirfoilFile, RunDir, CommandFile)
            self.ScreenOutputFile = ScreenOutputFile
            
#===============================================================================
        def SetFlap(self, Fc, Ht, da):
            """
            Add the commands to set the flap location
            
            Inputs
                Fc - Fraction of chord of the flap
                Ht - Fraction of thickness location of the hinge (0 - bottom, 0.5 - mid, 1 - top)
                da - Flap deflection in degrees
            """
            
            self.AddCommand('\n\n\n\n')     #Get back to the root directory
            self.AddCommand('gdes')         #Enter the airfoil design
            self.AddCommand('flap')         #Set the flap
            self.AddCommand(str(1-Fc))      #Set the hinge x-location
            self.AddCommand("999")          #Desire hinge thickness
            self.AddCommand(str(Ht))        #Set the hinge y-location
            self.AddCommand(str(da/ARCDEG)) #Set the hinge y-location
            self.AddCommand('x')            #Set the buffered airfoil as the current airfoil
            self.AddCommand('\n')           #Return to the main menu

#===============================================================================
        def ToggleVisc(self, Re):
            """
            Toggle viscosity and sets the Reynolds number
            
            Inputs
                Re - Reynolds number
            """
            
            self.AddCommand('\n\n\n\n')        #Get back to the root directory
            self.AddCommand('oper')            #Enter oper
            self.AddCommand('visc ' + str(Re)) #Set the Re
            self.AddCommand('\n')              #Return to the main menu

#===============================================================================
        def AddHingeMoment(self, Alpha2d, Re):
            """
            Add the commands to calculate a hinge moment
            
            Inputs
                Alpha2d - 2D airfoil Angle of attack
                Re - Reynolds number
            """
            
            self.AddCommand('\n\n\n\n')                 #Get back to the root directory
            self.AddCommand('oper')                     #Enter oper
            self.AddCommand('Re ' + str(Re))            #Set the Re
            self.AddCommand('iter 800')                 #Set the max iterations
            self.AddCommand('init')                     #Initialize BL's
            self.AddCommand('a ' + str(Alpha2d/ARCDEG)) #Compute the angle of attack
            self.AddCommand('fmom')                     #Display the hinge moment
            self.AddCommand('\n')                       #Return to the main menu

#===============================================================================
        def Exit(self):
            """
            Exits
            """
            self.AddCommand('\n\n\n\n\n\n')
            self.AddCommand('quit')           
                                
#===============================================================================
    def AddRun(self, RunName, AirfoilFile, ScreenOutputFile, CommandFile = None):
        """
        Creates a run for executing xfoil and dumping the results into a single file
        """
        
        if CommandFile is None:
            CmdFile = self.RunDir + self.CommandFile
        else:
            CmdFile = self.RunDir + CommandFile
            
        self.Runs.append(ACXFoil.XFoilRun(AirfoilFile, self.RunDir + ScreenOutputFile, self.RunDir, CmdFile))
        self.__dict__[RunName] = self.Runs[-1]
        
#===============================================================================
    def ExecuteXFoil(self):
        """
        
        """
        ACBatchRun.Execute(self)

#===============================================================================
    def ReadXFoilFiles(self):
        """
        
        """
        
        HingeMoments = []
        
        for run in self.Runs:
            File = run.ScreenOutputFile
            HingeMoments.append( self._ReadHingeMoment(File) )
                
        return HingeMoments   

#===============================================================================
    def _ReadHingeMoment(self, OutputFile):
        """
        Reads the XFoil output file

        Inputs
            Output file
        """
        Chm = []
        key = 'Hinge moment/span'
        try:
            f = open(OutputFile,'r')
            line = f.readline()
            #
            # Find the line containing the derivative
            #
            while True:
                while line.find(key + ' =') == -1:
                    line = f.readline()
                    if len(line) == 0: #This means it hit the end of the file 
                        f.close()
                        return npy.array(Chm)
    
                line  = line[line.find(key):]
                line  = line[line.find('=')+1:]
                line  = line.lstrip(' ')
                line  = line[:line.find(' ')]
                Chm.append( float(line) )
            
        except IOError, (errno, strerror):
            print "While reading '", OutputFile, "'"
            raise
        
        