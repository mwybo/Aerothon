"""
A class for performing AVL calculations for a lifting surface
"""
from ACAVL import ACAVL

################################################################################
class ACAVLWing():
    """
    A class for some general AVL operations for lifting surfaces.
    
    This class is used for inheritance purposes only
    """

#===============================================================================
    def GetAVLCL(self, Alpha, CaseName, RunDir = '', Execute = True):
        """
        Sets up and runs an AVL run for a range of angles of attack
        
        Inputs:
            Alpha   - a range of angle of attacks in degrees (as floats, not scalars)
            CaseName- the case name string of the output files
            RunDir  - the directory in which to output the AVL output files
            Execute - True or False for executing the AVL run
        """
        
        if self.dirty: self.Refresh()
        
        AVL = ACAVL()
        AVL.RunDir = RunDir

        filename = 'AVL' + self.name + '.txt'
        
        if Execute:
            self.WriteAVLWing(RunDir + filename)
        
        AVL.AddRun('AlphaRun', filename)
        
        for a in Alpha:           
            AVL.AlphaRun.AddCommand('a a %1.1f' % a)
            AVL.AlphaRun.DumpStability('BiWing_' + CaseName + '_%1.1f.txt' % a)
            
        AVL.AlphaRun.Exit()
        
        class CLClass:
            _Coefficients = ['CLtot']
            
        if Execute:
            AVL.ExecuteAVL()
        
        CLCs = AVL.ReadAVLFiles(CLClass)
        
        CL = []
        for cl in CLCs:
            CL.append(cl.CLtot)
            
        return CL
    
#===============================================================================
    def GetAVLOswald(self, AVLOutFile, RunDir = '', Execute = True):
        """
        Sets up and runs an AVL run to find the Oswald efficiency for the current wing configuration
        
        Inputs:
            AVLoutfile - The name of the AVL st output
        """
        
        if self.dirty: self.Refresh()
        
        AVL = ACAVL()
        AVL.RunDir = RunDir
        
        filename = 'AVL' + self.name + '.txt'
        
        if Execute:
            self.WriteAVLWing(RunDir + filename)
        AVL.AddRun('O_effRun', filename)
        
        AVL.O_effRun.AddCommand('a a 0')
        AVL.O_effRun.DumpStability(AVLOutFile)
            
        AVL.O_effRun.Exit()
        
        class OeffClass:
            _Coefficients = ['e']
            
        if Execute:
            AVL.ExecuteAVL()
        oeff = AVL.ReadAVLFiles(OeffClass)
            
        return oeff[0].e