from ACBatchRun import ACBatchRun
import os
import sys

################################################################################
class ACAVL(ACBatchRun):
    """
    A class for writing, running and reading AVL files.
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
            self.__dict__['Exec'] = SoftwareDir + 'avl.exe'
        elif sys.platform.startswith('darwin'):
            self.__dict__['Exec'] = SoftwareDir + 'avl.mac'
        else:
            raise "An executable for AVL is not available for this Operating system. Add it to the 'Software' folder and modify 'ACAVL.py'"
            
        self.__dict__['CommandFile'] = 'AVLCommands.dat'
        
#===============================================================================
    class AVLRun(ACBatchRun.BatchRun):
        """
        
        """

#===============================================================================
        def __init__(self, InputFile, RunDir, CommandFile):
            """
            
            """
            ACBatchRun.BatchRun.__init__(self, InputFile, RunDir, CommandFile)
            self.AddCommand('oper')
            
#===============================================================================
        def DumpStability(self, FileName):
            """
            Dumps a stability file from AVL.
            """
            
            self.STFiles.append(self.RunDir + FileName)
            self.AddCommand('x')
            self.AddCommand('st')
            self.AddCommand(self.RunDir + FileName)

#===============================================================================
        def Exit(self):
            """
            Exits
            """
            self.AddCommand('\n\n\n')
            self.AddCommand('q')           
                                
#===============================================================================
    def AddRun(self, RunName, InputFile, CommandFile = None):
        """
        
        """
        
        if CommandFile is None:
            CmdFile = self.RunDir + self.CommandFile
        else:
            CmdFile = self.RunDir + CommandFile
            
        self.Runs.append(ACAVL.AVLRun(self.RunDir + InputFile, self.RunDir, CmdFile))
        self.__dict__[RunName] = self.Runs[-1]
        
#===============================================================================
    def ExecuteAVL(self):
        """
        
        """
        ACBatchRun.Execute(self)

#===============================================================================
    def ReadAVLFiles(self, CoeffClass, **args):
        """
        
        """
        
        Coeffs = []
        
        for run in self.Runs:
            for File in run.STFiles:
                Coeffs.append(CoeffClass(**args))
                self._ReadAVLOutput(Coeffs[-1], File)
                
        return Coeffs   

#===============================================================================
    def _ReadAVLOutput(self, Derivatives, AVLOutputFile):
        """
        Reads the AVL output file

        Inputs
            Derivatives - ACDerivaties class with a list of derivatives to read in
        """
        for key in Derivatives._Coefficients:
            try:

                f = open(AVLOutputFile,'r')
                line = f.readline()
                #
                # Find the line containing the derivative
                #
                while line.find(key + ' =') == -1:
                    line = f.readline()
                    if len(line) == 0: #This means it hit the end of the file and we have a programming error as key is not in the file
                        raise "Could not find " + key + " in " + AVLOutputFile
                               
                line  = line[line.find(key):]
                line  = line[line.find('=')+1:]
                line  = line.lstrip(' ')
                line  = line[:line.find(' ')]
                Derivatives.__dict__[key] = float(line)
                f.close()
            except IOError, (errno, strerror):
                print "While reading '", AVLOutputFile, "'"
                raise
                    