"""
A class for writing, running and reading AVL files.
"""

import os
import subprocess as subp


################################################################################
class ACBatchRun():
    """
    A class for writing, running and reading batch files.
    """

#===============================================================================        
    def __init__(self):
        """
        
        """
        
        #
        # The software directory is one folder up 
        #
        path = __file__
        path = path[:path.rfind(os.sep)] #Remove this file name
        path = path[:path.rfind(os.sep)] #Remove the Aerothon folder
        
        SoftwareDir = path + os.sep + 'Software' + os.sep
        
        self.__dict__['SoftwareDir'] = SoftwareDir
        self.__dict__['Exec'] = SoftwareDir
        self.__dict__['Runs']   = []
        self.__dict__['RunDir'] = ''
        self.__dict__['CommandFile'] = 'Commands.dat'

#===============================================================================
    class BatchRun():
        """
        
        """

#===============================================================================
        def __init__(self, InputFile, RunDir, CommandFile):
            """
            
            """
            self.Commands    = []
            self.STFiles     = []
            self.InputFile   = InputFile
            self.CommandFile = CommandFile
            self.RunDir      = RunDir
            self.ScreenOutputFile = None

#===============================================================================
        def WriteCommandFile(self):
            """
            
            """
            cFile = open(self.CommandFile,'w')
            
            for command in self.Commands:
                cFile.write(command + '\n')
 
            cFile.write('\n')
            cFile.write('\n')
            cFile.write('quit\n')         
           
            cFile.close()

#===============================================================================
        def AddCommand(self, Command):
            """
            Adds a command to the list for AVL to run.
            """
            self.Commands.append(Command)
                                    
#===============================================================================
    def AddRun(self, RunName, InputFile, CommandFile = None):
        """
        
        """
        
        if CommandFile is None:
            CmdFile = self.RunDir + self.CommandFile
        else:
            CmdFile = self.RunDir + CommandFile
            
        self.Runs.append(ACBatchRun.BatchRun(self.RunDir + InputFile, self.RunDir, CmdFile))
        self.__dict__[RunName] = self.Runs[-1]

#===============================================================================
    def Execute(self):
        """
        
        """
        
        for run in self.Runs:
            for File in run.STFiles:
                if os.access(File, os.F_OK): #Check that the file exists before deleting it
                    os.remove(File)
            run.WriteCommandFile()
            
            stdout = open(run.ScreenOutputFile, 'w') if run.ScreenOutputFile is not None else None
            cmd = open(run.CommandFile,'r')
            
            subp.check_call([self.Exec, run.InputFile], stdin = cmd, stdout = stdout)
            
            cmd.close()
            if stdout is not None:
                stdout.close()


            
#===============================================================================

                