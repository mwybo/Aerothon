import numpy as npy

class ACTableBase:
    """
    A class for useful table functions
    """
#===============================================================================
    def _GetTableColWidths(self, rowLabels, cellText, colLabels = None): 
        """
        Calculates column widths for a table
        
        Inputs:
            rowLabels - The row lables
            cellText  - The cell text
            colLabels - The column lables
        """                        
        #
        # Compute some 'nice' column widths for the table
        #
        rlableWidth = 0
        if colLabels is not None:
            colWidths = [len(colLabels[i]) for i in range(len(colLabels))]
        else:
            colWidths = [len(cellText[0]) for i in range(len(cellText[0]))]
            
        for ri in range(len(cellText)):
            row = cellText[ri]
            r = []
            for rr in row:
                while rr.find('$') != -1: #Remove any latex code from the length of the string
                    r1 = rr[:rr.find('$')]
                    rr = r1 + rr[rr.find('$',len(r1)+1):]
                    rr = rr.replace( '$', '#', 1) #Replace just one '$' with '#' 
                r.append(rr)
            
            rlableWidth = max(len(rowLabels[ri]),rlableWidth)
            colWidths = npy.array([float(max(colWidths[i],len(r[i]))) for i in range(len(r))])
        colWidths /= sum(colWidths) + rlableWidth
        
        return colWidths    