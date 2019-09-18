""""
A collection of unit tests to guarantee that the ACPropeller class is functioning properly
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
import unittest
from scalar.units import FT, SEC, LBF, ARCDEG, RPM, HP, IN
from Aerothon.ACPropeller import ACPropeller
from Aerothon.ACBase import g

################################################################################
class ACPropellerTest(unittest.TestCase):
    """
    
    """
#===============================================================================
    def test_Power(self):
        """
        Tests that the engine class is working properly 
        """
        # Create the propeller defined in the mathcad pdf documents
        Prop = ACPropeller()
        Prop.D          = 14.2*IN
        Prop.PitchAngle = 12*ARCDEG
    #    Prop.Pitch      = 7.1117*IN 
        
        Prop.dAlpha     = 0*ARCDEG
        Prop.Solidity   = 0.0136
        Prop.RD         = 3/8
        Prop.AlphaStall = 14*ARCDEG
        N = 11400 * RPM
 
        self.assertAlmostEqual(Prop.P(N, 0*FT/SEC, 0*FT) / HP, 1.2950, 4)
        

################################################################################
if __name__ == '__main__':    
    suite = unittest.TestLoader().loadTestsFromTestCase(ACPropellerTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
