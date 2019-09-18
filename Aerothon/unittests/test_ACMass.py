""""
A collection of unit tests to guarantee that the ACMass class is functioning properly
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
import unittest
import pylab as pyl
from scalar.units import FT, SEC, LBF, ARCDEG, RPM, HP, IN
from Aerothon.ACMass import ACMassCyl
from Aerothon.ACBase import g

################################################################################
class ACMassTest(unittest.TestCase):
    """
    
    """
        

################################################################################
if __name__ == '__main__':    
    suite = unittest.TestLoader().loadTestsFromTestCase(ACMassTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

    Cylinder = ACMassCyl()
    
    Cylinder.LenDi = (2*IN, 1*IN)
    Cylinder.Xat   = (1, 0, 0)
    
    Cylinder.Draw(1)
    pyl.show()