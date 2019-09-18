""""
A collection of unit tests to guarantee that the ACMaterial class is functioning properly
"""

import unittest
from unum.units import FT, IN, SEC, LBF, SEC, LBM
from Aerothon.ACMaterial import ACMaterial
from Aerothon.ACBase import g

################################################################################
class ACMaterialTest(unittest.TestCase):
    """
    
    """
#===============================================================================
    def test_Density1(self):
        """
        Calculates ForceDensity from Density
        """
        Mat = ACMaterial()

        #
        # Simply convert force density into density
        #
        Mat.ForceDensity = LBF/IN**3

        Density = (LBF/g/IN**3).asNumber(LBM/IN**3)
        
        self.assertAlmostEqual(Mat.Density.asNumber(LBM/IN**3), Density, 5, 
                               "Filed to calculate Density from ForceDensity")
        
#===============================================================================
    def test_Density2(self):
        """
        Calculates ForceDensity from AreaDensity and Thickness
        """
        Mat = ACMaterial()
        
        Mat.AreaForceDensity = LBF/IN**2
        Mat.Thickness   = IN
        
        Density = (LBF/g/IN**3).asNumber(LBM/IN**3)

        self.assertAlmostEqual(Mat.Density.asNumber(LBM/IN**3), Density, 5, 
                               "Filed to calculate Density from AreaForceDensity and Thickness")

#===============================================================================
    def test_ForceDensity1(self):
        """
        Calculates ForceDensity from Density
        """
        Mat = ACMaterial()

        #
        # Simply convert density into force density
        #
        Mat.Density = LBM/IN**3

        ForceDensity = (LBM*g/IN**3).asNumber(LBF/IN**3)
        
        self.assertAlmostEqual(Mat.ForceDensity.asNumber(LBF/IN**3), ForceDensity, 5, 
                               "Filed to calculate ForceDensity from Density")
        
#===============================================================================
    def test_ForceDensity2(self):
        """
        Calculates ForceDensity from AreaDensity and Thickness
        """
        Mat = ACMaterial()
        
        Mat.AreaDensity = LBM/IN**2
        Mat.Thickness   = IN
        
        ForceDensity = (LBM*g/IN**3).asNumber(LBF/IN**3)

        self.assertAlmostEqual(Mat.ForceDensity.asNumber(LBF/IN**3), ForceDensity, 5, 
                               "Filed to calculate ForceDensity from AreaDensity and Thickness")

################################################################################
if __name__ == '__main__':    
    suite = unittest.TestLoader().loadTestsFromTestCase(ACMaterialTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
