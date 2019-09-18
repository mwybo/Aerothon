""""
A collection of unit tests to guarantee that the ACWing class never becomes corrupted
"""

import unittest
from unum.units import FT, IN, SEC, LBF, SEC, ARCDEG
from Aerothon.ACWing import ACMainWing

################################################################################
class ACMainWingAreaTest(unittest.TestCase):
    """
    A set of tests design to test the various ways to calculate that are of a main wing.
    """
#===============================================================================
    def setUp(self):
        """
        This will setup a main wing with some default parameters. The remaining tests will calculate
        the area and other related properties in numerous ways.
        
        The two equations are:
        
        S = f(Lift_LO, V_Stall)
        S = f(AR, b)
        
        If S is specified, on variable from each of the two equations must be specified.
        If S is not specified, three other variables must be specified.
        
        """
        Wing = self.Wing = ACMainWing(1)
        Wing.Airfoil = ('../../Airfoils/S1223/S1223.dat','../../Airfoils/S1223')
        
        Wing.TR = [1, 0.7, 0.7]
        Wing.Fb = [0.3, 0.8, 1]
        Wing.Gam = [0.0 * ARCDEG, 0.0 * ARCDEG, 0.0 * ARCDEG]
        Wing.CEdge = 'LE'
        Wing.ConstUpper = True
        Wing.FullWing = True
        Wing.Alt_LO   = 600 * FT
        
        #
        # Setup consistent numbers for area calculations
        #
        self.Lift_LO = 27.00654 * LBF
        self.V_Stall = 41.2699923511 * FT/SEC
        self.AR      = 16
        self.b       = 10*FT
        self.S       = 900.*IN**2
    
#===============================================================================
    def ClearWing(self):
        Wing = self.Wing
        Wing.Lift_LO = None
        Wing.V_Stall = None
        Wing.S       = None
        Wing.AR      = None
        Wing.b       = None

#===============================================================================
    def test_S_and_V_Stall(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.Lift_LO = self.Lift_LO
        Wing.AR      = self.AR
        Wing.b       = self.b
        self.assertAlmostEqual(Wing.S.asNumber(IN**2), self.S.asNumber(IN**2), 0, 
                               "Filed to calculate Area from Lift of lift (Lift_LO), Aspectratio (AR) and Span (b)")

        self.assertAlmostEqual(Wing.V_Stall.asNumber(FT/SEC), self.V_Stall.asNumber(FT/SEC), 0, 
                               "Filed to calculate V_Stall from Lift of lift (Lift_LO), Aspectratio (AR) and Span (b)")

#===============================================================================
    def test_S_and_Lift_LO(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.V_Stall = self.V_Stall
        Wing.AR      = self.AR
        Wing.b       = self.b
        self.assertAlmostEqual(Wing.S.asNumber(IN**2), self.S.asNumber(IN**2), 0, 
                               "Filed to calculate Area from Stall Speed (V_Stall), Aspectratio (AR) and Span (b)")

        self.assertAlmostEqual(Wing.Lift_LO.asNumber(LBF), self.Lift_LO.asNumber(LBF), 0, 
                               "Filed to calculate Lift off Lift from Stall Speed (V_Stall), Aspectratio (AR) and Span (b)")

#===============================================================================
    def test_S_and_b(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.Lift_LO = self.Lift_LO
        Wing.V_Stall = self.V_Stall
        Wing.AR      = self.AR
        self.assertAlmostEqual(Wing.S.asNumber(IN**2), 900.7, 0, 
                               "Filed to calculate Area from Aspectratio (AR), Lift off Lift (Lift_LO) and Stall speed (V_Stall)")

        self.assertAlmostEqual(Wing.b.asNumber(FT), self.b.asNumber(FT), 0, 
                               "Filed to calculate span from Aspectratio (AR), Lift off Lift (Lift_LO) and Stall speed (V_Stall)")

#===============================================================================
    def test_S_and_AR(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.Lift_LO = self.Lift_LO
        Wing.V_Stall = self.V_Stall
        Wing.b       = self.b
        self.assertAlmostEqual(Wing.S.asNumber(IN**2), 901.5, 0, 
                               "Filed to calculate Area from Span (b), Lift off Lift (Lift_LO) and Stall speed (V_Stall)")

        self.assertAlmostEqual(Wing.AR, self.AR, 0, 
                               "Filed to calculate Aspectratio from Span (b), Lift off Lift (Lift_LO) and Stall speed (V_Stall)")

#===============================================================================
    def test_b_and_Lift_LO(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.Lift_LO = self.Lift_LO
        Wing.S       = self.S
        Wing.AR      = self.AR
        self.assertAlmostEqual(Wing.b.asNumber(FT), self.b.asNumber(FT), 0, 
                               "Filed to calculate span from Lift of Lift (Lift_LO), Area (S) and Aspectratio (AR)")

        self.assertAlmostEqual(Wing.V_Stall.asNumber(FT/SEC), self.V_Stall.asNumber(FT/SEC), 0, 
                               "Filed to calculate span from Lift of Lift (Lift_LO), Area (S) and Aspectratio (AR)")

#===============================================================================
    def test_b_and_V_Stall(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.V_Stall = self.V_Stall
        Wing.S       = self.S
        Wing.AR      = self.AR
        self.assertAlmostEqual(Wing.b.asNumber(FT), self.b.asNumber(FT), 0, 
                               "Filed to calculate span from Stall speed (V_Stall), Area (S) and Aspectratio (AR)")

        self.assertAlmostEqual(Wing.Lift_LO.asNumber(LBF), self.Lift_LO.asNumber(LBF), 0, 
                               "Filed to calculate Lift of Lift from Stall speed (V_Stall), Area (S) and Aspectratio (AR)")
#===============================================================================
    def test_AR_and_Lift_LO(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.Lift_LO = self.Lift_LO
        Wing.S       = self.S
        Wing.b       = self.b
        self.assertAlmostEqual(Wing.AR, self.AR, 0, 
                               "Filed to calculate span from Lift of Lift (Lift_LO), Area (S) and Aspectratio (AR)")

        self.assertAlmostEqual(Wing.V_Stall.asNumber(FT/SEC), self.V_Stall.asNumber(FT/SEC), 0, 
                               "Filed to calculate span from Lift of Lift (Lift_LO), Area (S) and Aspectratio (AR)")

#===============================================================================
    def test_AR_and_V_Stall(self):
        self.ClearWing()
        Wing = self.Wing
        Wing.V_Stall = self.V_Stall
        Wing.S       = self.S
        Wing.b       = self.b
        self.assertAlmostEqual(Wing.AR, self.AR, 0, 
                               "Filed to calculate Aspectratio from Stall speed (V_Stall), Area (S) and Aspectratio (AR)")

        self.assertAlmostEqual(Wing.Lift_LO.asNumber(LBF), self.Lift_LO.asNumber(LBF), 0, 
                               "Filed to calculate Lift of Lift from Stall speed (V_Stall), Area (S) and Aspectratio (AR)")

################################################################################
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ACMainWingAreaTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

