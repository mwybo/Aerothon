import unittest
import math
import scalar.units as un
import scalar.scalar as sc
import numpy as npy

one = npy.float64(1.)*un.M
sqrt3 = npy.sqrt(3.)*un.M
two = npy.float64(2.)*un.M
three = npy.float64(3.)*un.M
four = npy.float64(4.)*un.M
five = npy.float64(5.)*un.M

class ScalarNumpyTestCases(unittest.TestCase):

    def setUp(self): pass
    def tearDown(self): pass


    def testNumProperty(self):
        self.failUnlessEqual(three.num, 3.)
        self.failUnlessEqual(four.num, 4.)
        self.failUnlessEqual(five.num, 5.)


    def testSqrtFailure(self):
        """Test that sqrt fails rather than silently screws up"""

        self.assertRaises(Exception, npy.sqrt, three**2 + four**2)


    def testHalfPowerPasses(self):
        """Tests to see that **.5 works"""

        soln = (three**2 + four**2)**.5
        self.failUnlessEqual(soln, five)


    def testFloat(self):
        """Checks for failure when trying to convert scalars to floats"""
        self.assertRaises(sc._scalar.InconsistentUnits, float, three)
        self.assertRaises(sc._scalar.InconsistentUnits, float, 42.*un.SLUG*un.FT/un.SEC**2)


    def testTrigFailure(self):
        self.assertRaises(Exception, npy.arctan2, three/2., 1.)


    def testTrigPasses(self):
        basesoln = npy.arctan2(2., 1.)
        mathsoln = math.atan2(2., 1.)
        self.failUnlessEqual(basesoln, mathsoln)

        npysoln = npy.arctan2(two.num, one.num)
        self.failUnlessEqual(basesoln, npysoln)

        soln = sc.atan2(two, one)
        self.failUnlessEqual(npysoln, soln)


if __name__=='__main__':
    unittest.main()
