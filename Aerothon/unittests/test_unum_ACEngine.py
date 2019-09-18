""""
A collection of unit tests to guarantee that the ACEngine class is functioning properly
"""

import unittest
from unum.units import FT, IN, SEC, LBF, SEC, LBM, PSFC, RPM, HP, PSI, HR
from Aerothon.ACEngine import ACEngine
from Aerothon.ACBase import g

################################################################################
class ACEngineTest(unittest.TestCase):
    """
    
    """
#===============================================================================
    def test_Engine1(self):
        """
        Tests that the engine class is working properly 
        """
        # Create the engine defined in the mathcad pdf documents
        Engine  = ACEngine()
        Engine.Rbs          = 1.1
        Engine.Rla          = 3.5
        Engine.NumCyl       = 1
        Engine.NumRev       = 1
        Engine.CompRatio    = 9
        Engine.Vd           = 0.607*IN**3
        Engine.PistonSpeedR = 38.27*FT/SEC
        Engine.BMEPR        = 50.038*LBF/IN**2
        Engine.Rnmt         = 0.01
        Engine.Rtmt         = 1.5
        Engine.MEPtlmt      = 10.1526*LBF/IN**2
        Engine.SFCmt        = 1*PSFC
        Engine.Rnsfc        = 0.8
        Engine.A_F          = 16
        Engine.PS           = 1

        Nr = Engine.Nr()
        
        self.assertAlmostEqual(Engine.B().asNumber(IN), 0.9473, 4)
        self.assertAlmostEqual(Engine.L().asNumber(IN), 0.8612, 4)

        self.assertAlmostEqual(Nr.asNumber(RPM)                            , 15997.7 , 1, "Rated crankshaft rotational speed (Nr)")
        self.assertAlmostEqual(Engine.Tr().asNumber(FT*LBF)                , 0.4028  , 4, "Rated crankshaft rotational speed (Nr)")
        self.assertAlmostEqual(Engine.Pr().asNumber(HP)                    , 1.2270  , 4, "Brake power at Nr")
        self.assertAlmostEqual(Engine.Nmt().asNumber(RPM)                  , 159.9766, 4, "RPM at max torque")
        self.assertAlmostEqual(Engine.Tm().asNumber(FT*LBF)                , 0.6043  , 4, "Max Torque")
        self.assertAlmostEqual(Engine.BMEPMT().asNumber(PSI)               , 75.057  , 3, "Max brake mean effective pressure")
        self.assertAlmostEqual(Engine.BMEPsl(Nr).asNumber(PSI)             , 50.038  , 3, "Quadratic fit of BMEP at sea level at Nr")
        self.assertAlmostEqual(Engine.IMEPMT().asNumber(PSI)               , 85.2096 , 4, "Max indicated mean effective pressure")
        self.assertAlmostEqual(Engine.MEPtlsl(Nr).asNumber(PSI)            , 35.1716 , 4, "Assumed total losses")
        self.assertAlmostEqual(Engine.MEPtf(Nr).asNumber(PSI)              , 16.8122 , 4, "Total friction losses")
        self.assertAlmostEqual(Engine.MEPtp(0*FT, Nr).asNumber(PSI)        , 18.3594 , 4, "Total pumping losses")
        self.assertAlmostEqual(Engine.IMEP(0*FT).asNumber(PSI)             , 85.2096 , 4, "Indicated (gas) mean effective pressure")
        self.assertAlmostEqual(Engine.BMEP(0*FT, Nr).asNumber(PSI)         , 50.038  , 3, "Brake mean effective pressure")
        self.assertAlmostEqual(Engine.Ti(0*FT).asNumber(FT*LBF)            , 0.686   , 3, "Indicated torque")
        self.assertAlmostEqual(Engine.Pi(0*FT, Nr).asNumber(HP)            , 2.0895  , 4, "Indicated power")
        self.assertAlmostEqual(Engine.Tb(0*FT, Nr).asNumber(FT*LBF)        , 0.4028  , 4, "Brake torque")
        self.assertAlmostEqual(Engine.Pb(0*FT, Nr).asNumber(HP)            , 1.227   , 3, "Brake power")
        
        self.assertAlmostEqual(Engine.Mdotfr().asNumber(LBM/HR)            , 1.1599  , 4, "Rated fuel flow at minimum sfc")
        self.assertAlmostEqual(Engine.Mdotf(0*FT, Nr).asNumber(LBM/HR)     , 1.4499  , 4, "Fuel flow (mass/hr)")
        self.assertAlmostEqual(Engine.Mdota(0*FT, Nr).asNumber(LBM/HR)     , 23.1977 , 4, "Air flow (mass/hr)")
        self.assertAlmostEqual(Engine.FuelFlow(0*FT, Nr).asNumber(FT**3/HR), 0.0323  , 4, "Fuel flow (volume/hr)")
        self.assertAlmostEqual(Engine.SFC(0*FT, Nr).asNumber(PSFC)         , 1.1816  , 4, "Specific fuel consumption")

#        print ''
#        print 'Rated crankshaft rotational speed (Nr)   = ',Nr.asUnit(RPM)
#        print 'Torque at Nr                             = ',Engine.Tr().asUnit(FT*LBF)
#        print 'Brake power at Nr                        = ',Engine.Pr().asUnit(HP)
#        print 'RPM at max torque                        = ',Engine.Nmt().asUnit(RPM)
#        print 'Max Torque                               = ',Engine.Tm().asUnit(FT*LBF)
#        print 'Max brake mean effective pressure        = ',Engine.BMEPMT().asUnit(PSI)
#        print 'Quadratic fit of BMEP at sea level at Nr = ',Engine.BMEPsl(Nr)
#        print 'Max indicated mean effective pressure    = ',Engine.IMEPMT().asUnit(PSI)
#        print 'Assumed total losses                     = ',Engine.MEPtlsl(Nr).asUnit(PSI)
#        print 'Total friction losses                    = ',Engine.MEPtf(Nr).asUnit(PSI)
#        print 'Total pumping losses                     = ',Engine.MEPtp(0*FT, Nr).asUnit(PSI)
#        print 'Indicated (gas) mean effective pressure  = ',Engine.IMEP(0*FT).asUnit(PSI)
#        print 'Brake mean effective pressure            = ',Engine.BMEP(0*FT, Nr).asUnit(PSI)
#        print ''
#        print 'Indicated torque = ',Engine.Ti(0*FT).asUnit(FT*LBF)
#        print 'Indicated power  = ',Engine.Pi(0*FT, Nr).asUnit(HP)
#        print 'Brake torque     = ',Engine.Tb(0*FT, Nr).asUnit(FT*LBF)
#        print 'Brake power      = ',Engine.Pb(0*FT, Nr).asUnit(HP)
#        print ''
#        print 'Rated fuel flow at minimum sfc = ',Engine.Mdotfr().asUnit(LBM/HR)
#        print 'Fuel flow (mass/hr)            = ',Engine.Mdotf(0*FT, Nr).asUnit(LBM/HR)
#        print 'Air flow  (mass/hr)            = ',Engine.Mdota(0*FT, Nr).asUnit(LBM/HR)
#        print 'Fuel flow (volume/hr)          = ',Engine.FuelFlow(0*FT, Nr).asUnit(FT**3/HR)
#        print 'Specific fuel consumption      = ',Engine.SFC(0*FT, Nr).asUnit(PSFC)
        

################################################################################
if __name__ == '__main__':    
    suite = unittest.TestLoader().loadTestsFromTestCase(ACEngineTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
