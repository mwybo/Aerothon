from __future__ import division # let 5/2 = 2.5 rather than 2
from scalar.units import LBF, FT, IN, M, SEC, KG, LBM, PSI, OZM, OZF, GRAM
from scalar.units import AsUnit
from ACMaterial import ACMaterial

#
# Composites
#

CarbonTube = ACMaterial()
CarbonTube.LinearForceDensity = 0.068 * LBF/FT
#
CarbonBar = ACMaterial()
CarbonBar.ForceDensity = 61.44 * LBF/FT**3
#
CarbonFabric = ACMaterial()
CarbonFabric.ForceDensity = 93.83 * LBF/FT**3 #Carbon cloth 0.01" thick with 3 times weight in resin

#
# Woods
#

Balsa = ACMaterial()
Balsa.ForceDensity = 12.96* LBF/FT**3 #Measured by Tommy Marks, Amanda McGee 9/2012
#previously 8.22* LBF/FT**3 #Measured by Seezan Prejapati 11/2011
#Previously 6.18* LBF/FT**3 #Measured by Matt Finke, Prakhar Aghamkar
#previously 9.216 * LBF/FT**3  #Measured 1/15/09, Bob Klenke, Robert Knapke
Balsa.E = 500000 * PSI

# source <http://www.engineeringtoolbox.com/wood-density-d_40.html>
Basswood = ACMaterial()
Basswood.ForceDensity = 25.92* LBF/FT**3 #Measured by Tommy Marks, Amanda McGee 9/2012
#previously 30.74* LBF/FT**3 #measured by Seezan Prejapati

Poplar = ACMaterial()
Poplar.ForceDensity = 29.376 * LBF/FT**3 #Measured 1/15/09 Bob Klenke, Robert Knapke

# for "Finish Birch aircraft plywood, source <http://www.rcgroups.com/forums/showthread.php?t=592951>
AircraftPly = ACMaterial()
AircraftPly.ForceDensity = 40.7* LBF/FT**3 #Measured by Tommy Marks, Amanda McGee 9/2012
#previously 40.6 * LBF/FT**3

# Aluminum Balsa Combo Spar, Prototype 1 11/20/2016
AluminumBalsa = ACMaterial()
AluminumBalsa.ForceDensity = (3.2144/0.194444) * LBF/FT**3 # lbsf spar weight / volume of spar

#
# Metals
#

Steel = ACMaterial()
Steel.ForceDensity = 0.284*LBF/IN**3
Steel.E = 29000000*PSI

#print "Measured Mass Density: ", AsUnit(47*GRAM/(2*3*1/16*IN**3), 'g/in**3')
#print "Steel Mass Density: ", AsUnit(Steel.Density, 'g/in**3')

Aluminum = ACMaterial()
Aluminum.ForceDensity = 0.098*LBF/IN**3
Aluminum.E = 10300000*PSI

AluminumTube = ACMaterial()
AluminumTube.LinearForceDensity = 0.008*LBF/IN
AluminumTube.E = 10300000*PSI

#
# Foams
#

PinkFoam = ACMaterial()
PinkFoam.ForceDensity = 1.56* LBF/FT**3 #Measured by Tommy Marks, Amanda McGee 9/2012
#previously 1.4 * LBF/FT**3

BlueFoam = ACMaterial()
BlueFoam.ForceDensity = 1.3* LBF/FT**3 #Measured by Tommy Marks, Amanda McGee 9/2012
#previously 2.13 * LBF/FT**3

SpyderFoam = ACMaterial()
SpyderFoam.ForceDensity = 2.3 * LBF/FT**3

#
# Coverings
#

Monokote = ACMaterial()
Monokote.AreaDensity = 0.05*OZM/(6*IN*5*IN) #Measured by Tommy Marks, Amanda McGee 9/2012
#Previously 1.5*OZM/(12*IN * 30*IN)

Ultracote = ACMaterial()
Ultracote.AreaForceDensity = 0.05*OZF/(5.75*IN*5.25*IN) #Measured by Tommy Marks, Amanda McGee 9/2012
#Previously 0.13*OZF/(FT**2)

#
# Plastics
#
