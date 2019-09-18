"""
This is the class for analyzing propellers
Used for Propulsion analysis
"""

from __future__ import division # let 5/2 = 2.5 rather than 2
from ACBase import Velocity, Length, Angle, Unitless, PerAngle, Force, Power, Torque
from scalar.units import RPM, IN, ARCDEG, RAD, FT, SEC, LBF, HP, OZF
from scalar.units import AsUnit
from ACMass import ACMassCyl
import AeroUtil as AUtil
import cmath as math
import numpy as npy
import pylab as pyl

################################################################################
class ACPropellerError(Exception):
    def __init__(self, message):
        self.msg = message
    def __str__(self):
        return self.msg

################################################################################
class ACPropeller(ACMassCyl):
    """
    Propeller class
    
    Von Mises'
    Representative Blade Propeller Model
    A Simple Parametric Propeller Model

    Introduction
    Our Problem is that we need to investigate a family of props in order to determine the best one for our airplane. This
    means we will need a map for each different prop.  We'll be lucky if we can get exactly what we need, so the
    following can be used to generate a representative family of maps using one propeller data set as a starting point and
    changing beta and solidity.  Von Mises' Theory of flight, representative blade theory was used as a reference,
    pgs 302-309. Two prop data sets were used to correlate and check the behavior of this model. They are listed below.
    You should correlate this model to a real propeller that is close to what you want. For instance if you need a small
    model aircraft propeller, you should find some data on that type and recorrelate this model to it.  You get the idea.
    
    
    Attributes:
        D            - Propeller diameter
        Thickness    - Propeller thickness at the hub (for drawing purposes)
        RD           - The location on the profile chord where the PitchAngle is defined (default 3/8) von Mises 306
        PitchAngle   - Propeller pitch angle (Also known as "Beta")
        Pitch        - Pitch at 75% along the blade
        dAlpha       - Difference between measured alpha and zero lift alpha
        AlphaZeroCL  - Angle of attack with zero lift (default 0*deg)
        AlphaStall   - Stall angle of attack
        AlphaDead    - Post stall angle of attack for zero lift (default AlphaStall + 4 deg)
        CLSlope      - 2D airfoil lift slope (default 0.068/deg)
        CDCurve      - 2D airfoil drag curvature (default 2.2)
        CDp          - Parasitic drag (default 0.01)
        Solidity     - Proportional to the blade disk area, similar to the activity factor (AreaBlades/(2*D**2))
        MaxTipSpeed  - Maximum desired tip speed (default 700 ft/s)

        --Units
        ThrustUnit     - Unit for thrust (default lbr)
        ThrustUnitName - Name of unit for thrust (default lbf)
        TorqueUnit     - Unit for thrust (default in*ozf)
        TorqueUnitName - Name of unit for thrust (default in*ozf)
        PowerUnit      - Unit for power (default hp)
        PowerUnitName  - Name of unit for power (default hp)

        ThrustData   - An array of static thrust data corrected for sea level and standard day
                    #            RPM,        Thrust
                    ThrustData = [(7000  *RPM, 9.1 *LBF),
                                  (8500  *RPM, 9.4 *LBF),
                                  (9500  *RPM, 9.8 *LBF)]

        TorqueData   - An array of torque data corrected for sea level and standard day
                    #              RPM,        torque
                    TorqueData = [(7000  *RPM, 107.35 *IN*OZF),
                                  (8500  *RPM, 104.24 *IN*OZF),
                                  (9500  *RPM, 101.13 *IN*OZF)]
        
        DynThrustData - An array of dynamic thrust data corrected for sea level and standard day
                    #            RPM,        Thrust,     Velocity
                    ThrustData = [(7000  *RPM, 9.1 *LBF, 20*FT/SEC),
                                  (8500  *RPM, 8.5 *LBF, 25*FT/SEC),
                                  (9500  *RPM, 8.1 *LBF, 30*FT/SEC)]
    """

#===============================================================================    
    def __init__(self):
        super(ACPropeller,self).__init__()
        
        UnitList = self.UnitList
        
        self.name = 'Propeller'
        
        self.NoneList['D']           = None          ; UnitList['D']           = Length
        self.__dict__['Thickness']   = 0.25*IN       ; UnitList['Thickness']   = Length
        self.NoneList['PitchAngle']  = None          ; UnitList['PitchAngle']  = Angle
        self.NoneList['Pitch']       = None          ; UnitList['Pitch']       = Length
        self.__dict__['dAlpha']      = 0*ARCDEG      ; UnitList['dAlpha']      = Angle
        self.__dict__['AlphaZeroCL'] = 0*ARCDEG      ; UnitList['AlphaZeroCL'] = Angle
        self.__dict__['AlphaStall']  = 0*ARCDEG      ; UnitList['AlphaStall']  = Angle
        self.NoneList['AlphaDead']   = None          ; UnitList['AlphaDead']   = Angle
        self.__dict__['CLSlope']     = 0.068/ARCDEG  ; UnitList['CLSlope']     = PerAngle
        self.__dict__['CDCurve']     = 2.2           ; UnitList['CDCurve']     = Unitless
        self.__dict__['CDp']         = 0.01          ; UnitList['CDp']         = Unitless
        self.__dict__['Solidity']    = 0             ; UnitList['Solidity']    = Unitless
        self.__dict__['RD']          = 3/8           ; UnitList['RD']          = Unitless
        self.__dict__['MaxTipSpeed'] = 700 * FT/SEC  ; UnitList['MaxTipSpeed'] = Velocity
        self.__dict__['ThrustData']  = None
        self.__dict__['TorqueData']  = None
        self.__dict__['DynThrustData']  = None

        self.__dict__['ThrustUnit'] = LBF            ; UnitList['ThrustUnit'] = Force
        self.__dict__['ThrustUnitName'] = 'lbf'

        self.__dict__['TorqueUnit'] = IN*OZF         ; UnitList['TorqueUnit'] = Torque
        self.__dict__['TorqueUnitName'] = 'in*ozf'

        self.__dict__['PowerUnit'] = HP              ; UnitList['PowerUnit'] = Power
        self.__dict__['PowerUnitName'] = 'hp'          

#===============================================================================
    def Refresh(self):
        """
        Virtual function used to recompute attributes

        All derived instances of this methods must set
        self.param.refreshing = True at the beginning of the method and
        self.param.refreshing = False when completed.

        Otherwise an infinite loop will be created!!!
        """
        self.param.refreshing = True
        
        self._CheckConsistent()
        
        self.Axis  = [-1,0,0]
        self.LenDi = [self.Thickness, self.D]
        
        self.param.refreshing = False
        
        #
        # The now that the diameter is updated, call the refresh for ACMassCyl
        #
        super(ACPropeller, self).Refresh()

#===============================================================================
    def Length(self):
        """
        Returns the length of the propeller and nose cone
        """
        return self.Thickness + self.NoseCone.LenDi[0]
        
#===============================================================================
    def CL(self,Alpha):
        """
        Computes the propeller lift coefficient
        
        Inputs:
            Alpha - Angle of attack
        """
        if self.dirty: self.Refresh()
        
        cla    = self.CLSlope
        astall = self.AlphaStall
        a0     = self.AlphaZeroCL
        aDead  = self.AlphaDead
        
        def cl(ap):
            ap = ap - a0
            return cla*ap - 0.11*(ap/RAD)**3  #Approximate Cl rolloff with alpha

        Alpha = self.MakeList(Alpha)
        
        CL = []
        for a in Alpha:
            if a < -aDead:
                CL.append( cl(-astall)/2 )
            elif a < -astall and a >= -aDead:
                CL.append(cl(-astall))
            elif a > astall and a <= aDead:
                CL.append(cl(astall))
            elif a > aDead:
                CL.append( cl(astall)/2 )
            else:
                CL.append(cl(a))
                
        return CL[0] if len(CL) == 1 else npy.array(CL)
       
#===============================================================================
    def CD(self,Alpha):
        """
        Computes the propeller drag coefficient
        
        Inputs:
            Alpha - Angle of attack
        """
        if self.dirty: self.Refresh()

        aDead  = self.AlphaDead
        
        cdp = self.CDp
        cdc = self.CDCurve
        
        def cd(ap):
            return cdp + cdc*(ap / RAD)**2
        
        Alpha = self.MakeList(Alpha)
        
        CD = []
        for a in Alpha:
            CD.append( cd(a) )
                
        return CD[0] if len(CD) == 1 else npy.array(CD)
        
#===============================================================================        
    def CT(self,Alpha):
        """
        Computes the propeller thrust coefficient
        
        Inputs:
            Alpha  - angle of attack of the propeller blade
        """
        if self.dirty: self.Refresh()

        betap = self.PitchAngle + self.dAlpha
        lam   = self.RD*2*math.pi
        CL    = self.CL
        CD    = self.CD
        mu    = self.Solidity
        
        def ct(a):
            gam   = (betap - a) / RAD
            sing  = math.sin(gam).real
            cosg  = math.cos(gam).real
            return mu*lam**2*(CL(a)*cosg-CD(a)*sing)/cosg**2

        Alpha = self.MakeList(Alpha)
        aDead = self.AlphaDead

        CT = []
        for a in Alpha:
            if a < -aDead:
                CT.append( ct( -aDead ) )
            elif a > aDead:
                CT.append( ct( aDead ) )
            else: 
                CT.append( ct(a) )
        
        return CT[0] if len(CT) == 1 else npy.array(CT)
    
#===============================================================================        
    def CP(self,Alpha):
        """
        Computes the propeller power coefficient
        
        Inputs:
            Alpha  - angle of attack of the propeller blade
        """
        if self.dirty: self.Refresh()
        
        betap = self.PitchAngle + self.dAlpha
        lam   = self.RD*2*math.pi
        CL    = self.CL
        CD    = self.CD
        mu    = self.Solidity
        
        def cp(a):
            gam   = (betap - a) / RAD
            sing  = math.sin(gam).real
            cosg  = math.cos(gam).real
            return mu*lam**3*(CL(a)*sing+CD(a)*cosg)/cosg**2
        
        Alpha = self.MakeList(Alpha)
        aDead = self.AlphaDead
        
        CP = []
        for a in Alpha:
            if a < -aDead:
                CP.append( cp( -aDead ) )
            elif a > aDead:
                CP.append( cp( aDead ) )
            else: 
                CP.append( cp(a) )
            
        return CP[0] if len(CP) == 1 else npy.array(CP)
        
#===============================================================================
    def J(self,Alpha=None,N=None,V=None):
        """
        Computes the Advance ratio (J)
            - if N and V are given, Alpha is also calculated
            - do not enter both N and V as lists (only one can be a list)
        
        Inputs:
            Alpha  - angle of attack of the propeller blade
                OR
            N      - RPM of the propeller blade
            v      - forward velocity of the aircraft 
        """
        if self.dirty: self.Refresh()
        
        if Alpha is not None and N is None and V is None:
            # If the user specifies only the blade Alpha
            
            Alpha = self.MakeList(Alpha)
            
            J = []
            for a in Alpha:
                J.append(self._Jalpha(a))
                
            if len(J) == 1:
                return J[0]
            else:
                return npy.array(J)
        
        elif Alpha is None and N is not None and V is not None:
            # If the user specifies N and V
            try:
                Nlen = len(N)
                
            except:
                Nlen = 1
                
            try:
                Vlen = len(V)
                
            except:
                Vlen = 1
                
            if Nlen > 1 and Vlen > 1:
                errormsg = 'Both N and V are lists.  Only N or V can be a list'
                raise ACPropellerError(errormsg)
            
            N = self.MakeList(N)
            V = self.MakeList(V)
            
            J= []; Alpha = []     
            for n in N:
                for v in V:
                    JVal, AlphaVal = self._Jnv(n,v)
                    J.append(JVal)
                    Alpha.append(AlphaVal)
                        
            return (J[0], Alpha[0]) if len(J) == 1 else (npy.array(J), npy.array(Alpha))
            
#===============================================================================
    def _Jalpha(self,a):
        """
        Computes the Advance ratio
        
        Inputs:
            Alpha  - angle of attack of the propeller blade 
        """

        betap = self.PitchAngle + self.dAlpha
        lam   = self.RD*2*math.pi
        
        gam   = (betap - a) / RAD
        return lam*math.tan(gam).real
        
#===============================================================================
    def _Jnv(self,N,V):
        """
        Computes the Advance ratio and propeller angle of attack
        von Mises pg 303, 306
        
        Inputs:
            N      - RPM of the propeller blade
            v      - forward velocity of the aircraft 
        """
        
        d     = self.D
        betap = self.PitchAngle + self.dAlpha
        rd    = self.RD
        
        J     = V/(d*N/(2.*math.pi))
        Alpha = betap - (math.atan( J/(2*math.pi*rd) )).real*RAD
        
        return J, Alpha

#===============================================================================
    def Eta(self,Alpha):
        """
        Computes the propeller efficiency given
        
        Inputs:
            Alpha  - angle of attack of the propeller blade
        """
        if self.dirty: self.Refresh()
        
        J  = self.J
        CT = self.CT
        CP = self.CP
        
        Eta = []
        for a in Alpha:
            if CP(a) > 0 and CT(a) > 0:
                Eta.append(CT(a)*J(Alpha=a)/CP(a))
            elif CP(a) < 0 and CT(a) < 0:
                Eta.append(-CP(a)/(CT(a)*J(Alpha=a)))
            else:
                Eta.append(0)
                
        if len(Eta) == 0:
            return Eta[0]
        else:        
            return npy.array(Eta)
            
#===============================================================================
    def T(self,N,V,h):
        """
        Computes the thrust produced by the propeller
            - do not enter both N and V as lists
                
        Inputs:
            N      - RPM of the propeller blade
            v      - forward velocity of the aircraft 
            h      - altitude
        """
        if self.dirty: self.Refresh()
        
        CT  = self.CT
        rho = AUtil.DensAlt(h)
        d   = self.D
        
        try:
            Nlen = len(N)
            
        except:
            Nlen = 1
            
        try:
            Vlen = len(V)
            
        except:
            Vlen = 1
            
        if Nlen > 1 and Vlen > 1:
            errormsg = 'Both N and V are lists.  Only N or V can be a list'
            raise ACPropellerError(errormsg)
        
        N = self.MakeList(N)
        V = self.MakeList(V)
        
        T=[]    
        for n in N:
            for v in V:
                J, Alpha = self.J(N=n,V=v)
                if( n < 10*RPM):
                    T.append(0) #Just set it to zero if we are at this low of RPM, the calculation failes for RPM = 0
                else:
                    T.append((CT(Alpha) * rho * (n/(2*math.pi))**2 * d**4) / LBF)

        return T[0]*LBF if len(T) == 1 else npy.array(T)*LBF
        
#===============================================================================
    def P(self,N,V,h):
        """
        Computes the power absorbed by the propeller
            - do not enter both N and V as lists (only one can be a list)
            
        Inputs:
            N      - RPM of the propeller blade
            v      - forward velocity of the aircraft 
            h      - altitude
        """
        if self.dirty: self.Refresh()
        
        CP  = self.CP
        rho = AUtil.DensAlt(h)
        d   = self.D
        
        try:
            Nlen = len(N)
            
        except:
            Nlen = 1
            
        try:
            Vlen = len(V)
            
        except:
            Vlen = 1
            
        if Nlen > 1 and Vlen > 1:
            errormsg = 'Both N and V are lists.  Only N or V can be a list'
            raise ACPropellerError(errormsg)
        
        N = self.MakeList(N)
        V = self.MakeList(V)
        
        P=[]    
        for n in N:
            for v in V:
                J, Alpha = self.J(N=n,V=v)
                if( n < 10*RPM):
                    P.append(0) #Just set it to zero if we are at this low of RPM, the calculation fails for RPM = 0
                else:
                    P.append((CP(Alpha) * rho * (n/(2*math.pi))**3 * d**5) / HP)

        return P[0]*HP if len(P) == 1 else npy.array(P)*HP

#===============================================================================
    def MaxRPM(self, ):
        """
        Calculates the maximum RPM for a maximum desired tip speed
            
        Inputs:
            MaxTipVel - Maximum tip velocity desired
        """
        if self.dirty: self.Refresh()
        
        if self.MaxTipSpeed is not None:
            return self.MaxTipSpeed / (self.D * math.pi) * 2. * math.pi #Times 2 pi to get it in RPM
    
#===============================================================================
    def __getattr__(self, key):
        """
        Attempts to calculate any variables set to None
        """
        ans = super(ACPropeller, self).__getattr__(key)
        if ans is not None:
            return ans

        # Try to calculate the value instead
        if key == 'Pitch':
            return self._CalcPitch()
        elif key == 'PitchAngle':
            return self._CalcPitchAngle()
        elif key == 'D':
            return self._CalcD()
        elif key == 'AlphaDead':
            return self._CalcAlphaDead()

        #
        # If no calculation exist return None
        #
        return None

#===============================================================================
    def _CheckConsistent(self):
        """
        Makes sure that a consistent set of parameters has been specified.
        """
        self._CheckEquation(['Pitch','PitchAngle','D'], Need = 2)
        
#===============================================================================
    def _CalcPitch(self):
        """
        Calculates the pitch at 75% along the blade
        """
        d    = self.D
        beta = self.PitchAngle / RAD
        
        return 2*self.RD*math.pi*d*math.tan(beta).real

#===============================================================================
    def _CalcPitchAngle(self):
        """
        Calculates the pitch angle
        """        
        d     = self.D
        pitch = self.Pitch
        
        return ( math.atan(pitch/(2*self.RD*math.pi*d)).real )*RAD

#===============================================================================
    def _CalcD(self):
        """
        Calculates the blade diameter
        """        
        beta  = self.PitchAngle / RAD
        pitch = self.Pitch
        
        return pitch/(2*self.RD*math.pi*math.tan(beta).real)

#===============================================================================
    def _CalcAlphaDead(self):
        """
        Calculates the post stall angle of attack with zero CL
        """                
        return self.AlphaStall + 4*ARCDEG

#===============================================================================
    def CoefPlot(self,Alpha, fig = 1):
        """
        Plots:
            CL, CD and CL/CD vs. Alpha
            CT, CP and Eta vs. the Advanced Ratio (J)
            
        Inputs:
            Alpha  - angle of attack of the propeller blade
            fig    - Figure number
        """
        if self.dirty: self.Refresh()

        CL  = self.CL(Alpha)
        CD  = self.CD(Alpha)
        CT  = self.CT(Alpha)
        CP  = self.CP(Alpha)
        J   = self.J(Alpha=Alpha)
        Eta = self.Eta(Alpha)
        
        Alpha = Alpha / ARCDEG
        
        pyl.figure(fig)
        pyl.suptitle(self.name)
        
        xlabel_offset = -0.1
        ylabel_offset = -0.125
        
        pyl.subplot(231)
        pyl.plot(Alpha, CL)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_L$')
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
    
        pyl.subplot(232)
        pyl.plot(Alpha, CD)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_D$')
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
    
        pyl.subplot(233)
        pyl.plot(Alpha, CL/CD)
        pyl.xlabel(r'$\alpha^o$'); pyl.ylabel(r'$C_L/C_D$')
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
    
        pyl.subplot(234)
        pyl.plot(J, CT)
        pyl.xlabel(r'J'); pyl.ylabel(r'$C_T$')
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
    
        pyl.subplot(235)
        pyl.plot(J, CP)
        pyl.xlabel(r'J'); pyl.ylabel(r'$C_P$')
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
    
        pyl.subplot(236)
        pyl.plot(J, Eta)
        pyl.xlabel('J'); pyl.ylabel(r'$\eta$')
        pyl.gca().xaxis.set_label_coords(0.5,xlabel_offset)
        pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
        

#===============================================================================
    def _RPMTitle(self):
        '''
        Creates convineant figure title when looking at RPM data
        '''
        maxN = self.MaxRPM()
        title = self.name
        if maxN is not None:
            title += " ( Max RPM " + AsUnit(maxN, 'rpm', '%3.0f') + ')'
#            title += " at tip speed of " + AsUnit(self.MaxTipSpeed, 'ft/s') + " ) "
        pyl.suptitle(title, x = 0.1)

#===============================================================================
    def PTPlot(self,N,V,h,xAxisVar, fig = 1):
        """
        Plots: Power and Thrust vs. the aircraft velocity or RPM
            
        Inputs:
            N        - RPM of the propeller blade
            v        - forward velocity of the aircraft 
            h        - altitude
            xAxisVar - 'N','V'; defines which variable thrust and power are plotted against (A STRING!!!)
            fig      - Figure number
        """
        if self.dirty: self.Refresh()
        
        pyl.figure(fig)
        self._RPMTitle()
        maxN = self.MaxRPM()
        legend = []
        lines = []
        
        def plot(X, xlabel, T, P, Torque, color = ''):
            pyl.subplot(131)
            lines.extend([ X, T, color ])
            pyl.xlabel(xlabel); pyl.ylabel('Thrust (' + self.ThrustUnitName + ')')

            pyl.subplot(132)
            pyl.plot(X, Torque, color)
            pyl.xlabel(xlabel); pyl.ylabel('Torque (' + self.TorqueUnitName + ')')
        
            pyl.subplot(133)
            pyl.plot(X, P, color)
            pyl.xlabel(xlabel); pyl.ylabel('Power (' + self.PowerUnitName + ')')
        
        if xAxisVar == 'N':
            X=N / RPM
            xlabel = 'RPM'
            pyl.figure(fig)
            pyl.hold(True)
            
            for v in V:
                T = self.T(N,v,h)
                P = self.P(N,v,h)
                Torque = P/N
                                
                T = T / self.ThrustUnit
                P = P / self.PowerUnit
                Torque = Torque / self.TorqueUnit
                
                legend.append( 'V = ' + AsUnit( v, 'ft/s', '%3.0f' ) )
                
                plot(X, xlabel, T, P, Torque)
            
            legend.append( 'VTip = ' + AsUnit(self.MaxTipSpeed, 'ft/s') )
            X = maxN / RPM
            for v in V:
                T = self.T(maxN,v,h)
                P = self.P(maxN,v,h)
                Torque = P/maxN
                
                T = T / self.ThrustUnit
                P = P / self.PowerUnit
                Torque = Torque / self.TorqueUnit

                plot(X, xlabel, T, P, Torque, 'or')
                
        
        elif xAxisVar == 'V':
            X = V / (FT/SEC)
            xlabel = 'V (ft/s)'
            pyl.figure(fig)
            pyl.hold(True)
            
            for n in N:
                T = self.T(n,V,h)
                P = self.P(n,V,h)
                Torque = P / n
                
                T = T / self.ThrustUnit
                P = P / self.PowerUnit
                Torque = Torque / self.TorqueUnit
                
                legend.append( 'N = ' + AsUnit(n, 'rpm', '%3.0f') )
                
                plot(X, xlabel, T, P, Torque)
                
        else:
            errormsg = "Incorrect xAxis variable specified.  Specify 'N' or 'V' (as strings)"
            raise ACPropellerError(errormsg)

        pyl.subplot(131)
        lines = pyl.plot( *lines )
        pyl.figlegend(lines, legend, loc = 'upper right')
        #pyl.subplot(131); pyl.legend(legend, loc = 'best')
        #pyl.subplot(132); pyl.legend(legend, loc = 'best')
        #pyl.subplot(133); pyl.legend(legend, loc = 'best')

#===============================================================================
    def PlotTestData(self, fig = 1, color = ''):
        """
        Plots measured force data with the propeller model
            
        Inputs:
            fig - Figure number
        """
        if self.ThrustData is None and self.TorqueData is None and self.DynThrustData:
            raise ACPropellerError("No test data has been given")
        
        if self.ThrustData is not None and self.TorqueData is not None and self.DynThrustData is not None:
            subplot1 = 221; subplot2 = 222;  subplot3 = 223;  subplot4 = 224;
        elif self.ThrustData is not None and self.TorqueData is not None and self.DynThrustData is None:
            subplot1 = 131; subplot2 = 132;  subplot3 = 133;
        else:
            subplot1 = 111; subplot2 = 121; subplot3 = 122;
        
        ylabel_offset = -0.1
        
        pyl.figure(fig)
        self._RPMTitle()
        maxN = self.MaxRPM()
    
        legend = ["Propeller Model", "Test Data"]
        lines = []
        if maxN is not None:
            legend += ["VTip = " + AsUnit( self.MaxTipSpeed, 'ft/s') ]

        if self.ThrustData is not None:
            ThrustData = npy.array(self.ThrustData)
            
            RPMData   = self.ToUnumList(ThrustData[:,0], RPM) / RPM
            ForceData = self.ToUnumList(ThrustData[:,1], LBF) / self.ThrustUnit
    
            N = npy.linspace(0,max(RPMData)*1.05, 20)*RPM 
            
            Thrust = self.T(N, V=0*FT/SEC, h=0*FT) / self.ThrustUnit
            N = N / RPM
            
            if maxN is not None:
                maxT = self.T(maxN, V=0*FT/SEC, h=0*FT) / self.ThrustUnit
                MaxN = maxN / RPM
                        
            ax = pyl.subplot(subplot1)
            ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
            lines  = [ N, Thrust, 'b']
            lines += [RPMData, ForceData,'og' ]
            if maxN is not None:
                lines += [ MaxN, maxT, 'or']
            
            lines = pyl.plot( *lines )
            pyl.figlegend(lines, legend, loc = 'upper right', ncol = len(lines), frameon = False)
            pyl.ylim(ymin = 0, ymax = max( pyl.ylim()[1], max(Thrust)*1.1 ) )
            pyl.xlabel('RPM'); pyl.ylabel('Thrust (' + self.ThrustUnitName + ')')
            pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
            pyl.grid()

        if self.TorqueData is not None:
            TorqueData = npy.array(self.TorqueData)
            
            RPMData    = self.ToUnumList(TorqueData[:,0], RPM)
            TorqueData = self.ToUnumList(TorqueData[:,1], IN*OZF)
            
            PowerData = (RPMData*TorqueData) / self.PowerUnit
            RPMData = RPMData / RPM
            TorqueData = TorqueData / self.TorqueUnit
            
            N = npy.linspace(0,max(RPMData)*1.05, 20)*RPM 
            
            Power = self.P(N, V=0*FT/SEC, h=0*FT) 
            Torque = Power / N

            if maxN is not None:
                maxP = self.P(maxN, V=0*FT/SEC, h=0*FT) 
                maxT = maxP / maxN
            
            Power = Power / self.PowerUnit
            Torque = Torque / self.TorqueUnit
            N = N / RPM
            
            if maxN is not None:
                MaxN = maxN / RPM
                maxP = maxP / self.PowerUnit
                maxT = maxT / self.TorqueUnit
                
            ax = pyl.subplot(subplot2)
            ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
            pyl.plot(N, Torque, color)
            pyl.plot(RPMData, TorqueData,'og')
            if maxN is not None:
                pyl.plot(MaxN, maxT, 'or')
            #pyl.legend(legend, loc = 'best')
            pyl.xlabel('RPM'); pyl.ylabel('Torque (' + self.TorqueUnitName + ')')
            pyl.ylim(ymin = 0)
            pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
            pyl.grid()
            
            ax = pyl.subplot(subplot3)
            ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
            pyl.plot(N, Power, color)
            pyl.plot(RPMData, PowerData,'og')
            if maxN is not None:
                pyl.plot(MaxN, maxP, 'or')
            if self.ThrustData is None:
                pyl.legend(legend, loc = 'best')
            pyl.xlabel('RPM'); pyl.ylabel('Power (' + self.PowerUnitName + ')')
            pyl.ylim(ymin = 0)
            pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
            pyl.grid()
        
        if self.DynThrustData is not None:
            DynThrustData = npy.array(self.DynThrustData)
            
            RPMData      = self.ToUnumList(DynThrustData[:,0], RPM)
            ForceData    = self.ToUnumList(DynThrustData[:,1], LBF) / self.ThrustUnit
            VelocityData = self.ToUnumList(DynThrustData[:,2], FT/SEC) / (FT/SEC)
            
            ThrustModel = npy.zeros(len(VelocityData))   
            for i in xrange(len(DynThrustData)):
                N, F, V = DynThrustData[i,:]
                ThrustModel[i] = self.T(N, V=V, h=0*FT) / self.ThrustUnit

                                                 
            ax = pyl.subplot(subplot4)
            #ax.xaxis.get_major_formatter().set_powerlimits((-3,3))
            pyl.plot( VelocityData, ForceData, 'og' )
#            pyl.plot( VelocityData, ThrustModel, 'b' )             
            pyl.ylim(ymin = 0, ymax = max( pyl.ylim()[1], max(ForceData)*1.1 ) )
            pyl.xlim(xmin = 0)
            pyl.xlabel('Velocity (ft/s)'); pyl.ylabel('Thrust (' + self.ThrustUnitName + ')')
            pyl.gca().yaxis.set_label_coords(ylabel_offset,0.5)
            pyl.grid()
            
################################################################################
if __name__ == '__main__':
    Prop = ACPropeller()
    Prop.D          = 14.2*IN
    Prop.PitchAngle = 12*ARCDEG
#    Prop.Pitch      = 7.1117*IN 
    
    Prop.dAlpha     = 4*ARCDEG
    Prop.Solidity   = 0.01361
    Prop.RD         = 3/8
    Prop.AlphaStall = 14*ARCDEG
    N = 11400 * RPM

    # These are junk numbers for debugging purposes
    #                 RPM,        Thrust
    Prop.ThrustData = [(7000  *RPM, 9.1 *LBF),
                       (8500  *RPM, 9.4 *LBF),
                       (9500  *RPM, 9.8 *LBF)]
    
    Thurst = Prop.T(N,0*FT/SEC,0*FT)
    Thurst = Prop.T(N,50*FT/SEC,0*FT)
    
    print 'Pitch:          ', Prop.Pitch
    print 'PitchAngle:     ', Prop.PitchAngle
    print 'Diameter:       ', Prop.D
    print 'Power Absorbed: ', AsUnit(Prop.P(N, 0*FT/SEC, 0*FT), "hp")
    
    h=0*FT
    N=npy.linspace(10000,16000,4)*RPM #11400*RPM
    
    Alpha = npy.linspace(-25,25,31)*ARCDEG
    V     = npy.linspace(0,100,11)*FT/SEC
    
    Prop.CoefPlot(Alpha,fig = 1)
    Prop.PTPlot(N,V,h,'V', fig = 2)
    
    Prop.PlotTestData(fig=3)

    J, Alpha = Prop._Jnv(N=11400 * RPM, V=40*FT/SEC)
    print J, Alpha/ARCDEG
    
    pyl.show()


#####OLD TEST#####

#    CL  = Prop.CL(Alpha)
#    CD  = Prop.CD(Alpha)
#    CT  = Prop.CT(Alpha)
#    CP  = Prop.CP(Alpha)
#    J   = Prop.J(Alpha=Alpha)
#    Eta = Prop.Eta(Alpha)
#    T   = Prop.T(N,V,h)
#    P   = Prop.P(N,V,h)
#    J2, Alp2 = Prop.J(N=N,V=V)
#    rho = AUtil.DensAlt(h)
#    CP2 = Prop.CP(Alp2)
#    CP3 = Prop.CP(16*ARCDEG)
#    print 'Thrust   : ', T
#    print 'Power    : ', P
#    print 'Velocity : ', V
#    print 'J        : ', J
##    print 'J2       : ', J2
##    print 'CP2      : ', CP2
##    print 'Alpha2   : ', Alp2
##    print 'CP3      : ', CP3   
#    
#    Alpha = Alpha / (ARCDEG)
#    V     = V / (FT/SEC)
#    T     = T / (LBF)
#    P     = P / (HP)
#    
#    pyl.figure(1)
#    pyl.subplot(231)
#    pyl.plot(Alpha, CL)
#    pyl.xlabel(r'$\alpha[^o]$'); pyl.ylabel(r'$C_L$')
#    
#    pyl.subplot(232)
#    pyl.plot(Alpha, CD)
#    pyl.xlabel(r'$\alpha[^o]$'); pyl.ylabel(r'$C_D$')
#    
#    pyl.subplot(233)
#    pyl.plot(Alpha, CL/CD)
#    pyl.xlabel(r'$\alpha[^o]$'); pyl.ylabel(r'$C_L/C_D$')
#    
#    pyl.subplot(234)
#    pyl.plot(J, CT)
#    pyl.xlabel(r'J'); pyl.ylabel(r'$C_T$')
#    
#    pyl.subplot(235)
#    pyl.plot(J, CP)
#    pyl.xlabel(r'J'); pyl.ylabel(r'$C_P$')
#    
#    pyl.subplot(236)
#    pyl.plot(J, Eta)
#    pyl.xlabel('J'); pyl.ylabel(r'$\eta$')
#    
#    pyl.figure(2)
#    pyl.subplot(121)
#    pyl.plot(V, T)
#    pyl.xlabel('V'); pyl.ylabel('Thrust')
#    
#    pyl.subplot(122)
#    pyl.plot(V, P)
#    pyl.xlabel('V'); pyl.ylabel('Power')
    