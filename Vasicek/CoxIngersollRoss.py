############################################
# CoxIngersolRoss.py
############################################
# Description:
# * Pricing using the CoxIngersolRoss model.

import math as m
from scipy.stats import chi2, ncx2, norm
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D
from Pricing import InterestRatePricing
from Vasicek import VasicekParam

class CoxIngersollRossPricing(InterestRatePricing):
    """
    * Singleton class containing pricing methods 
    obeying the mean-varying and bounded model:
    dr = alpha (mu - r) + sigma * sqrt(r) * dW.
    """
    #####################
    # Constructors:
    #####################
    def __init__(self, params):
        if not isinstance(params, CIRParams):
            raise Exception('params must be a CIRParams object.')
        self.Params = params
        self.ModelName = 'CoxIngersollRoss'

    #####################
    # Methods:
    #####################
    def ZeroCouponBondFV(self, today, bondMaturity):
        """ (Equation 7.57 & 7.63)
        * Return price of zero coupon bond obeying the CoxIngersollRoss model.
        Inputs:
        * today: Years from today (numeric, non-negative).
        * bondMaturity: Years from today until bond matures (numeric, non-negative).
        """
        errMsgs = []
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric and non-negative.')
        if not InterestRatePricing.ValidTenor(bondMaturity):
            errMsgs.append('bondMaturity must be numeric and non-negative.')
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))

        a = self.Params.A(today, bondMaturity)
        b = self.Params.B(today, bondMaturity)
        r = self.Params.InstantaneousRate

        return a * m.exp(-r * b)

    def ZCBFuturesFV(self, today, futureExp, bondMaturity):
        """ (Equation 7.63)
        * Return fair strike of futures contract.
        Inputs:
        * today: Years from present (numeric, non-negative).
        * futureExp: Years until futures expires from present (numeric, non-negative).
        * bondMaturity: Years until bond matures from present (numeric, non-negative).
        """
        errMsgs = []
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric, non-negative.')
        if not InterestRatePricing.ValidTenor(futureExp):
            errMsgs.append('futureExp must be numeric, non-negative.')
        if not InterestRatePricing.ValidTenor(bondMaturity):
            errMsgs.append('bondMaturity must be numeric, non-negative.')
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))

        params = self.Params
        t = today    
        T = futureExp
        s = bondMaturity
        r = params.InstantaneousRate
        c = params.C(t, T, s)
        d = params.D(t, T, s)

        return c * m.exp(-r * d)

    def ZCBOptionFV(self, today, optionExp, bondMaturity, strike):
        """ (Equation 7.68)
        * Calculate fair value of option on zero coupon bond.
        Inputs:
        * today: Years from present to price option (numeric, non-negative).
        * optionExp: Option expiry, years from present (numeric, non-negative).
        * bondMaturity: Years from present until bond matures (numeric, non-negative)/
        * strike: Strike price on option (numeric, positive).
        """
        errMsgs = []
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric, non-negative.')
        if not InterestRatePricing.ValidTenor(optionExp):
            errMsgs.append('optionExp must be numeric, non-negative.')
        if not InterestRatePricing.ValidTenor(bondMaturity):
            errMsgs.append('bondMaturity must be numeric, non-negative.')
        if not InterestRatePricing.IsNumeric(strike):
            errMsgs.append('strike must be numeric.')
        elif strike <= 0:
            errMsgs.append('strike must be positive.')

        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))

        params = self.Params
        t = today
        T = optionExp
        s = bondMaturity
        mu = self.Params.Mu
        alpha = self.Params.Alpha
        a = self.Params.A(T, s)
        b = self.Params.B(T, s)
        r_inst = self.Params.InstantaneousRate 
        lam = self.Params.Lambda
        gam = self.Params.Gamma
        sig = self.Params.Sigma
        zero_t_s = self.ZeroCouponBond(t, s)
        zero_t_T = self.ZeroCouponBond(t, T)
        phi = 2 * gam / (sig * sig)
        phi /= m.exp(gam * (T - t)) - 1
        psi = (alpha + lam + gam) / (sig * sig)
        r = m.log(a / strike) / b
        d_1 = 2 * (phi + psi + b) * r
        d_2 = 2 * (phi + psi) * r
        df = 4 * alpha * mu / (sig * sig)
        df_nc = 2 * phi * phi * r_inst * m.exp(gam  * (s - t)) / (phi + psi)

        fv = zero_t_s * ncx2.ppf(d_1, df, df_nc)
        fv -= zero_t_T * strike * ncx2.ppf(d_2, df, df_nc)

        return fv

    def GenerateZeroCurve(self, tStart, tEnd, tStep):
        """
        * Generate zero coupon yield curve.
        Inputs:
        * tStart: Start year (numeric, non-negative).
        * tEnd: End year (numeric, non-negative).
        * tStep: Fraction of year step (numeric, positive).
        """
        return InterestRatePricing._InterestRatePricing__GenZeroCurve(self, tStart, tEnd, tStep)

    def CouponBondFV(self, coupon, cStart, cEnd, bMaturity, freq):
        """ (Equation 7.73)
        * Calculate present value of coupon bond using discount factor curve
        adhering to CIR model.
        Inputs:
        * coupon: Fixed % of face value paid at each coupon date (numeric, non-negative).
        * cStart: Coupon start date, in years from present (numeric, non-negative).
        * cEnd: Final coupon payment date, in years from present (numeric, non-negative, >= cStart).
        * bMaturity: Bond maturity date, when face value is returned, years from present (numeric, non-negative).
        * freq: Coupon payment frequency per year, in fraction of year (numeric, positive, <= 1).
        """
        errMsgs = []
        if not InterestRatePricing.IsNumeric(coupon):
            errMsgs.append('coupon must be numeric.')
        elif coupon < 0:
            errMsgs.append('coupon must be non-negative.')
        if not InterestRatePricing.ValidTenor(cStart):
            errMsgs.append('cStart must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(cEnd):
            errMsgs.append('cEnd must be numeric non-negative.')
        elif cEnd < cStart:
            errMsgs.append('cStart must be >= cEnd.')
        if not InterestRatePricing.ValidTenor(bMaturity):
            errMsgs.append('bMaturity must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(freq):
            errMsgs.append('freq must be numeric non-negative.')
        elif freq > 1:
            errMsgs.append('freq must be <= 1.')

        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))
        # Generate discount factor curve using model:
        curve = self.GenerateZeroCurve(cStart, bondMaturity, freq)
        # Calculate present value of bond according to set schedule:
        pv = 0
        while cStart <= cEnd:
            pv += curve[cStart]
            cStart += freq
        pv *= coupon
        pv += curve[bondMaturity]

        return pv

    def CouponBondFutureFV(self, today, coupon, cStart, cEnd, freq, futureMaturity, bondMaturity):
        """ (Equation 7.73)
        * Calculate fair strike of coupon bond future contract using CIR model.
        Inputs:
        * today: Years from present to value futures contracts (numeric, non-negative).
        * coupon: Fixed % of face value paid at each coupon date (numeric, non-negative).
        * cStart: Coupon start date, in years from present (numeric, non-negative).
        * cEnd: Final coupon payment date, in years from present (numeric, non-negative, >= cStart).
        * freq: Coupon payment frequency per year, in fraction of year (numeric, positive, <= 1).
        * futureMaturity: Maturity date of futures contract, years from present (numeric, non-negative).
        * bMaturity: Bond maturity date, when face value is returned, years from present (numeric, non-negative).
        """
        errMsgs = []
        if not InterestRatePricing.IsNumeric(coupon):
            errMsgs.append('coupon must be numeric.')
        elif coupon < 0:
            errMsgs.append('coupon must be non-negative.')
        if not InterestRatePricing.ValidTenor(cStart):
            errMsgs.append('cStart must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(cEnd):
            errMsgs.append('cEnd must be numeric non-negative.')
        elif cEnd < cStart:
            errMsgs.append('cStart must be >= cEnd.')
        if not InterestRatePricing.ValidTenor(bondMaturity):
            errMsgs.append('bondMaturity must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(futureMaturity):
            errMsgs.append('futureMaturity must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(freq):
            errMsgs.append('freq must be numeric non-negative.')
        elif freq > 1:
            errMsgs.append('freq must be <= 1.')
        
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))
        # Calculate present value of bond according to set schedule:
        val = 0
        while cStart <= cEnd:
            val += self.ZCBFuturesFV(today, futureMaturity, cStart)
            cStart += freq

        val *= coupon
        val += self.ZCBFuturesFV(today, futureMaturity, bondMaturity)
        
        return val

    def CouponBondOptionFV(self, today, coupon, freq, optionExpiry, bondMaturity):
        """
        * 
        """
        if not InterestRatePricing.IsNumeric(coupon):
            errMsgs.append('coupon must be numeric.')
        elif coupon < 0:
            errMsgs.append('coupon must be non-negative.')
        if not InterestRatePricing.ValidTenor(cStart):
            errMsgs.append('cStart must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(cEnd):
            errMsgs.append('cEnd must be numeric non-negative.')
        elif cEnd < cStart:
            errMsgs.append('cStart must be >= cEnd.')
        if not InterestRatePricing.ValidTenor(bondMaturity):
            errMsgs.append('bondMaturity must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(futureMaturity):
            errMsgs.append('futureMaturity must be numeric non-negative.')
        if not InterestRatePricing.ValidTenor(freq):
            errMsgs.append('freq must be numeric non-negative.')
        elif freq > 1:
            errMsgs.append('freq must be <= 1.')

        pass
class CIRParams(VasicekParam):
    """
    * Object stores parameters used in the pricing methods in
    CoxIngersolRoss object.
    """
    __args = VasicekParam.ArgDict()
    def __init__(self, params):
        """
        * Initiate new parameter object for CIRPricing object.
        Ensures that all parameter values are valid.
        Possible parameters are [alpha, lambda, mu, r, T, t].
        Call VasicekParams.ArgDict() to get copy of default kwargs dictionary.
        """
        self._VasicekParam__ValidateAndSet(params)
    ####################
    # Interface Methods:
    ####################
    def c(self, t_1, t_2):
        """
        * Return 'c' value useful in pricing.
        """
        a = self.Alpha
        sig = self.Sigma
        exp = m.exp(-a * (t_2 - t_1))
        c = 2 * a / (sig * sig)
        c /= (1 - exp)
        return c
    def c_RN(self, t_1, t_2):
        """
        * Return risk-neutral version of 'c' value useful in pricing bond futures.
        """
        a = self.Alpha_RN
        sig = self.Sigma
        exp = m.exp(-a * (t_2 - t_1))
        c = 2 * a / (sig * sig)
        c /= (1 - exp)
        return c
    def epsilon(self, t_1, t_2):
        """
        * Return 
        """
        pass
    def D(self, t_1, t_2, t_3):
        """
        * Return D value useful in pricing bond futures.
        """
        a = self.Alpha_RN 
        b = self.B(t_2, t_3)
        c = self.c_RN(t_1, t_2)
        exp = m.exp(- a * (t_2 - t_1))
        out = b * c * exp
        out /= c + b

        return out
    def C(self, t_1, t_2, t_3):
        """
        * Return C value useful in pricing bond futures.
        """
        a = self.A(t_2, t_3)
        b = self.B(t_2, t_3)
        c = self.c_RN(t_1, t_2)
        sig = self.Sigma
        mu = self.Mu
        al = self.Alpha
        out = a * (c / (c + b)) ** (2 * al * mu / (sig * sig))
        return out 

    def B(self, t_1, t_2):
        """
        * Return B useful in ZCB price calculation.
        """
        return 2 * (m.exp(self.Gamma * (t_2 - t_1)) - 1) / self.__Denom(t_1, t_2)
    def A(self, t_1, t_2):
        """
        * Return A useful in ZCB price calculation.
        """
        num = 2 * self.Gamma * m.exp((self.Alpha + self.Lambda + self.Gamma) * (t_2 - t_1) / 2)
        return (num / self.__Denom(t_1, t_2)) ** (2 * self.Alpha * self.Mu / (self.Sigma * self.Sigma))
    def __Denom(self, t_1, t_2):
        """
        * Return denominator used in several formulas.
        """
        return (self.Alpha + self.Lambda + self.Gamma) * (m.exp(self.Gamma * (t_2 - t_1)) - 1) + 2 * self.Gamma
    ####################
    # Properties:
    ####################
    @property
    def Alpha_RN(self):
        return self.Alpha + self.Lambda
    @property
    def Mu_RN(self):
        out = self.Alpha * self.Mu
        return out / (self.Alpha + self.Mu)
    @property
    def Gamma(self):
        return m.sqrt(self.Alpha_RN ** 2 + 2 * self.Sigma * self.Sigma)
    ####################
    # Static Methods:
    ####################
    @staticmethod
    def __StaticConstructor(args):
        """
        * Fill __args with the base class' args.
        """
        args = VasicekParam._VasicekParam__args.copy()

    @staticmethod
    def ArgDict():
        """
        * Return copy of default params dictionary used in 
        static pricing methods.
        """
        return CIRParams.__args.copy()