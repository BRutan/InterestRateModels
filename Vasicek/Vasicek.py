############################################
# Vasicek.py
############################################
# Description:
# * Define VasicekPricing class (derived from InterestRatePricing abstract class) to
# price interest rate products and derivatives following Vasicek's mean reverting interest rate SDE
# dr = alpha (mu - r) + sigma * dW.

import math as m
from scipy.stats import norm
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D
from Pricing import InterestRatePricing, ParamsObj

__all__ = [ 'VasicekParam', 'VasicekPricing' ]

class VasicekPricing(InterestRatePricing):
    """
    * Class that implements pricing obeying the Vasicek interest rate model
    dr = alpha (mu - r) + sigma * dW.
    """
    ###############################
    # Constructors:
    ###############################
    def __init__(self, params):
        """
        * Initialize class with VasicekParam object.
        Inputs:
        * params: Expecting a VasicekParam object.
        """
        self.Params = params
        self.ModelName = 'Vasicek'
    ###############################
    # Properties:
    ###############################
    @property
    def Params(self):
        return InterestRatePricing.Params
    @Params.setter
    def Params(self, params):
        if not isinstance(params, VasicekParam):
            raise Exception('params must be a VasicekParam object.')
        InterestRatePricing.Params = params
    ###############################
    # Interface Methods:
    ###############################
    def GenerateZeroCurve(self, tStart, tEnd, tStep):
        """
        * Generate zero coupon yield curve.
        Inputs:
        * tStart: Start year (numeric, non-negative).
        * tEnd: End year (numeric, non-negative).
        * tStep: Fraction of year step (numeric, positive).
        """
        return InterestRatePricing._InterestRatePricing__GenZeroYieldCurve(self, tStart, tEnd, tStep)

    def GenerateDiscountFactors(self, tStart, tEnd, tStep):
        """
        * Generate discount factor curve.
        Inputs:
        * tStart: Start year (numeric, non-negative).
        * tEnd: End year (numeric, non-negative).
        * tStep: Fraction of year step (numeric, positive).
        """
        return InterestRatePricing._InterestRatePricing__GenDiscountFactorCurve(self, tStart, tEnd, tStep)

    def ZeroCouponBondFV(self, today, bondMaturity):
        """ (From Equation 7.30)
        * Calculate fair value of zero coupon bond obeying Vasicek's mean reversion model.
        dr = alpha (mu - r) + sigma * dW
        Inputs:
        * today: Years from present (numeric, non-negative).
        * bondMaturity: Years until bond matures from present (numeric, non-negative).
        """
        errMsgs = []
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric and non-negative.')
        if not InterestRatePricing.ValidTenor(bondMaturity):
            errMsgs.append('bondMaturity must be numeric and non-negative.')
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))

        params = self.Params
        r = params.InstantaneousRate
        sig = params.Sigma
        f = params.F(today, bondMaturity)
        g = params.G(today, bondMaturity)

        return m.exp(-r * f - g)

    def ZCBOptionFV(self, strike, today, T_option, bondMaturity):
        """ (From Equation 7.43)
        * Calculate price of option on zero coupon bond obeying Vasicek's mean reversion model.
        Inputs:
        * strike: Bond option strike price (numeric, non-negative).
        * today: Today's date from present in years (numeric, non-negative).
        * T_option: Year that option expires (numeric, positive).
        * bondMaturity: Maturity of bond in years (numeric, non-negative).
        We note in the formula that s corresponds to the expiration of the bond (from params object), T the expiration of the option.
        """
        errMsgs = []
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append("today must be numeric and non-negative.")
        if not InterestRatePricing.ValidTenor(T_option):
            errMsgs.append("T_option must be numeric and non-negative.")
        if not InterestRatePricing.IsNumeric(strike):
            errMsgs.append("strike must be numeric.")
        if len(errMsgs) > 0:
            raise Exception(''.join(errMsgs))

        params = self.Params
        t = today
        s = bondMaturity
        T = T_option
        zero_t_s = self.ZeroCouponBondFV(t, s)
        zero_t_T = self.ZeroCouponBondFV(t, T)
        f = params.F(T, s)
        sig = params.Sigma
        alpha = params.Alpha
        # Calculate v_p, m_p_1, m_p_2:
        v_p = 1 - m.exp(-2 * alpha * (T - t))
        v_p /= 2 * alpha
        v_p *= (sig * f) ** 2
        m_p_2 = m.log(zero_t_s / zero_t_T) - .5 * v_p
        m_p_1 = m_p_2 + v_p
        # Calculate option value:
        d_1 = (m_p_1 - m.log((m_p_1 - m.log(strike))) / m.sqrt(v_p))
        d_2 = (m_p_2 - m.log((m_p_1 - m.log(strike))) / m.sqrt(v_p))
        optVal = zero_t_s * norm.cdf(d_1) - zero_t_T * strike * norm.cdf(d_2)

        return optVal

    def ZCBFuturesFV(self, today, futureExp, bondMaturity):
        """
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
        r = params.InstantaneousRate
        t = today    
        s = bondMaturity
        T = futureExp
        # Calculate X(t, T, s) and Y(t, T, s):
        x = params.X(t, T, s)
        y = params.Y(t, T, s)

        return m.exp(-r * x - y)

    def ZCBFuturesOptionFV(self, strike, today, optionExp, bondExp, futureExp):
        """ (Equation 7.46)
        * Calculate price of option on future.
        Inputs:
        * strike: Strike on option (numeric, positive).
        * today: years from present to price derivative (numeric, non-negative).
        * optionExp: option expiry, years from present (numeric, non-negative).
        * bondExp: bond maturity, years from present (numeric, non-negative).
        * futureExp: future contract expiration, years from present (numeric, non-negative).
        """
        errMsgs = []
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric, non-negative.')
        if not InterestRatePricing.ValidTenor(optionExp):
            errMsgs.append('optionExp must be numeric, non-negative.')
        if not InterestRatePricing.ValidTenor(bondExp):
            errMsgs.append('bondExp must be numeric, non-negative.')
        if not InterestRatePricing.ValidTenor(futureExp):
            errMsgs.append('futureExp must be numeric, non-negative.')
        if not InterestRatePricing.IsNumeric(strike):
            errMsgs.append('strike must be numeric.')
        elif strike <= 0:
            errMsgs.append('strike must be positive.')
        
        params = self.Params
        t = today
        T = optionExp
        s = futureExp
        w = bondExp
        a = params.Alpha
        sig = params.Sigma
        
        x_T_s_w = params.X(T, s, w)
        zero_t_T = self.ZeroCouponBondFV(today, s)
        fut_t_s_w = self.ZCBFuturesFV(today, s, w)
        f_t_T = params.F(t, T)

        h_t = zero_t_T * fut_t_s_w * m.exp(sig * sig / 2 * f_t_T * f_t_T * x_T_s_w) 
        v_h =  x_T_s_w * x_T_s_w
        v_h *= sig * sig * (1 - m.exp(-2 * a * (T - t)))
        v_h /= 2 * a
        d_h = 1 / m.sqrt(v_h)
        d_h *= (m.log(h_t / (zero_t_T * strike)) + v_h / 2)

        price = h_t * norm.cdf(d_h) - zero_t_T * strike * norm.cdf(d_h - m.sqrt(v_h))

        return price


class VasicekParam(ParamsObj):
    """
    * Object serves as parameter to Vasicek pricing static functions.
    """
    __args = {"alpha" : .001, "lambda" : 0, "mu" : 0, "r" : 0, "sigma" : 0}
    def __init__(self, argDict):
        """
        * Initiate new parameter object for VasicekModel object.
        Ensures that all parameter values are valid.
        Possible parameters are [alpha, lambda, mu, r, T, t].
        Call VasicekParams.ArgDict() to get copy of default kwargs dictionary.
        """
        self.__ValidateAndSet(argDict)

    #################
    # Getters:
    #################
    @property
    def Alpha(self):
        return self.__alpha
    @property
    def InstantaneousRate(self):
        return self.__instRate
    @property
    def Lambda(self):
        return self.__lambda
    @property
    def Mu(self):
        return self.__mu
    @property
    def Sigma(self):
        return self.__sigma
    @property
    def Params(self):
        """
        * Return dictionary containing this object's parameters.
        """
        return {"alpha" : self.Alpha, "lambda" : self.Lambda, "mu" : self.Mu, "r" : self.InstantaneousRate, "sigma" : self.Sigma}
    #################
    # Setters:
    #################
    @Alpha.setter
    def Alpha(self, alpha):
        if not isinstance(alpha, int) and not isinstance(alpha, float):
            raise Exception("Alpha must be numeric.")
        elif alpha <= 0:
            raise Exception("Alpha must be positive.")
        self.__req['alpha'] = True
        self.__alpha = alpha
    @InstantaneousRate.setter
    def InstantaneousRate(self, rate):
        if not isinstance(rate, int) and not isinstance(rate, float):
            raise Exception("Instantaneous rate must be numeric.")
        elif rate < 0:
            raise Exception("Instantaneous rate must be positive.")
        self.__req['r'] = True
        self.__instRate = rate
    @Lambda.setter
    def Lambda(self, lambd):
        if not isinstance(lambd, int) and not isinstance(lambd, float):
            raise Exception("Lambda must be numeric.")
        self.__req['lambda'] = True
        self.__lambda = lambd
    @Mu.setter
    def Mu(self, mu):
        if not isinstance(mu, int) and not isinstance(mu, float):
            raise Exception("Mu must be numeric.")
        self.__req['mu'] = True
        self.__mu = mu
    @Sigma.setter
    def Sigma(self, sig):
        if not isinstance(sig, int) and not isinstance(sig, float):
            raise Exception("Sigma must be numeric.")
        elif sig < 0:
            raise Exception("Sigma must be non-negative.")
        self.__req['sigma'] = True
        self.__sigma = sig
    @Params.setter
    def Params(self, params):
        if not isinstance(params, dict):
            raise Exception("Params must be a dictionary.")
        for key in params.keys():
            self.__SetAttribute(key, params[key])
    #################
    # Interface Functions:
    #################
    def F(self, t_1, t_2):
        """ 
        * Return F(t, T) used in several pricing formulae.
        Inputs:
        * t_1: Expecting numeric with t_2 > t_1.
        * t_2: Expecting numeric with t_2 > t_1.
        """
        f = 1 - m.exp(-self.Alpha * (t_2 - t_1))
        f /= self.Alpha
        return f

    def G(self, t_1, t_2):
        """
        * Return G(t, T) used in several pricing formulae.
        """
        mu = self.Mu
        sig = self.Sigma
        alpha = self.Alpha
        lambd = self.Lambda
        T = t_2
        t = t_1
        f = self.F(t_1, t_2)
        g = (mu - sig * sig / (2 * alpha * alpha) - sig * lambd / alpha)
        g *= T - t - f
        g += sig * sig * f * f / (4 * alpha)
    
        return g

    def X(self, t_1, t_2, t_3):
        """
        * Calculate X(t, T, s) function used in many pricing formulas.
        """
        f_t_s = self.F(t_1, t_3)
        f_t_T = self.F(t_1, t_2)
        
        return f_t_s - f_t_T

    def Y(self, t_1, t_2, t_3):
        """
        * Calculate Y(t, T, s) function used in many pricing formulas.
        """
        a = self.Alpha
        mu = self.Mu
        lam = self.Lambda
        sig = self.Sigma
        f_T_s = self.F(t_2, t_3)
        x = self.X(t_1, t_2, t_3)
        y = mu - lam * sig / a - sig * sig / (2 * a * a)
        y *= t_3 - t_2 - x
        y -= sig * sig / (2 * a * a) * (x - a / 2 * x * x - f_T_s)
        
        return y
        
    #################
    # Helper Methods:
    #################
    def __SetAttribute(self, paramStr, val):
        """
        * Set attribute of object parameter to value.
        """
        setattr(self, VasicekParam.__AttrKeyToAttr(paramStr), val)
    #################
    # Static Methods:
    #################
    @staticmethod
    def __AttrKeyToAttr(paramStr):
        """
        * Return attribute associated with key in dictionary.
        """
        if paramStr == 'alpha':
            return 'Alpha'
        if paramStr == 'lambda':
            return 'Lambda'
        if paramStr == 'r':
            return 'InstantaneousRate'
        if paramStr == 'mu':
            return 'Mu'
        if paramStr == 'sigma':
            return 'Sigma'
    @staticmethod
    def ArgDict():
        """ 
        * Return default dictionary containing all possible 
        arguments to the constructor.
        """
        return VasicekParam.__args.copy()
    #################
    # Private Helpers:
    #################
    def __ValidateAndSet(self, argDict):
        """
        * Ensure all parameters are valid.
        """
        # Ensure that all arguments have been passed:
        self.__req = {"alpha" : False, "lambda" : False, "mu" : False, "r" : False, "sigma" : False}
        errMsgs = []
        for arg in argDict.keys():
            if arg != 'T':
                tempArg = str(arg).lower() 
                try:
                    self.__SetAttribute(arg, argDict[arg])
                except Exception as ex:
                    errMsgs.append(ex.message)

        # Set T last (ensure that t is not greater than T):
        if "T" in argDict.keys() and "t" in argDict.keys():
            try:
                self.T = argDict['T']
            except Exception as ex:
                errMsgs.append(ex.message)

        # Ensure that all required arguments were set:
        unsetArgs = []
        for arg in self.__req.keys():
            if self.__req[arg] == False:
                unsetArgs.append(arg)
        
        if len(unsetArgs) > 0:
            errMsgs.append("The following required arguments were not set:" + ','.join(unsetArgs))

        # Raise exception message if any issues occurred:
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))
