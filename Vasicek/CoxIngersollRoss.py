############################################
# CoxIngersolRoss.py
############################################
# Description:
# * Pricing using the CoxIngersolRoss model.

import math as m
from scipy.stats import norm
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
    def ZeroCouponBond(self, today, bondMaturity):
        """
        * Return price of zero coupon bond obeying the CoxIngersollRoss model.
        Inputs:
        * today: 
        """
        errMsgs = []
        if not InterestRatePricing.ValidTenor(today):
            errMsgs.append('today must be numeric and non-negative.')
        if not InterestRatePricing.ValidTenor(bondMaturity):
            errMsgs.append('bondMaturity must be numeric and non-negative.')
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))

        return self.Params.A(today, bondMaturity)

    def GenerateZeroCurve(self, tStart, tEnd, tStep):
        """
        * Generate zero coupon yield curve.
        Inputs:
        * tStart: Start year (numeric, non-negative).
        * tEnd: End year (numeric, non-negative).
        * tStep: Fraction of year step (numeric, positive).
        """
        return InterestRatePricing._InterestRatePricing__GenZeroCurve(self, tStart, tEnd, tStep)

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
        Inputs:
        * t_1: Expecting duratio 
        """
        return 2 * self.Alpha / (self.Sigma * self.Sigma * (1 - m.exp(-self.Alpha * (t_2 - t_1))))
    def epsilon(self, t_1, t_2):
        """
        * Return 
        """
        pass
    def B(self, t_1, t_2):
        """
        * Return B useful in ZCB price calculation.
        """
        return 2 * (m.exp(self.Gamma * (t_2 - t_1)) - 1) / self.__Denom(t_1, t_2)
    def A(self, t_1, t_2):
        """
        * Return 
        """
        num = 2 * self.Gamma * m.exp((self.Alpha + self.Lambda + self.Gamma) * (t_2 - t_1) / 2)
        return (num / self.__Denom(t_1, t_2)) ** (2 * self.Alpha * self.Mu / (self.Sigma * self.Sigma))
    def __Denom(self, t_1, t_2):
        """
        * Return denominator used in several formulas.
        """
        return (self.Alpha + self.Lambda + self.G(t_1, t_2)) * (m.exp(self.G(t_1, t_2) * (t_2 - t_1)) - 1) + 2 * self.Gamma
    ####################
    # Properties:
    ####################
    @property
    def Alpha_RN(self):
        return self.Alpha + self.Lambda
    @property
    def Mu_RN(self):
        out = self.Alpha * self.Mu
        return out / self.Alpha + self.Mu
    @property
    def Gamma(self):
        return m.sqrt(self.Alpha_RN * self.Alpha_RN + 2 * self.Sigma * self.Sigma)
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