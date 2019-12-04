############################################
# Vasicek.py
############################################
# Description:
# * 

import math as m
from scipy.stats import norm
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D

class VasicekPricing(object):
    """
    * Singleton class that implements pricing obeying the Vasicek interest rate model.
    """
    def __init__(self):
        """
        * Do nothing since object is a singleton.
        """
        pass
    ###############################
    # Static Methods:
    ###############################
    @staticmethod
    def PlotZeroCurve(params, termStruct):
        """
        * Plot the zero coupon yield curve (generated using GenerateZeroCurve()).
        Inputs:
        * params: Expecting VasicekParam object.
        * termStruct: Expecting dictionary mapping { T - t -> zero_coupon_yield }.
        """
        errMsgs = []
        if not isinstance(params, VasicekParam):
            errMsgs.append('')
        if not isinstance(termStruct, dict):
            errMsgs.append('termStruct must be a dictionary.')
        elif len(termStruct.keys()) == 0:
            errMsgs.append('termStruct must have at least one yield.')
        else:
            keyErr = False
            valErr = False
            for key in termStruct.keys():
                if not keyErr and not VasicekPricing.__IsNumeric(key):
                    errMsgs.append('All keys in termStruct must be numeric.')
                    keyErr = True
                if not valErr and not VasicekPricing.__IsNumeric(termStruct[key]):
                    errMsgs.append('All values in termStruct must be numeric.')
                    valErr = True
                if keyErr and valErr:
                    break
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))

        X = termStruct.keys()
        Y = termStruct.values()
        tEnd = max(list(X)) - 1
        title = ''.join(['Zero Curve (Vasicek Model)', params.ParamsString])
        tStart = min(list(X))
        fig = plotter.figure()
        fig .suptitle(title, fontsize = 10)
        axis = fig.add_subplot('111')
        axis.plot(X, Y)
        axis.set_ylabel('y(T - t)')
        axis.set_xlabel('T - t')

        plotObj.show()

        return plotObj

    @staticmethod
    def GenerateZeroCurve(params, tStart, tEnd, tStep):
        """
        * Generate zero coupon yield curve.
        Inputs:
        * tStart: Start year (numeric, non-negative).
        * tEnd: End year (numeric, non-negative).
        * tStep: Fraction of year step (numeric, positive).
        """
        errs = []
        if not isinstance(params, VasicekParam):
            errs.append('params must be a VasicekParam object.')
        if not VasicekPricing.__IsNumeric(tStart) and not tStart >= 0:
            errs.append("tStart must be a non-negative numeric value.")
        if not VasicekPricing.__IsNumeric(tEnd) and not tEnd >= 0 and not tEnd > tStart:
            errs.append("tEnd must be a non-negative numeric value greater than tStart.")
        if not VasicekPricing.__IsNumeric(tStep) and not tStep > 0:
            errs.append("tStep must be a positive numeric value.")
        if len(errs) > 0:
            raise Exception('\n'.join(errs))
        # Generate term structure { T - t -> zero_coupon_yield }:
        termStruct = {}
        origParams = params.Params.copy()
        while tStart < tEnd:
            params.t = tStart
            bondPrice = VasicekPricing.ZeroCouponBond(params)
            _yield = m.log(1 / bondPrice) / (tEnd - tStart)
            termStruct[tStart] = _yield
            tStart += tStep
        # Return passed parameters to original state:
        params = origParams
        
        return termStruct

    @staticmethod
    def ZeroCouponBond(params):
        """
        * Calculate price of zero coupon bond obeying Vasicek's mean reversion model.
        dr = alpha (mu - r) + sigma * dW
        Inputs:
        * params: Expecting a VasicekParam object.
        """
        VasicekPricing.__Validate(params)

        alpha = params.Alpha
        lambd = params.Lambda
        mu = params.Mu
        t = params.t
        T = params.T
        r = params.InstantaneousRate
        sig = params.Sigma
        # Calculate F(t, T): 
        f = 1 - m.exp(-alpha * (T - t))
        f /= alpha
        # Calculate G(t, T):
        g = (mu - sig * sig / (2 * alpha * alpha) - sig * lambd / alpha)
        g *= T - t - f
        g += sig * sig * f * f / (4 * alpha)

        return m.exp(-r * f - g)

    @staticmethod
    def ZeroCouponBondOption(p_s_t, p_t_T, strike, params):
        """
        * Calculate price of option on zero coupon bond obeying Vasicek's mean reversion model
        dr = alpha (mu - r) + sigma * dW
        Inputs:
        * p_s_t: Expecting price of bond at ... (numeric, positive).
        * p_t_T: Expecting price of bond at present (numeric, positive).
        * params: Expecting a VasicekParam object.
        """
        VasicekPricing.__Validate(params)
        alpha = params.Alpha
        lambd = params.Lambda
        mu = params.Mu
        t = params.t
        T = params.T
        r = params.InstantaneousRate
        sig = params.Sigma
        
        


    ###########################
    # Private Static Helpers:
    ###########################
    @staticmethod
    def __IsNumeric(val):
        """
        * Determine if value is numeric.
        """
        return isinstance(val, int) or isinstance(val, float)

    @staticmethod
    def __Validate(params):
        if not isinstance(params, VasicekParam):
            raise Exception("params must be a VasicekParam object.")
       
    @staticmethod
    def __StringToFunction(str):
        """
        * Map string to function.
        """
        pass
class VasicekParam(object):
    """
    * Object serves as parameter to Vasicek pricing static functions.
    """
    __args = {"alpha" : .001, "lambda" : 0, "mu" : 0, "r" : 0, "sigma" : 0, "T" : 0, "t" : 0}
    def __init__(self, argDict):
        """
        * Initiate new parameter object for VasicekModel object.
        Ensures that all parameter values are valid.
        Possible parameters are [alpha, lambda, mu, r, T, t].
        Call VasicekParams.ArgDict() to get copy of default kwargs dictionary.
        """
        # Ensure that all arguments have been passed:
        self.__req = {"alpha" : False, "lambda" : False, "mu" : False, "r" : False, "sigma" : False, "T" : False, "t" : False}
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
    def T(self):
        return self.__T
    @property
    def t(self):
        return self.__t
    @property
    def ParamsString(self):
        """
        * Return a string detailing the parameters.
        """
        params = self.Params
        return ','.join([key + ':' + str(params[key]) for key in params])
    @property
    def Params(self):
        """
        * Return dictionary containing this object's parameters.
        """
        return {"alpha" : self.Alpha, "lambda" : self.Lambda, "mu" : self.Mu, "r" : self.InstantaneousRate, "sigma" : self.Sigma, "T" : self.T, "t" : self.t}
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
    @T.setter
    def T(self, T_in):
        if not isinstance(T_in, int) and not isinstance(T_in, float):
            raise Exception("T must be numeric.")
        elif T_in < 0:
            raise Exception("T must be non-negative.")
        elif T_in < self.t:
            raise Exception("T must be >= t.")
        self.__req['T'] = True
        self.__T = T_in
    @t.setter
    def t(self, t_in):
        if not isinstance(t_in, int) and not isinstance(t_in, float):
            raise Exception("t must be numeric.")
        elif t_in < 0:
            raise Exception("t must be non-negative.")
        self.__req['t'] = True
        self.__t = t_in
    @Params.setter
    def Params(self, params):
        if not isinstance(params, dict):
            raise Exception("Params must be a dictionary.")
        for key in params.keys():
            self.__SetAttribute(key, params[key])
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
        if paramStr == 'T' or paramStr == 't':
            return paramStr
        
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
