############################################
# InterestRatePricing.py
############################################
# Description:
# * Define abstract base class with methods 
# to use in derived pricing models.

from abc import ABC, abstractmethod, abstractproperty
import math as m
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D

class InterestRatePricing(ABC):
    """
    * Abstract base class with abstract methods
    for use in derived interest rate pricing 
    classes.
    """
    ########################
    # Constructors:
    ########################
    def __init__(self):
        pass
    
    ########################
    # Properties:
    ########################
    @property
    def ModelName(self):
        return self.__modelName
    @ModelName.setter
    def ModelName(self, name):
        if not isinstance(name, str):
            raise Exception('ModelName must be a string.')
        self.__modelName = name.strip()
    @property
    def Params(self):
        return self.__params
    @Params.setter
    def Params(self, value):
        self.__params = value
    ########################
    # Methods:
    ########################
    @abstractmethod
    def ZeroCouponBond(self, today, bondMaturity):
        pass

    def PlotZeroCurve(self, termStruct):
        """
        * Plot the zero coupon yield curve (generated using GenerateZeroCurve()).
        Inputs:
        * params: Expecting VasicekParam object.
        * termStruct: Expecting dictionary mapping { T - t -> zero_coupon_yield }.
        """
        errMsgs = []
        if not isinstance(termStruct, dict):
            errMsgs.append('termStruct must be a dictionary.')
        elif len(termStruct.keys()) == 0:
            errMsgs.append('termStruct must have at least one yield.')
        else:
            keyErr = False
            valErr = False
            for key in termStruct.keys():
                if not keyErr and not InterestRatePricing.IsNumeric(key):
                    errMsgs.append('All keys in termStruct must be numeric.')
                    keyErr = True
                if not valErr and not InterestRatePricing.IsNumeric(termStruct[key]):
                    errMsgs.append('All values in termStruct must be numeric.')
                    valErr = True
                if keyErr and valErr:
                    break
        if len(errMsgs) > 0:
            raise Exception('\n'.join(errMsgs))

        params = self.Params
        X = list(termStruct.keys())
        Y = list(termStruct.values())
        tEnd = max(X) - 1
        title = ''.join(['Zero Curve (', self.ModelName, '){', params.ParamsString, '}'])
        tStart = min(X)
        fig = plotter.figure()
        fig .suptitle(title, fontsize = 8)
        axis = fig.add_subplot('111')
        axis.plot(X, Y)
        axis.set_ylabel('y(T - t)')
        axis.set_xlabel('T - t')

        fig.show()

        return fig

    ###########################
    # Static Helpers:
    ###########################
    @staticmethod
    def IsNumeric(val):
        """
        * Determine if value is numeric.
        """
        return isinstance(val, int) or isinstance(val, float)

    @staticmethod
    def ValidTenor(val):
        """
        * Determine if value corresponds to a valid tenor (numeric, non-negative).
        """
        return InterestRatePricing.IsNumeric(val) and val >= 0

    @staticmethod
    def __GenZeroCurve(model, tStart, tEnd, tStep):
        """
        * Generate zero coupon yield curve.
        Inputs:
        * model: Expecting object derived from InterestRatePricing.
        * tStart: Start year (numeric, non-negative).
        * tEnd: End year (numeric, non-negative).
        * tStep: Fraction of year step (numeric, positive).
        """
        errs = []
        if not isinstance(model, InterestRatePricing):
            errs.append('model must be derived from InterestRatePricing.')
        if not InterestRatePricing.IsNumeric(tStart) and not tStart >= 0:
            errs.append("tStart must be a non-negative numeric value.")
        if not InterestRatePricing.IsNumeric(tEnd) and not tEnd >= 0 and not tEnd > tStart:
            errs.append("tEnd must be a non-negative numeric value greater than tStart.")
        if not InterestRatePricing.IsNumeric(tStep) and not tStep > 0:
            errs.append("tStep must be a positive numeric value.")
        if len(errs) > 0:
            raise Exception('\n'.join(errs))

        # Generate term structure { T - t -> zero_coupon_yield }:
        termStruct = {}
        while tStart < tEnd:
            bondPrice = model.ZeroCouponBond(tStart, tEnd)
            _yield = m.log(1 / bondPrice) / (tEnd - tStart)
            termStruct[tEnd - tStart] = _yield
            tStart += tStep
        
        return termStruct

class ParamsObj(ABC):
    """
    * Abstract base class for all parameter objects for
    each model.
    """
    def __init__(self, params):
        self.Params = params
    @abstractproperty 
    def Params(self):
        pass
    @abstractmethod
    def Params(self, params):
        pass
    @property
    def ParamsString(self):
        """
        * Return a string detailing the parameters.
        """
        params = self.Params
        return ','.join([key + ' :{0:.2f}'.format(params[key]) for key in params])
