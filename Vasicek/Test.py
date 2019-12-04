####################################
# Test.py (Vasicek)
####################################
# Description:
# * Test the singleton Vasicek
# class.

from Vasicek import VasicekParam, VasicekPricing

if __name__ == '__main__':
    """
    * Test all functions in Vasicek class.
    """
    argDict = VasicekParam.ArgDict()
    argDict['alpha'] = .2456
    argDict['mu'] = .0648
    argDict['sigma'] = .0289
    argDict['lambda'] = -.2718
    argDict['r'] = .06
    params = VasicekParam(argDict)
    curve = VasicekPricing.GenerateZeroCurve(params, 1, 30, 1)
    plot = VasicekPricing.PlotZeroCurve(params, curve)
    optPrice = VasicekPricing.ZCBOptionPrice(params, .6, 0, .5, 5)
    futPrice = VasicekPricing.ZCBFuturesPrice(params, 0, 1, 5)
    futOptPrice = VasicekPricing.ZCBFuturesOptionPrice(params, .5, 0, 1, 5, 6)
    print("Vasicek Model Results:")
    print("Parameters: {" + params.ParamsString + "}")
    print("Option Price [T_Bond = 5, K = .6, T_Opt = .5]: " + str(optPrice))
    print("Future Price [T_Bond = 5, T_Fut = 1]: " + str(futPrice))
    print("Futures Option Price [T_Bond = 5, K = .5, T_Opt = 1, T_Fut = 5]: " + str(futOptPrice))