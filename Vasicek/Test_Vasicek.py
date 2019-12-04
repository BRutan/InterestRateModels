####################################
# Test_Vasicek.py
####################################
# Description:
# * Test the Vasicek class.

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
    model = VasicekPricing(params)
    curve = model.GenerateZeroCurve(0, 30, 1)
    plot = model.PlotCurve(curve)
    optFV = model.ZCBOptionFV(.6, 0, .5, 5)
    futFV = model.ZCBFuturesFV(0, 1, 5)
    futOptFV = model.ZCBFuturesOptionFV(.5, 0, 1, 5, 6)
    print("Vasicek Model Results:")
    print("Parameters: {" + params.ParamsString + "}")
    print("Option Fair Value [T_Bond = 5, K = .6, T_Opt = .5]: " + str(optFV))
    print("Future Fair Value[T_Bond = 5, T_Fut = 1]: " + str(futFV))
    print("Futures Option Fair Value [T_Bond = 5, K = .5, T_Opt = 1, T_Fut = 5]: " + str(futOptFV))