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
    argDict['T'] = 30
    argDict['t'] = 0
    params = VasicekParam(argDict)
    curve = VasicekPricing.GenerateZeroCurve(params, 1, 30, 1)
    plot = VasicekPricing.PlotZeroCurve(params, curve)


