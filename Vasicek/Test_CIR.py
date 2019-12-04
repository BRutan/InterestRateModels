####################################
# Test_CIR.py
####################################
# Description:
# * Test the CoxIngersollRoss class.

from CoxIngersollRoss import CIRParams, CoxIngersollRossPricing

if __name__ == '__main__':

    args = CIRParams.ArgDict()
    args['alpha'] = .2456
    args['mu'] = .0648
    args['sigma'] = .14998
    args['lambda'] = -.129
    args['r'] = .06
    params = CIRParams(args)
    model = CoxIngersollRossPricing(params)
    curve = model.GenerateZeroCurve(0, 30, 1)
    plot = model.PlotCurve(curve)

    futuresFV = model.ZCBFuturesFV(0, 1, 5)
    optFV = model.ZCBOptionFV(0, .5, 5, .6)
    couponBond = model.CouponBondFV(.05, 1, 10, 10, 1)
    coupBondOption = model.CouponBondOptionFV(0, .05, 1, 1, .5, 5)
    coupFuture = model.CouponBondFutureFV(0, .05, 1, 5, 1, 3, 5)
    
    print("Cox Ingersoll Ross Model")
    print(params.ParamsString)
    print("Futures FV: " + str(futuresFV))
    print("Annual Coupon Bond FV (5%, $1 Face Value, 10 Years To Maturity): " + str(couponBond))
    print("Coupon Bond Option FV: " + str(coupBondOption))
    print("Coupon Bond Future FV: " + str(coupFuture))
