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


