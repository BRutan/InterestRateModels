############################################
# CoxIngersolRoss.py
############################################
# Description:
# * Pricing using the CoxIngersolRoss model.

class CoxIngersolRoss(object):
    """
    * Singleton class containing pricing methods 
    obeying the CIR model.
    """
    pass



class CIRParams(object):
    """
    * Object stores parameters used in the pricing methods in
    CoxIngersolRoss object.
    """
    __args = {'alpha' : .001, 'mu' : 0, 'r' : 0, 'sigma' : 0}
    def __init__(self, params):
        pass
    

    ####################
    # Static Methods:
    ####################
    @staticmethod
    def ArgDict():
        """
        * Return copy of default params dictionary used in 
        static pricing methods.
        """
        return CIRParams.__args.copy()

