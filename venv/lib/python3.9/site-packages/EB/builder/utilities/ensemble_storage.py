__author__ = 'robswift'

class EnsembleStorage:
    """
    Stores an ensemble's composition and performance
    """

    def __init__(self):
        """
        class attributes
        :return:
        """
        self.ensemble = None    # a tuple: ('receptor_x', ...., 'receptor_y')
        self.auc = None         # a float: auc value
        self.ef = {}            # a dictionary: {'fpf value' : float(enrichment factor)}

    def set_prop(self, prop, value, ef=None):
        """
        set attributes values
        :param prop:
        :param value:
        :param ef:
        :return:
        """
        if ef:
            # prop should be restricted to n_decoys, an int, the no. of decoys corresponding to a given FPF.
            # value is restricted to the corresponding enrichment factor and should be a float
            self.ef[prop] = value
        else:
            if prop == 'ensemble':
                # value is a tuple of strings that gives the ensemble composition
                self.ensemble = value
            elif prop == 'auc':
                # value is a float that gives the auc value
                self.auc = value

    def get_prop(self, prop, ef=None):
        """
        return property attributes
        :param prop:
        :return:
        """
        if ef:
            return self.ef[prop]
        else:
            if prop == 'auc':
                return self.auc
            elif prop == 'ensemble':
                return self.ensemble