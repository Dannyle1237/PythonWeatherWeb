__author__ = 'robswift'

class Molecule:
    """
	csv molecule class
	"""

    def __init__(self):
        """
		class attributes
		"""
        self.id = None
        self.status = None
        self.bm = None
        self.graph = None
        self.scores = {}
        self.properties = {}
        self.tags = ('id', 'status', 'bm', 'graph')

    def SetProp(self, prop, value, score=None):
        """
		set attribute
		"""
        if score:
            self.scores[prop] = value
        else:
            if prop == 'id':
                self.id = value
            elif prop == 'status':
                self.status = value
            elif prop == 'bm':
                self.bm = format(value)
            elif prop == 'graph':
                self.graph = format(value)
            else:
                self.properties[prop] = value

    def GetProp(self, prop):
        """
		get attribute
		"""

        if prop == 'scores':
            return [x[1] for x in self.scores.items()]
        elif prop:
            if prop == 'id':
                return self.id
            elif prop == 'status':
                return self.status
            elif prop == 'bm':
                return self.bm
            elif prop == 'graph':
                return self.graph
            elif prop in self.scores.keys():
                return self.scores[prop]
            elif prop in self.properties.keys():
                return self.properties[prop]
            else:
                return None