class Ensmember:
    """A class representing an ensemble member which consists of multiple modeldata objects."""

    def __init__(self, data=[], experiment=None , institution=None, model=None):
        """Initialize an ensemble member

        Parameters
        ----------
        experiment : str
            The experiment that the ensemble member belongs to. e.g. 'historical', 'rcp85', 'ssp585'
        data : list of Modeldata objects
            The modeldata objects that make up the ensemble member
        institution : str
            The institution that the ensemble member belongs to. e.g. "COSMO-CLM"
        model : str
            The model name of the ensemble member. e.g. "ALARO-1"
        """
        self.data = data
        self.experiment = experiment
        self.institution = institution
        self.model = model

    def __str__(self):
        return f"Ensemble member: \n model: {self.model} \n institution: {self.institution} \n experiment: {self.experiment} \n data: {self.data}"

    def __repr__(self):
        return self.__str__()
    
    def add_modeldata(self, modeldata):
        """Add a modeldata object to the ensemble member."""
        self.data.append(modeldata)

    def _is_consistent(self, other):
        """Check if the ensemble member is consistent with another ensemble member."""
        #TODO: Implement this method to check if the ensemble members are consistent - i.e. cover the same time period, domain and variables
        return True

    def _check_modeldata(self, modeldata):
        """Check if the modeldata objects in the ensemble member are consistent."""
        if not self.data:
            return False
        return modeldata._is_consistent(self.data[0])

