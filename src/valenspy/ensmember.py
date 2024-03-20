class Ensmember:
    """A class representing one ensemble member which consists of multiple modeldata objects."""

    def __init__(self, data=[], experiment=None, institution=None, model=None, run_id=None):
        """Initialize an ensemble member

        Parameters
        ----------
        data : list of Modeldata objects
            The modeldata objects that make up the ensemble member
        experiment : str
            The experiment that the ensemble member belongs to. e.g. 'historical', 'rcp85', 'ssp585'
        institution : str
            The institution that the ensemble member belongs to. e.g. "COSMO-CLM"
        model : str
            The model name of the ensemble member. e.g. "ALARO-1"
        run_id : str
            The run id of the ensemble member. e.g. "r1i1p1f1"
        """
        self.data = data
        self.experiment = experiment
        self.institution = institution
        self.model = model
        self.run_id = run_id

    def __str__(self):
        return f"Ensemble member: \n model: {self.model} \n institution: {self.institution} \n experiment: {self.experiment} \n data: {self.data}"

    def __repr__(self):
        return self.__str__()
    
    
    
    def add_modeldata(self, modeldata):
        """Add a modeldata object to the ensemble member."""
        self.data.append(modeldata)

    def _is_consistent(self):
        #To be improved (maybe using dask and xarray to try to concatenate the data and see if it works) i.e. no overlapping time periods (for different vars) or different domains
        """Check if the ensemble member is consistent, i.e. the modeldata covers the same time period, domain and have the same resolution."""

        
        return True

