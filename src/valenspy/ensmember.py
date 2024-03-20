import xarray as xr
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
        self.load_data()

    def __str__(self):
        return f"Ensemble member: \n model: {self.model} \n institution: {self.institution} \n experiment: {self.experiment} \n run_id: {self.run_id} \n time period: {self.time_period} \n resolution: {self.resolution} \n domain: {self.get_domain}"

    def __repr__(self):
        return self.__str__()
    
    def load_data(self):
        """Load the data from the modeldata objects in the ensemble member."""
        self.ds = xr.open_mfdataset([modeldata.file_location for modeldata in self.data], combine='by_coords', chunks='auto')
    
    @property
    def time_period(self):
        """Return the time period of the ensemble member."""
        return self.ds.time.min().values, self.ds.time.max().values
    
    @property
    def get_domain(self):
        """Return the domain of the ensemble member."""
        return self.data[0].domain_bound
    
    @property
    def resolution(self):
        """Return the resolution of the ensemble member as the nominal resolution attribute of the ds."""
        return self.ds.attrs['nominal_resolution']
    
    def add_modeldata(self, modeldata):
        """Add a modeldata object to the ensemble member."""
        self.data.append(modeldata)

    def _is_consistent(self):
        #Check if the xarray contains data for each variable, for the full time period and the same domain (lat, lon and resolution)
        return True

