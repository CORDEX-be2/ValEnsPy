import xarray as xr
import dask
import shapely
import os

class Modeldata:
    """
    A class to represent a single netCDF file which is in CF convention. 
    The file is loaded as a dask xarray dataset.
    """

    def __init__(self, file_location):
        """
        Initialize modeldata object from a location of a netCDF file.
        
        Parameters
        ----------
        file_location : str
            The location of the netCDF file.
        """
        self.file_location = file_location
        self.load_dataset()

    def __str__(self):
        """Return a string representation of the dataset."""
        if not self.ds:
            return f"{self.file_location} is not a valid netCDF file."
        start_dt = self.ds.time.min().values
        end_dt = self.ds.time.max().values
        variables = ", ".join(self.ds.data_vars)
        domain = self.get_domain()
        return f"Modeldata object: {start_dt} - {end_dt}, variables: {variables}, domain: {domain.bounds}"
    
    def __repr__(self):
        """Return a string representation of the dataset."""
        return self.__str__()
    
    def load_dataset(self):
        """
        Load the dataset as a dask xarray dataset and store it in the ds attribute.
        A check is performed to see if the file is a valid netCDF file in CF convention.
        
        Raises
        ------
        FileNotFoundError
            If the file_location is not a valid file location.
        ValueError
            If the file is not in CF convention.
        """
        if not os.path.isfile(self.file_location):
            raise FileNotFoundError(f"{self.file_location} is not a valid file location.")
        else:
            self.ds = xr.open_dataset(self.file_location, chunks={'time': 'auto'})
            if not self._is_CF_convention():
                raise ValueError(f"{self.file_location} is not in the CF convention.")

    @property
    def domain_bound(self):
        """Property: The outer bounds of the domain of the dataset as a shapely polygon."""
        if not self.ds:
            return None
        return shapely.geometry.box(self.ds.lon.min(), self.ds.lat.min(), self.ds.lon.max(), self.ds.lat.max())
    
    def _is_same_ensmember(self, other):
        """Check if the modeldata object is from the same ens_member with another modeldata object."""
        #TODO: Implement this method to check if the modeldata objects are possibly from the same ensemble member
        #Based on??
        #No overlapping time periods (for the same variable), the same domain, ... ?
        #These two modeldatasets should (at least) be concatenatable along the time dimension.
        return True
    
    def _is_CF_convention(self):
        """Check if the dataset is in CF convention."""
        #TODO: Expand this method to check if the dataset is in CF convention
        return True
