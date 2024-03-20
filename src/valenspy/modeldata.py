import xarray as xr
import dask
import shapely

class Modeldata:
    """A class to represent a single xarray dataset in CF convention."""

    def __init__(self, file_location=None):
        """Initialize an empty dataset."""
        self.file_location = file_location
        self.load_dataset()

    def __str__(self):
        """Return a string representation of the dataset."""
        if self.ds:
            return "Empty dataset"
        start_dt = self.ds.time.min().values
        end_dt = self.ds.time.max().values
        variables = ", ".join(self.ds.data_vars)
        domain = self.get_domain()
        return f"Modeldata object: {start_dt} - {end_dt}, variables: {variables}, domain: {domain.bounds}"
    
    def __repr__(self):
        """Return a string representation of the dataset."""
        return self.__str__()
    
    def load_dataset(self):
        """Load the dataset as a dask xarray dataset."""
        if self.file_location:
            self.ds = xr.open_dataset(self.file_location, chunks={'time': 'auto'})
        else:
            self.ds = None

    def get_domain(self):
        """Returns a shapely polygon of the domain of the dataset."""
        if not self.ds:
            return None
        return shapely.geometry.box(self.ds.lon.min(), self.ds.lat.min(), self.ds.lon.max(), self.ds.lat.max())
    
    def _is_consistent(self, other):
        """Check if the modeldata object is consistent with another modeldata object."""
        if self.ds.empty or other.ds.empty:
            return False
        same_domain = self.get_domain().equals(other.get_domain())
        return same_domain
    
    def _check_CF_convention(self):
        """Check if the dataset is in CF convention."""
        #TODO: Expand this method to check if the dataset is in CF convention

        # Check if the dataset has the required variables
        required_vars = ['time', 'lat', 'lon']
        for var in required_vars:
            if var not in self.ds.variables:
                return False

        # Check if the dataset has the required attributes
        required_attrs = ['standard_name', 'units']
        for var in self.ds.variables:
            for attr in required_attrs:
                if attr not in self.ds[var].attrs:
                    return False

        return True
