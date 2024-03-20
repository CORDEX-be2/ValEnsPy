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

    @property
    def domain_bound(self):
        """Returns a shapely polygon of the domain of the dataset."""
        if not self.ds:
            return None
        return shapely.geometry.box(self.ds.lon.min(), self.ds.lat.min(), self.ds.lon.max(), self.ds.lat.max())
    
    def _is_same_ensmember(self, other):
        """Check if the modeldata object is from the same ens_member with another modeldata object."""
        if not self.ds or not other.ds:
            return False
        if not self.domain_bound.equals(other.domain_bound):
            return False
        for var in self.ds.data_vars:
            if var in other.ds.data_vars:
                if self.ds[var].time.min() > other.ds[var].time.max() or self.ds[var].time.max() < other.ds[var].time.min():
                    return False
        return True
    
    def _is_CF_convention(self):
        """Check if the dataset is in CF convention."""
        #TODO: Expand this method to check if the dataset is in CF convention
        return True
