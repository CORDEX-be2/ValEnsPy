from valenspy.preprocessing_tasks.task import PreprocessingTask
import xarray as xr

class Set_Domain(PreprocessingTask):
    """Selection of a common area for all grids."""

    #TODO extend this to accept the target area as a dataset, or any other format that can be used to select the area.
    def __init__(self, target_domain: xr.Dataset, name="", description=None):
        """Initialize the Regrid task."""
        super().__init__("regrid_" + name, description)
        #ToDo: Extend target grid to more options - kwargs should be passed and saved here.
        self.target_domain = target_domain

    def apply(self, data: xr.Dataset):
        """Apply the sub-selection of a common area on the data.

        Parameters
        ----------
        data : xr.Dataset
            The data to apply the task to.
        target_grid : xr.Dataset
            The target dataset to which the data should be selected.

        Returns
        -------
        xr.Dataset
            The regridded dataset.
        """
        #ToDo: Implement selection of common area here!

        xmin, ymin, xmax, ymax = self.target_domain.lon.min(), self.target_domain.lat.min(), self.target_domain.lon.max(), self.target_domain.lat.max()
        return data.sel(lon=slice(xmin, xmax), lat=slice(ymin, ymax))