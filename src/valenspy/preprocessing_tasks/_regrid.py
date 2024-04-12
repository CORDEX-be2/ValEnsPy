from valenspy.preprocessing_tasks.task import PreprocessingTask
import xarray as xr


class Regrid(PreprocessingTask):
    """A regrid preprocessing task."""

    def __init__(self, target_grid: xr.Dataset, name="", description=None):
        """Initialize the Regrid task."""
        super().__init__("regrid_" + name, description)
        # ToDo: Extend target grid to more options - kwargs should be passed and saved here.
        self.target_grid = target_grid

    def apply(self, data: xr.Dataset):
        """Apply the regridding task to an xr.Dataset.

        Parameters
        ----------
        data : xr.Dataset
            The data to apply the task to.
        target_grid : xr.Dataset
            The target grid to regrid the data to.

        Returns
        -------
        xr.Dataset
            The regridded dataset.
        """
        # ToDo: Implement regridding here
        # https://github.com/pangeo-data/xESMF from Pangeo for regridding capabilities (including conservative regridding)!

        return data.interp(
            lon=self.target_grid.lon, lat=self.target_grid.lat, method="linear"
        )
