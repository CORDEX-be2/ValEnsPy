from valenspy.preprocessing_tasks.task import PreprocessingTask
import xarray as xr


class Regrid(PreprocessingTask):
    REGRID_METHODS_XARRAY = {"linear", "nearest"} #Native xarray regridding methods (for n-dimensional data)
    """A regrid preprocessing task."""

    def __init__(self, target_grid: xr.Dataset, name="", description=None):
        """Initialize the Regrid task."""
        super().__init__("regrid_" + name, description)
        
        self.target_grid = target_grid

    def apply(self, data: xr.Dataset, method="linear") -> xr.Dataset:
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

        if method not in self.REGRID_METHODS_XARRAY:
            raise NotImplementedError(
                f"{method} regridding is not yet implemented. Please choose from {self.REGRID_METHODS_XARRAY}."
            )

        return data.interp(
            lon=self.target_grid.lon, lat=self.target_grid.lat, method=method
        )
