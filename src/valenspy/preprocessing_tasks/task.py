import xarray as xr

class Task():
    """Base class for all tasks."""

    def __init__(self, name=None, description=None):
        """Initialize the Task."""
        self.name = name
        self.description = description

    def apply(self, data: xr.Dataset):
        """Apply the task to the data.

        Parameters
        ----------
        data : xr.Dataset
            The data to apply the task to.

        Returns
        -------
        xr.Dataset
            The data after applying the task.
        """
        pass

class PreprocessingTask(Task):
    """A class representing a preprocessing task."""

    def __init__(self, name=None, description=None):
        """Initialize the PreprocessingTask."""
        super().__init__("pre_processing_" + name, description)

    def apply(self, data: xr.Dataset):
        """Apply the preprocessing task to the data.

        Parameters
        ----------
        data : xr.Dataset
            The data to apply the task to.

        Returns
        -------
        xr.Dataset
            The data after applying the task.
        """
        pass

    def _update_metadata(self, data: xr.Dataset):
        """Update the metadata of the data after applying the task.

        Parameters
        ----------
        data : xr.Dataset
            The data to update the metadata of.

        Returns
        -------
        xr.Dataset
            The data with updated metadata.
        """
        pass