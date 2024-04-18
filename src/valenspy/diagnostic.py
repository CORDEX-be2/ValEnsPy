from datatree import DataTree
import xarray as xr

from abc import abstractmethod


class Diagnostic:
    """An abstract class representing a diagnostic."""

    def __init__(
        self, diagnostic_function, visualization_function, name=None, description=None
    ):
        """Initialize the Diagnostic."""
        self.name = name
        self.description = description
        self.diagnostic_function = diagnostic_function
        self.visualization_function = visualization_function

    @abstractmethod
    def apply(self, data: DataTree):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data : DataTree
            The data to apply the diagnostic to.

        Returns
        -------
        Dataset
            The data after applying the diagnostic either as a DataTree, Dataset, DataArray, Scalar, or a pandas DataFrame.
        """
        pass

    @abstractmethod
    def visualize(self):
        """Visualize the diagnostic.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        pass


class Model2Ref(Diagnostic):
    """A class representing a diagnostic that compares a model to a reference."""

    def __init__(
        self, diagnostic_function, visualization_function, name=None, description=None
    ):
        """Initialize the Model2Ref diagnostic."""
        super().__init__(diagnostic_function, visualization_function, name, description)

    def apply(self, data: xr.Dataset, ref: xr.Dataset, **kwargs):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data : xr.Dataset
            The data to apply the diagnostic to.
        ref : xr.Dataset
            The reference data to compare the data to.

        Returns
        -------
        xr.Dataset
            The data after applying the diagnostic.
        """
        return self.diagnostic_function(data, ref, **kwargs)

    def visualize(self, data: xr.Dataset, ref: xr.Dataset):
        """Visualize the diagnostic.

        Parameters
        ----------
        data : xr.Dataset
            The data to visualize.
        ref : xr.Dataset
            The reference data to compare the data to.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        return self.visualization_function(self.apply(data, ref))


class Ensemble2Ref(Diagnostic):
    """A class representing a diagnostic that compares an ensemble to a reference."""

    def __init__(
        self, diagnostic_function, visualization_function, name=None, description=None
    ):
        """Initialize the Ensemble2Ref diagnostic."""
        super().__init__(diagnostic_function, visualization_function, name, description)

    def apply(self, data: xr.Dataset, ref: xr.Dataset):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data : xr.Dataset
            The data to apply the diagnostic to.
        ref : xr.DataArray
            The reference data to compare the data to.

        Returns
        -------
        xr.Dataset
            The data after applying the diagnostic.
        """
        return self.diagnostic_function(data, ref)

    def visualize(self, data: xr.Dataset, ref: xr.Dataset):
        """Visualize the diagnostic.

        Parameters
        ----------
        data : xr.Dataset
            The data to visualize.
        ref : xr.Dataset
            The reference data to compare the data to.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        return self.visualization_function(self.apply(data, ref))
