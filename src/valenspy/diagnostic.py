from datatree import DataTree
import xarray as xr

from abc import abstractmethod

class Diagnostic:
    """An abstract class representing a diagnostic."""

    def __init__(self,diagnostic_function, visualization_function, name=None, description=None):
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
        

class Model2Obs(Diagnostic):
    """A class representing a diagnostic that compares a model to observations."""

    def __init__(self,diagnostic_function, visualization_function, name=None, description=None):
        """Initialize the Model2Obs diagnostic."""
        super().__init__(diagnostic_function,visualization_function, name, description)

    def apply(self, data: xr.Dataset, obs: xr.Dataset):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data : DataTree
            The data to apply the diagnostic to.
        obs : xr.Dataset
            The observations to compare the data to.

        Returns
        -------
        Dataset
            The data after applying the diagnostic either as a DataTree, Dataset, DataArray, Scalar, or a pandas DataFrame.
        """
        return self.diagnostic_function(data, obs)

    def visualize(self, data: xr.Dataset, obs: xr.Dataset):
        """Visualize the diagnostic.

        Parameters
        ----------
        data : DataTree
            The data to visualize.
        obs : xr.Dataset
            The observations to compare the data to.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        return self.visualization_function(self.apply(data, obs))

    
class Ensemble2Obs(Diagnostic):
    """A class representing a diagnostic that compares an ensemble to observations."""

    def __init__(self,diagnostic_function, name=None, description=None):
        """Initialize the Ensemble2Obs diagnostic."""
        super().__init__(diagnostic_function,visualization_function, name, description)

    def apply(self, data: DataTree, obs: xr.Dataset):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data : DataTree
            The data to apply the diagnostic to.
        obs : xr.Dataset
            The observations to compare the data to.

        Returns
        -------
        Dataset
            The data after applying the diagnostic either as a DataTree, Dataset, DataArray, Scalar, or a pandas DataFrame.
        """
        return self.diagnostic_function(data, obs)
    
