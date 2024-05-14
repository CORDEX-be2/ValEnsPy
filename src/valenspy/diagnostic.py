from datatree import DataTree
import xarray as xr
import matplotlib.pyplot as plt

from abc import abstractmethod


class Diagnostic:
    """An abstract class representing a diagnostic."""

    def __init__(
        self, diagnostic_function, visualization_function, name=None, description=None
    ):
        """Initialize the Diagnostic.

        Parameters
        ----------
        diagnostic_function
            The function that applies a diagnostic to the data.
        visualization_function
            The function that visualizes the results of the diagnostic.
        name : str
            The name of the diagnostic.
        description : str
            The description of the diagnostic.
        """
        self.name = name
        self.description = description
        self.diagnostic_function = diagnostic_function
        self.visualization_function = visualization_function

    @abstractmethod
    def apply(self, data):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data
            The data to apply the diagnostic to.

        Returns
        -------
        Results
            The data after applying the diagnostic either as a DataTree, Dataset, DataArray, Scalar, or a pandas DataFrame.
        """
        pass

    @abstractmethod
    def visualize(self, result):
        """Visualize the diagnostic.
        Parameters
        ----------
        Result
            The output of the diagnostic function.

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

    def visualize(self, result, ax=None, **kwargs):
        """Visualize the diagnostic.

        Parameters
        ----------
        result :
            The output of the diagnostic function.

        Returns
        -------
        Figure :
            The figure representing the diagnostic.
        """
        if ax is None:
            ax = plt.gca()
        if isinstance(result, tuple):
            ax = self.visualization_function(*result, ax=ax, **kwargs)
        else:
            ax = self.visualization_function(result, ax=ax, **kwargs)
        return ax


class Ensemble2Ref(Diagnostic):
    """A class representing a diagnostic that compares an ensemble to a reference."""

    def __init__(
        self, diagnostic_function, visualization_function, name=None, description=None
    ):
        """Initialize the Ensemble2Ref diagnostic."""
        super().__init__(diagnostic_function, visualization_function, name, description)

    def apply(self, data: xr.Dataset, ref: xr.Dataset, **kwargs):
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
        return self.diagnostic_function(data, ref, **kwargs)

    def visualize(self, result, axes=None, facetted=True, **kwargs):
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
        if axes is None:
            if facetted:
                fig, axes = plt.subplots(1, len(result), figsize=(5 * len(result), 5))
            else:
                ax = plt.gca()
        return self.visualization_function(
            result, axes=axes, facetted=facetted, **kwargs
        )

    @classmethod
    def from_model2ref(cls, model2ref: Model2Ref, facetted=True):
        """Create an Ensemble2Ref diagnostic from a Model2Ref diagnostic.

        Parameters
        ----------
        model2ref : Model2Ref
            The Model2Ref diagnostic to convert.

        Returns
        -------
        Ensemble2Ref
            The Ensemble2Ref diagnostic.
        """

        def diagnostic_function(dt: DataTree, ref, **kwargs):
            ensemble_results = {}
            if isinstance(ref, DataTree):
                for data_node, ref_node in zip(dt.leaves, ref.leaves):
                    ensemble_results[data_node.path] = model2ref.diagnostic_function(
                        data_node.ds, ref_node.ds, **kwargs
                    )
            else:
                for data_node in dt.leaves:
                    ensemble_results[data_node.path] = model2ref.diagnostic_function(
                        data_node.ds, ref, **kwargs
                    )
            return ensemble_results

        def visualization_function(results, axes, facetted=facetted, **kwargs):
            if facetted:
                for path, result, ax in zip(
                    results.keys(), results.values(), axes.flatten()
                ):
                    model2ref.visualize(result, ax=ax, **kwargs)
                    ax.set_title(path.replace("/", " "))
            else:
                for path, result in results.items():
                    model2ref.visualize(
                        result, ax=axes, label=f'{path.replace("/", " ")}', **kwargs
                    )
            return axes

        return Ensemble2Ref(
            diagnostic_function,
            visualization_function,
            model2ref.name,
            model2ref.description,
        )
