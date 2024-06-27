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
        self._description = description
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

    @property
    def description(self):
        """Return the description of the diagnostic a combination of the name, the type and the description and the docstring of the diagnostic and plot functions."""
        return f"{self.name} ({self.__class__.__name__})\n{self._description}\n Diagnostic function: {self.diagnostic_function.__name__}\n {self.diagnostic_function.__doc__}\n Visualization function: {self.visualization_function.__name__}\n {self.visualization_function.__doc__}"

class Model2Self(Diagnostic):
    """A class representing a diagnostic that compares a model to itself."""

    def __init__(
        self, diagnostic_function, visualization_function, name=None, description=None
    ):
        """Initialize the Model2Self diagnostic."""
        super().__init__(diagnostic_function, visualization_function, name, description)
    
    def apply(self, data: xr.Dataset, **kwargs):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data : xr.Dataset
            The data to apply the diagnostic to.

        Returns
        -------
        xr.Dataset
            The data after applying the diagnostic.
        """
        return self.diagnostic_function(data, **kwargs)

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

# =============================================================================
# Pre-made diagnostics
# =============================================================================

from valenspy.diagnostic_functions import daily_cycle, spatial_bias, time_series_spatial_mean, temporal_bias, daily_cycle_bias
from valenspy.diagnostic_visualizations import plot_daily_cycle, plot_spatial_bias, plot_time_series

#Model2Self diagnostics
vp_DailyCycle = Model2Self(daily_cycle, plot_daily_cycle, "Daily Cycle", "The daily cycle of the data.")
vp_TimeSeriesSpatialMean = Model2Self(time_series_spatial_mean, plot_time_series, "Time Series Spatial Mean", "The time series of the spatial mean of the data.")
#Model2Ref diagnostics
vp_SpatialBias = Model2Ref(spatial_bias, plot_spatial_bias, "Spatial Bias", "The spatial bias of the data compared to the reference.")
vp_TemporalBias = Model2Ref(temporal_bias, plot_time_series, "Temporal Bias", "The temporal bias of the data compared to the reference.")
vp_DailyCycleBias = Model2Ref(daily_cycle_bias, plot_daily_cycle, "Daily Cycle Bias", "The daily cycle bias of the data compared to the reference.")
