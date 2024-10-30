from datatree import DataTree
import xarray as xr
import matplotlib.pyplot as plt

from abc import abstractmethod


class Diagnostic:
    """An abstract class representing a diagnostic."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None
    ):
        """Initialize the Diagnostic.

        Parameters
        ----------
        diagnostic_function
            The function that applies a diagnostic to the data.
        plotting_function
            The function that visualizes the results of the diagnostic.
        name : str
            The name of the diagnostic.
        description : str
            The description of the diagnostic.
        """
        self.name = name
        self._description = description
        self.diagnostic_function = diagnostic_function
        self.plotting_function = plotting_function

    @abstractmethod
    def apply(self, data):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        data
            The data to apply the diagnostic to. Data can be an xarray DataTree, Dataset or DataArray.

        Returns
        -------
        Results
            The data after applying the diagnostic either as a DataTree, Dataset, DataArray, Scalar, or a pandas DataFrame.
        """
        pass

    def plot(self, result, ax=None, **kwargs):
        """Plot the diagnostic.

        Parameters
        ----------
        result : xr.Dataset or xr.DataArray or DataTree
            The output of the diagnostic function.

        Returns
        -------
        ax : matplotlib.axis.Axis
            The axis (singular) of the plot.
        """
        return self.plotting_function(result, ax=ax, **kwargs)

    @property
    def description(self):
        """Return the description of the diagnostic a combination of the name, the type and the description and the docstring of the diagnostic and plot functions."""
        return f"{self.name} ({self.__class__.__name__})\n{self._description}\n Diagnostic function: {self.diagnostic_function.__name__}\n {self.diagnostic_function.__doc__}\n Visualization function: {self.plotting_function.__name__}\n {self.plotting_function.__doc__}"


class Model2Self(Diagnostic):
    """A class representing a diagnostic that compares a model to itself."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None
    ):
        """Initialize the Model2Self diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description)

    def apply(self, ds: xr.Dataset, **kwargs):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        ds : xr.Dataset
            The data to apply the diagnostic to.

        Returns
        -------
        xr.Dataset
            The data after applying the diagnostic.
        """
        return self.diagnostic_function(ds, **kwargs)


class Model2Ref(Diagnostic):
    """A class representing a diagnostic that compares a model to a reference."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None
    ):
        """Initialize the Model2Ref diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description)

    def apply(self, ds: xr.Dataset, ref: xr.Dataset, **kwargs):
        """Apply the diagnostic to the data. Only the common variables between the data and the reference are used.

        Parameters
        ----------
        ds : xr.Dataset
            The data to apply the diagnostic to.
        ref : xr.Dataset
            The reference data to compare the data to.

        Returns
        -------
        xr.Dataset
            The data after applying the diagnostic.
        """
        ds, ref = _select_common_vars(ds, ref)
        return self.diagnostic_function(ds, ref, **kwargs)


class Ensemble2Self(Diagnostic):
    """A class representing a diagnostic that compares an ensemble to itself."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None
    ):
        """Initialize the Ensemble2Self diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description)

    def apply(self, dt: DataTree, **kwargs):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        dt : DataTree
            The data to apply the diagnostic to.

        Returns
        -------
        DataTree or dict
            The data after applying the diagnostic as a DataTree or a dictionary of results with the tree nodes as keys.
        """
        return self.diagnostic_function(dt, **kwargs)

    def plot(self, result, facetted=True, **kwargs):
        """Plot the diagnostic.

        If axes are provided, the diagnostic is plotted facetted. If ax is provided, the diagnostic is plotted non-facetted. 
        If neither axes nor ax are provided, the diagnostic is plotted on the current axis and no facetting is applied.

        Parameters
        ----------
        result : DataTree
            The result of applying the ensemble diagnostic to a DataTree.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        if "ax" in kwargs and "axes" in kwargs:
            raise ValueError("Either ax or axes can be provided, not both.")
        elif "ax" not in kwargs and "axes" not in kwargs:
            ax = plt.gca()
            return self.plotting_function(result, ax=ax, **kwargs)
        else:
            return self.plotting_function(result, **kwargs)


    @classmethod
    def from_model2self(cls, model2self: Model2Self, facetted=True):
        """Create an Ensemble2Self diagnostic from a Model2Self diagnostic.

        Parameters
        ----------
        model2self : Model2Self
            The Model2Self diagnostic to convert.

        Returns
        -------
        Ensemble2Self
            The Ensemble2Self diagnostic.
        """

        def diagnostic_function(dt: DataTree, **kwargs):
            return dt.map_over_subtree(model2self.diagnostic_function, **kwargs)

        def plotting_function(
            dt: DataTree, variable=None, **kwargs
        ):
            if "axes" in kwargs:
                axes = kwargs.pop("axes")
                for ds, ax in zip(dt.leaves, axes.flatten()):
                    model2self.plot(ds[variable], ax=ax, **kwargs)
                    ax.set_title(ds.path.replace("/", " "))
                return axes
            if "ax" in kwargs:
                ax = kwargs.pop("ax")
                for ds in dt.leaves:
                    model2ref.plot(ds[variable], ax=ax, label=f'{ds.path.replace("/", " ")}', **kwargs)
                return ax

        return Ensemble2Self(
            diagnostic_function,
            plotting_function,
            model2self.name,
            model2self.description,
        )


class Ensemble2Ref(Diagnostic):
    """A class representing a diagnostic that compares an ensemble to a reference."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None
    ):
        """Initialize the Ensemble2Ref diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description)

    def apply(self, dt: DataTree, ref, **kwargs):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        dt : DataTree
            The data to apply the diagnostic to.
        ref : xr.DataSet or DataTree
            The reference data to compare the data to.

        Returns
        -------
        DataTree or dict
            The data after applying the diagnostic as a DataTree or a dictionary of results with the tree nodes as keys.
        """
        # TODO: Add some checks to make sure the reference is a DataTree or a Dataset and contain common variables with the data.
        return self.diagnostic_function(dt, ref, **kwargs)

    def plot(self, result, facetted=True, **kwargs):
        """Plot the diagnostic.

        If axes are provided, the diagnostic is plotted facetted. If ax is provided, the diagnostic is plotted non-facetted. 
        If neither axes nor ax are provided, the diagnostic is plotted on the current axis and no facetting is applied.

        Parameters
        ----------
        result : DataTree
            The result of applying the ensemble diagnostic to a DataTree.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        if "ax" in kwargs and "axes" in kwargs:
            raise ValueError("Either ax or axes can be provided, not both.")
        elif "ax" not in kwargs and "axes" not in kwargs:
            ax = plt.gca()
            return self.plotting_function(result, ax=ax, **kwargs)
        else:
            return self.plotting_function(result, **kwargs)


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
            if isinstance(ref, DataTree):
                ensemble_results = {}
                for data_node, ref_node in zip(dt.leaves, ref.leaves):
                    ds, ref = _select_common_vars(data_node.ds, ref_node.ds)
                    ensemble_results[data_node.path] = model2ref.diagnostic_function(
                        ds, ref, **kwargs
                    )
                return DataTree.from_dict(ensemble_results)
            else:
                return dt.map_over_subtree(
                    model2ref.diagnostic_function, ref=ref, **kwargs
                )

        def plotting_function(
            dt: DataTree, variable=None, **kwargs
        ):
            if "axes" in kwargs:
                axes = kwargs.pop("axes")
                for ds, ax in zip(dt.leaves, axes.flatten()):
                    model2ref.plot(ds[variable], ax=ax, **kwargs)
                    ax.set_title(ds.path.replace("/", " "))
                return axes
            if "ax" in kwargs:
                ax = kwargs.pop("ax")
                for ds in dt.leaves:
                    model2ref.plot(ds[variable], ax=ax, label=f'{ds.path.replace("/", " ")}', **kwargs)
                return ax

        return Ensemble2Ref(
            diagnostic_function,
            plotting_function,
            model2ref.name,
            model2ref.description,
        )



def _common_vars(ds1, ds2):
    """Return the common variables in two datasets."""
    return set(ds1.data_vars).intersection(set(ds2.data_vars))


def _select_common_vars(ds1, ds2):
    """Select the common variables in two datasets."""
    common_vars = _common_vars(ds1, ds2)
    return ds1[common_vars], ds2[common_vars]


# =============================================================================
# Pre-made diagnostics
# =============================================================================

from valenspy.diagnostic.functions import *
from valenspy.diagnostic.visualizations import *

# Model2Self diagnostics
DiurnalCycle = Model2Self(
    diurnal_cycle, plot_diurnal_cycle, "Diurnal Cycle", "The diurnal cycle of the data."
)
TimeSeriesSpatialMean = Model2Self(
    time_series_spatial_mean,
    plot_time_series,
    "Time Series Spatial Mean",
    "The time series of the spatial mean of the data.",
)
# Model2Ref diagnostics
SpatialBias = Model2Ref(
    spatial_bias,
    plot_spatial_bias,
    "Spatial Bias",
    "The spatial bias of the data compared to the reference.",
)
TemporalBias = Model2Ref(
    temporal_bias,
    plot_time_series,
    "Temporal Bias",
    "The temporal bias of the data compared to the reference.",
)
DiurnalCycleBias = Model2Ref(
    diurnal_cycle_bias,
    plot_diurnal_cycle,
    "Diurnal Cycle Bias",
    "The diurnal cycle bias of the data compared to the reference.",
)

# Ensemble2Ref diagnostics
MetricsRankings = Ensemble2Ref(
    calc_metrics_dt,
    plot_metric_ranking,
    "Metrics Rankings",
    "The rankings of ensemble members with respect to several metrics when compared to the reference."
)