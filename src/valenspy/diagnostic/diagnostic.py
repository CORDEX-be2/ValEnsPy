from xarray import DataTree
import xarray as xr
import matplotlib.pyplot as plt
from valenspy.processing.mask import add_prudence_regions
from valenspy.diagnostic.plot_utils import _augment_kwargs
from valenspy._utilities import generate_parameters_doc
import numpy as np
import inspect
import textwrap

#Import get_axis from xarray
from xarray.plot.utils import get_axis

from abc import abstractmethod
import warnings

class Diagnostic():
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

        self.__signature__ = inspect.signature(self.diagnostic_function)
        self.__doc__ = self.description

    def __call__(self, data, *args, **kwargs):
        return self.apply(data, *args, **kwargs)
        

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

    def plot(self, result, **kwargs):
        """Plot the diagnostic.

        Parameters
        ----------
        result : xr.Dataset or xr.DataArray or DataTree
            The output of the diagnostic function.
        title : str
            The title of the plot.
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        ax : matplotlib.axis.Axis
            The axis (singular) of the plot.
        """
        return self.plotting_function(result, **kwargs)
    
    def plot_dt_single(self, dt, var, ax, label="name", colors=None, **kwargs):
        """
        Plot the diagnostic by iterating over the leaves of a DataTree.
        
        Parameters
        ----------
        dt : DataTree
            The DataTree to plot.
        var : str
            The variable to plot.
        ax : matplotlib.axis.Axis
            The axis to plot on.
        label : str
            The attribute of the DataTree nodes to use as a title for the plots.
        colors : dict or list
            The colors to use for the different leaves of the DataTree.
            Either a dictionary with the colors as values and the DataTree paths as keys or a list of colors.
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        ax : matplotlib.axis.Axis
            The axis of the plot.
        """
        if colors:
            if isinstance(colors, list):
                colors = {dt_leave.path: color for dt_leave, color in zip(dt.leaves, colors)}

        for dt_leave in dt.leaves:
            if label:
                kwargs["label"] = getattr(dt_leave, label)
            if colors:
                kwargs["color"] = colors[dt_leave.path]
            self.plot(dt_leave[var], ax=ax, **kwargs)

        return ax
        
    def plot_dt_facetted(self, dt, var, axes, label="name", shared_cbar=None, **kwargs):
        """
        Plot the diagnostic by iterating over the leaves of a DataTree.
        
        Parameters
        ----------
        dt : DataTree
            The DataTree to plot.
        var : str
            The variable to plot.
        axes : np.ndarray
            The axes to plot on.
        label : str
            The attribute of the DataTree nodes to use as a title for the plots.
        shared_cbar : str
            How to handle the vmin and vmax of the plot. Options are None, "min_max", "abs".
            If None, the vmin and vmax are not automatically set. Passing the vmin and vmax as kwargs will still result in shared colorbars. 
            If "min_max", the vmin and vmax are set respectively to the minimum and maximum over all the leaves of the DataTree. 
            If "abs", the vmin and vmax are set to the maximum of the absolute value of the minimum and maximum over all the leaves of the DataTree.
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        axes : np.ndarray
            The axes of the plot.
        """
        #Flatten the axes if needed
        #Add option if axes is not provided to create new axes

        if shared_cbar:
            max = np.max([ds[var].values for ds in dt.max().leaves])
            min = np.min([ds[var].values for ds in dt.min().leaves])
            if shared_cbar == "min_max":
                kwargs = _augment_kwargs({"vmin": min, "vmax": max}, **kwargs)
            elif shared_cbar == "abs":
                abs_max = np.max([np.abs(min), np.abs(max)])
                kwargs = _augment_kwargs({"vmin": -abs_max, "vmax": abs_max}, **kwargs)

        for ax, dt_leave in zip(axes, dt.leaves):
            if label:
                kwargs["title"] = getattr(dt_leave, label)
            self.plot(dt_leave[var], ax=ax, **kwargs)
        return axes

    @property
    def description(self):
        """Generate the docstring for the diagnostic."""
        name_no_spaces = self.name.replace(" ", "")
        title = f"{self.name} - {self.__class__.__name__}\n\n"
        description = f"{self._description}\n\n"
        params = generate_parameters_doc(self.diagnostic_function)
        see_also = f"See also\n--------\n:py:class:`{self.__class__.__name__}`, :func:`{self.diagnostic_function.__module__}.{self.diagnostic_function.__name__}`,:func:`{self.plotting_function.__module__}.{self.plotting_function.__name__}` : Plotting function\n\n"
        examples = f"Examples\n--------\n>>> from valenspy.diagnostic import {name_no_spaces}\n>>> result = {name_no_spaces}(ds)\n>>> {name_no_spaces}.plot(result)\n\n"
        docstring = f"{title}{description}{params}{see_also}{examples}"
        return textwrap.dedent(docstring)

class DataSetDiagnostic(Diagnostic):
    """A class representing a diagnostic that operates on the level of single datasets."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="single"
    ):
        """
        Initialize the DataSetDiagnostic.
        
        Parameters
        ----------
        plot_type : str
            The type of plot to create. Options are "single" or "facetted".
            If "single", plot_dt will plot all the leaves of the DataTree on the same axis.
            If "facetted", plot_dt will plot all the leaves of the DataTree on different axes.
        """
        if plot_type not in ["single", "facetted"]:
            raise ValueError("Invalid plot_type provided. Options are 'single' or 'facetted'.")
        self.plot_type = plot_type
        super().__init__(diagnostic_function, plotting_function, name, description)
        

    def __call__(self, data, *args, **kwargs):
        if isinstance(data, DataTree):
            return self.apply_dt(data, *args, **kwargs)
        else:
            return self.apply(data, *args, **kwargs)
        
    def apply_dt(self, dt: DataTree, *args, **kwargs):
        """Apply the diagnostic to a DataTree by iterating over the each dataset in the tree.

        Parameters
        ----------
        dt : DataTree
            The data to apply the diagnostic to.
        *args
            Positional arguments to pass to the diagnostic function.
        **kwargs
            Keyword arguments to pass to the diagnostic function.

        Returns
        -------
        DataTree
            The data after applying the diagnostic.
        """
        #Bug fix needed until https://github.com/pydata/xarray/issues/9693 is resolved
        def apply(ds, *args, **kwargs):
            if not ds:
                return ds
            return self.apply(ds, *args, **kwargs)
        return dt.map_over_datasets(apply, *args, **kwargs)

    def apply(self, ds: xr.Dataset, *args, **kwargs):
        """Apply the diagnostic to a single dataset.

        Parameters
        ----------
        ds : xr.Dataset
            The data to apply the diagnostic to.
        *args
            Positional arguments to pass to the diagnostic function.
        **kwargs
            Keyword arguments to pass to the diagnostic function.

        Returns
        -------
        xr.Dataset
            The data after applying the diagnostic.
        """
        return self.diagnostic_function(ds, *args, **kwargs)

    def plot(self, result, title=None, **kwargs):
        """Plot the diagnostic. Single ax plots.

        Parameters
        ----------
        result : xr.Dataset or xr.DataArray or DataTree
            The output of the diagnostic function.
        title : str
            The title of the plot.
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        ax : matplotlib.axis.Axis
            The axis (singular) of the plot.
        """
        ax = super().plot(result, **kwargs)
        if not title:
            title = self.name
        ax.set_title(title)
        return ax

    def plot_dt(self, dt, *args, **kwargs):
        if self.plot_type == "single":
            return self.plot_dt_single(dt, *args, **kwargs)
        elif self.plot_type == "facetted":
            return self.plot_dt_facetted(dt, *args, **kwargs)

class DataTreeDiagnostic(Diagnostic):
    """A class representing a diagnostic that operates on the level of DataTrees."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type=None
    ):
        """Initialize the DataTreeDiagnostic.
        Parameters
        ----------
        plot_type : str, optional
            The type of plotting function to use. Default is None, which means the plotting function will be used as is. 
            Options are "single" or "facetted".
            If "single", plot_dt will plot all the leaves of the DataTree on the same axis.
            If "facetted", plot_dt will plot all the leaves of the DataTree on different axes.

        """
        if plot_type not in [None, "single", "facetted"]:
            raise ValueError("Invalid plot_type provided. Options are None, 'single', or 'facetted'.")
        self.plot_type = plot_type
        super().__init__(diagnostic_function, plotting_function, name, description)
        
    def __call__(self, data, *args, **kwargs):
        if not isinstance(data, DataTree):
            raise ValueError("Data must be a DataTree.")
        return self.apply(data, *args, **kwargs)
    
    def apply(self, dt: DataTree, *args, **kwargs):
        """Apply the diagnostic to a DataTree.

        Parameters
        ----------
        dt : DataTree
            The data to apply the diagnostic to.
        *args
            Positional arguments to pass to the diagnostic function.
        **kwargs
            Keyword arguments to pass to the diagnostic function.

        Returns
        -------
        DataTree or dict
            The data after applying the diagnostic as a DataTree or a dictionary of results with the tree nodes as keys.
        """
        return self.diagnostic_function(dt, *args, **kwargs)

    def plot_dt(self, dt, *args, **kwargs):
        """Plot the diagnostic by iterating over the leaves of a DataTree.

        Parameters
        ----------
        dt : DataTree
            The DataTree to plot.
        *args
            Positional arguments to pass to the plotting function.
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        #Check if the dt is a DataTree and if not raise an error
        if not isinstance(dt, DataTree):
            raise ValueError("Data must be a DataTree. Use self.plot to plot non-DataTree data results.")
        if not self.plot_type:
            warnings.warn("No plot_type specified, using the default plotting function. It is recommended to use self.plot instead of self.plot_dt when no plot_type is specified.")
            return self.plotting_function(dt, *args, **kwargs)
        elif self.plot_type == "single":
            return self.plot_dt_single(dt, *args, **kwargs)
        elif self.plot_type == "facetted":
            return self.plot_dt_facetted(dt, *args, **kwargs)
        else:
            raise ValueError("Invalid plot_type specified. Options are 'single', 'facetted', or None.")

class Model2Self(DataSetDiagnostic):
    """A class representing a diagnostic that compares a model to itself."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="single"
    ):
        """Initialize the Model2Self diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type)


class Model2Ref(DataSetDiagnostic):
    """A class representing a diagnostic that compares a model to a reference."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="facetted"
    ):
        """Initialize the Model2Ref diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type)

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

        return super().apply(ds, ref, **kwargs)

class Ensemble2Self(DataTreeDiagnostic):
    """A class representing a diagnostic that compares an ensemble to itself."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type=None
    ):
        """Initialize the Ensemble2Self diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type)

class Ensemble2Ref(DataTreeDiagnostic):
    """A class representing a diagnostic that compares an ensemble to a reference."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type=None
    ):
        """Initialize the Ensemble2Ref diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type)

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

def _common_vars(ds1, ds2):
    """Return the common variables in two datasets."""
    return set(ds1.data_vars).intersection(set(ds2.data_vars))

def _select_common_vars(ds1, ds2):
    """Select the common variables in two datasets."""
    common_vars = _common_vars(ds1, ds2)
    return ds1[common_vars], ds2[common_vars]

def _initialize_multiaxis_plot(n, subplot_kws={}):
    """Initialize a multi-axis plot."""
    fig, axes = plt.subplots(
            nrows=n//2+1, ncols=2, figsize=(10, 5 * n), subplot_kw=subplot_kws
        )
    return fig, axes
