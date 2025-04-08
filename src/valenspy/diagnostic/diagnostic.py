from datatree import DataTree
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

    def plot(self, result, title=None, **kwargs):
        """Plot the diagnostic. Single ax plots.

        Parameters
        ----------
        result : xr.Dataset or xr.DataArray or DataTree
            The output of the diagnostic function.

        Returns
        -------
        ax : matplotlib.axis.Axis
            The axis (singular) of the plot.
        """
        ax = self.plotting_function(result, **kwargs)
        if not title:
            title = self.name
        ax.set_title(title)
        return ax

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
        super().__init__(diagnostic_function, plotting_function, name, description)
        self.plot_type = plot_type

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

        Returns
        -------
        DataTree
            The data after applying the diagnostic.
        """
        return dt.map_over_subtree(self.apply, *args, **kwargs)

    def plot_dt(self, dt, *args, **kwargs):
        if self.plot_type == "single":
            return self.plot_dt_single(dt, *args, **kwargs)
        elif self.plot_type == "facetted":
            return self.plot_dt_facetted(dt, *args, **kwargs)

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

class Model2Self(DataSetDiagnostic):
    """A class representing a diagnostic that compares a model to itself."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="single"
    ):
        """Initialize the Model2Self diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type)

    def apply(self, ds: xr.Dataset, mask=None, **kwargs):
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
        if mask == "prudence":
            ds = add_prudence_regions(ds)
        return self.diagnostic_function(ds, **kwargs)
    
    def apply_dt(self, dt: DataTree, mask=None, **kwargs):
        """
        Apply the diagnostic to a DataTree.
        
        Parameters
        ----------
        dt : DataTree
            The DataTree to apply the diagnostic to.
            
        Returns
        -------
        DataTree
            The DataTree after applying the diagnostic.
        """
        if mask == "prudence":
            dt = dt.map_over_subtree(add_prudence_regions)
        return dt.map_over_subtree(self.diagnostic_function, **kwargs)


class Model2Ref(DataSetDiagnostic):
    """A class representing a diagnostic that compares a model to a reference."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="facetted"
    ):
        """Initialize the Model2Ref diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type)

    def apply(self, ds: xr.Dataset, ref: xr.Dataset, mask=None, **kwargs):
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
        if mask == "prudence":
            ds = add_prudence_regions(ds)
            ref = add_prudence_regions(ref)

        ds, ref = _select_common_vars(ds, ref)

        return self.diagnostic_function(ds, ref, **kwargs)

    def apply_dt(self, dt: DataTree, ref: xr.Dataset, mask=None, **kwargs):
        """
        Apply the diagnostic to a DataTree.
        
        Parameters
        ----------
        dt : DataTree
            The DataTree to apply the diagnostic to.
        ref : xr.Dataset
            The reference data to compare the data to.
            
        Returns
        -------
        DataTree
            The DataTree after applying the diagnostic.
        """
        if mask == "prudence":
            dt = dt.map_over_subtree(add_prudence_regions)
            ref = add_prudence_regions(ref)

        return dt.map_over_subtree(self.diagnostic_function, ref=ref, **kwargs)

class Ensemble2Self(Diagnostic):
    """A class representing a diagnostic that compares an ensemble to itself."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, iterative_plotting=False
    ):
        """Initialize the Ensemble2Self diagnostic."""
        self.iterative_plotting = iterative_plotting
        super().__init__(diagnostic_function, plotting_function, name, description)
        

    def apply(self, dt: DataTree, mask=None, **kwargs):
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
        if mask == "prudence":
            dt = dt.map_over_subtree(add_prudence_regions)

        return self.diagnostic_function(dt, **kwargs)

    def plot(self, result, variables=None, title=None, facetted=None, **kwargs):
        """Plot the diagnostic.

        If facetted multiple plots on different axes are created. If not facetted, the plots are created on the same axis.

        Parameters
        ----------
        result : DataTree
            The result of applying the ensemble diagnostic to a DataTree.

        Returns
        -------
        Figure
            The figure representing the diagnostic.
        """
        if not self.iterative_plotting:
            if facetted is not None:
                warnings.warn("facetted is ignored when using a non-iterative plotting function.")
            return self._plot_non_iterative(result, title=title, **kwargs)
        else:
            if variables is None:
                raise ValueError("variables must be provided when using an iterative plotting function. The variables can be a list of variables to plot or a single variable to plot.")
            return self._plot_iterative(result, title=title, variables=variables, facetted=facetted, **kwargs)

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
