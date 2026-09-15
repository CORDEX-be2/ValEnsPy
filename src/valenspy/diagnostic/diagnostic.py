from xarray import DataTree
import xarray as xr
import matplotlib.pyplot as plt
from valenspy.processing.mask import add_prudence_regions, mask_to_reference_coverage, _mask_ds_to_reference_coverage
from valenspy.diagnostic.plot_utils import _augment_kwargs
from valenspy._utilities import generate_parameters_doc
import numpy as np
import inspect
import re
import textwrap
from pathlib import Path

#Import get_axis from xarray
from xarray.plot.utils import get_axis

from abc import abstractmethod
import warnings

class Diagnostic():
    """An abstract class representing a diagnostic."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, short_name=None
    ):
        """Initialize the Diagnostic.

        Parameters
        ----------
        diagnostic_function
            The function that applies a diagnostic to the data.
        plotting_function
            The function that visualizes the results of the diagnostic.
        name : str
            The name of the diagnostic. Used as-is in `description`'s generated docstring, and
            as the fallback for `id`/`title` when `short_name` isn't given.
        description : str
            The description of the diagnostic.
        short_name : str, optional
            A brief label to use for `id`/`title` instead of `name`, for a diagnostic whose
            full `name` reads too long for a filename/title (e.g. name "Ensemble quantiles of
            closest member spatial mean", short_name "Closest Member Quantiles"). `name`
            itself is unaffected - `description`'s generated docstring always uses the full
            `name`. Default None (fall back to `name` everywhere).
        """
        self.name = name
        self.short_name = short_name
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
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        ax : matplotlib.axis.Axis
            The axis (singular) of the plot.
        """
        return self.plotting_function(result, **kwargs)
    
    def plot_dt_single(self, dt, var, ax=None, label="name", colors=None, subplot_kw=None, **kwargs):
        """
        Plot the diagnostic by iterating over the leaves of a DataTree.

        Parameters
        ----------
        dt : DataTree
            The DataTree to plot.
        var : str
            The variable to plot.
        ax : matplotlib.axis.Axis, optional
            The axis to plot on. If None, a new one is created (using `subplot_kw` if given).
        label : str
            The attribute of the DataTree nodes to use as a title for the plots.
        colors : dict or list
            The colors to use for the different leaves of the DataTree.
            Either a dictionary with the colors as values and the DataTree paths as keys or a list of colors.
        subplot_kw : dict, optional
            Passed to `plt.subplots` when `ax` is None. Ignored if `ax` is given.
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        ax : matplotlib.axis.Axis
            The axis of the plot.
        """
        if ax is None:
            _, ax = plt.subplots(subplot_kw=subplot_kw or {})

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
        
    def plot_dt_facetted(self, dt, var, axes=None, label="name", shared_cbar=None, subplot_kw=None, **kwargs):
        """
        Plot the diagnostic by iterating over the leaves of a DataTree.

        Parameters
        ----------
        dt : DataTree
            The DataTree to plot.
        var : str
            The variable to plot.
        axes : np.ndarray, optional
            The axes to plot on, one per leaf of `dt`. If None, new ones are created (using
            `subplot_kw` if given) - one row, one column per leaf.
        label : str
            The attribute of the DataTree nodes to use as a title for the plots.
        shared_cbar : str
            How to handle the vmin and vmax of the plot. Options are None, "min_max", "abs".
            If None, the vmin and vmax are not automatically set. Passing the vmin and vmax as kwargs will still result in shared colorbars.
            If "min_max", the vmin and vmax are set respectively to the minimum and maximum over all the leaves of the DataTree.
            If "abs", the vmin and vmax are set to the maximum of the absolute value of the minimum and maximum over all the leaves of the DataTree.
        subplot_kw : dict, optional
            Passed to `plt.subplots` when `axes` is None. Ignored if `axes` is given.
        **kwargs
            Keyword arguments to pass to the plotting function.

        Returns
        -------
        axes : np.ndarray
            The axes of the plot.
        """
        #Check how to deal with shared_cbar (shared vmin and vmas - should this be named differently?) and should the cbar really be shared?

        if axes is None:
            _, axes = plt.subplots(1, len(list(dt.leaves)), subplot_kw=subplot_kw or {})
        axes = np.atleast_1d(axes).ravel()

        if shared_cbar:
            max = np.max([ds[var].values for ds in dt.max().leaves])
            min = np.min([ds[var].values for ds in dt.min().leaves])
            if shared_cbar == "min_max":
                kwargs = _augment_kwargs({"vmin": min, "vmax": max}, **kwargs)
            elif shared_cbar == "abs":
                abs_max = np.max([np.abs(min), np.abs(max)])
                kwargs = _augment_kwargs({"vmin": -abs_max, "vmax": abs_max}, **kwargs)

        for ax, dt_leave in zip(axes, dt.leaves):
            self.plot(dt_leave[var], ax=ax, **kwargs)
            if label:
                title = getattr(dt_leave, label)
                ax.set_title(title)
        
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

    @property
    def id(self):
        """Short filename-safe slug of `short_name` (or `name`), e.g. "spatial-bias". Used as
        the top-level folder in `filename`.
        """
        return re.sub(r"[^0-9a-zA-Z]+", "-", (self.short_name or self.name).strip()).strip("-").lower()

    @staticmethod
    def _format_param_value(value):
        """Stringify one parameter value; a list/tuple is "-"-joined."""
        if isinstance(value, (list, tuple)):
            return "-".join(Diagnostic._format_param_value(v) for v in value)
        return str(value)

    def _filename_sections(self, var=None, **kwargs):
        """[id, "var_<value>", "<key>_<value>", ...] for whichever of var/kwargs are given -
        one path segment per parameter, built by `filename` via `_build_path`. Overridden by
        subclasses (e.g. `_ReferenceComparisonNaming`) adding their own named segments.
        """
        sections = [self.id]
        if var is not None:
            sections.append(f"var_{self._format_param_value(var)}")
        for key, value in kwargs.items():
            if value is not None:
                sections.append(f"{key}_{self._format_param_value(value)}")
        return sections

    @staticmethod
    def _build_path(sections, ext):
        """[a, b, c] -> "a/b/c.ext" - every section but the last becomes a subfolder."""
        *dirs, last = sections
        return "/".join(dirs + [f"{last}.{ext}"])

    def _detail(self, **kwargs):
        """" (<key>=<value>, ...)" from whichever kwargs are given, or "" if none are."""
        given = {k: v for k, v in kwargs.items() if v is not None}
        if not given:
            return ""
        return " (" + ", ".join(f"{k}={self._format_param_value(v)}" for k, v in given.items()) + ")"

    def filename(self, ext="png", var=None, **kwargs):
        """Build a relative output path: one subfolder per given parameter (`var` first, then
        `**kwargs` in order), with the last one doubling as the file itself. `id` is always the
        top-level folder; parameters left as None produce no folder level.

        Parameters
        ----------
        ext : str, optional
            File extension, no leading dot. Default "png".
        var : str, optional
            The variable this output is for, e.g. "tas".
        **kwargs
            Any other parameter distinguishing this output, e.g. `region="belgium"`. A
            list/tuple value is "-"-joined. Subclasses (e.g. `_ReferenceComparisonNaming`) may
            add their own named folder level.

        Returns
        -------
        str
            e.g. "spatial-bias/var_tas/region_belgium.png".
        """
        return self._build_path(self._filename_sections(var=var, **kwargs), ext)

    @property
    def _title_name(self):
        """`short_name` if given, else `name`."""
        return self.short_name or self.name

    @staticmethod
    def _var_label(var, long_name):
        """`long_name` if given, else `var`, else None - the value `title` shows for the
        variable. `filename` has no equivalent; it always uses `var`, never `long_name`.
        """
        return long_name if long_name is not None else var

    def title(self, var=None, long_name=None, **kwargs):
        """Build a readable title: "<name> of <var>", plus any other given parameter as a
        trailing "(key=value, ...)". `long_name`, if given, replaces `var` in the sentence -
        valenspy has no variable-name lookup of its own, so this is only ever what the caller
        passes. Nothing beyond the diagnostic's own name is required.

        Parameters
        ----------
        var : str, optional
            The variable this output is for, e.g. "tas".
        long_name : str, optional
            A more readable label to show instead of `var`, e.g. "Near-Surface Air Temperature".
        **kwargs
            Any other parameter worth noting, e.g. `region="belgium"`.

        Returns
        -------
        str
            e.g. "Spatial Bias of Near-Surface Air Temperature (region=belgium)".
        """
        var_label = self._var_label(var, long_name)
        base = f"{self._title_name} of {self._format_param_value(var_label)}" if var_label is not None else self._title_name
        return base + self._detail(**kwargs)

    def _plotting_function_wants_var(self):
        """Whether `plotting_function` itself takes a `var` parameter (e.g.
        plot_reference_future_periods_grid(result, var, ...)) as opposed to expecting an
        already-`var`-selected DataArray (e.g. plot_map(da, ...)) - both conventions exist
        across valenspy's own plotting functions, so `run` checks this rather than guessing.
        """
        try:
            return "var" in inspect.signature(self.plotting_function).parameters
        except (TypeError, ValueError):
            return False

    def run(self, data, ref=None, var=None, out_dir=None, ext="png", save_result=None,
            compute_kwargs=None, plot_kwargs=None, filename_kwargs=None, title_kwargs=None):
        """Compute this diagnostic's result and plot (+ optionally save) it in one call - the
        all-in-one convenience form of calling the diagnostic directly and then `plot`/
        `plot_dt` separately. Purely additive: computing, saving the raw result, and plotting
        by hand remain fully supported and are exactly what `run` is built from - nothing here
        requires going through `run`. Just `compute, then delegate to render` - see `render`
        for everything after the compute step (plot/title/save), split out for a caller that
        wants to compute ONCE and plot several times without recomputing - e.g. a diagnostic
        whose result already covers every variable, rendered once per variable:

        >>> result = diagnostic(data, ref, **compute_kwargs)   # compute once
        >>> for var in ("tas", "pr"):
        ...     diagnostic.render(result, var=var, out_dir="figures")   # no recompute

        Chaining diagnostics needs no separate object either - call an earlier diagnostic
        directly for its raw result, then `run` (or `render`, if a var-independent stage was
        already computed once - see above) the last one on that result:

        >>> intermediate = diagnostic_a(data, ref)
        >>> result, ax = diagnostic_b.run(intermediate, out_dir="figures")

        Parameters
        ----------
        data
            Passed to the diagnostic - `self(data, ref, **compute_kwargs)` if `ref` is given,
            else `self(data, **compute_kwargs)`.
        ref : optional
            The reference data, for a diagnostic whose `apply` takes one (Model2Ref,
            Ensemble2Ref). Leave None for one that doesn't.
        save_result : str or Path, optional
            If given, the raw result is also saved here via `save` - `.to_netcdf` for an
            xarray object, `.to_csv` for a DataFrame (e.g. a table-shaped diagnostic's
            output).
        compute_kwargs : dict, optional
            Extra keyword arguments for the compute step.
        var, out_dir, ext, plot_kwargs, filename_kwargs, title_kwargs
            See `render` - forwarded to it unchanged.

        Returns
        -------
        result, ax
            The diagnostic's result, and the axis (or array of axes) it was plotted on.
        """
        compute_kwargs = compute_kwargs or {}
        result = self(data, ref, **compute_kwargs) if ref is not None else self(data, **compute_kwargs)

        if save_result is not None:
            self.save(result, save_result)

        ax = self.render(
            result, var=var, out_dir=out_dir, ext=ext,
            plot_kwargs=plot_kwargs, filename_kwargs=filename_kwargs, title_kwargs=title_kwargs,
        )
        return result, ax

    def render(self, result, var=None, out_dir=None, ext="png",
               plot_kwargs=None, filename_kwargs=None, title_kwargs=None):
        """Plot (+ optionally save) an ALREADY-COMPUTED result - the plot/title/save half of
        `run`, on its own, for a caller that computed `result` itself (directly, or via a
        previous `run`/`render` call) and wants to render it without recomputing - see `run`'s
        own docstring for the "compute once, render per variable" example this exists for.

        Parameters
        ----------
        result
            An already-computed diagnostic result - whatever `self(data, ref, ...)` (or
            `run`'s own first return value) produces.
        var : str, optional
            The variable to plot, and passed on to `filename`/`title` for output naming. For a
            DataTree result, forwarded to `plot_dt` (which requires it). For any other result,
            forwarded to `plot` only if `plotting_function` itself takes a `var` parameter;
            otherwise `result[var]` is selected first if `result` is a Dataset containing it.
        out_dir : str or Path, optional
            If given, the figure is saved to `out_dir / self.filename(ext=ext, var=var,
            **filename_kwargs)`, creating any missing parent folders.
        ext : str, optional
            Figure file extension. Default "png".
        plot_kwargs, filename_kwargs, title_kwargs : dict, optional
            Extra keyword arguments for the plot/filename/title steps respectively - e.g.
            `title_kwargs={"reference": "ERA5"}` for a Model2Ref/Ensemble2Ref diagnostic (see
            `_ReferenceComparisonNaming`).

        Returns
        -------
        ax
            The axis (or array of axes) `result` was plotted on.
        """
        plot_kwargs = dict(plot_kwargs or {})
        filename_kwargs = filename_kwargs or {}
        title_kwargs = title_kwargs or {}

        if isinstance(result, DataTree):
            ax = self.plot_dt(result, var=var, **plot_kwargs)
        elif var is not None and self._plotting_function_wants_var():
            ax = self.plot(result, var=var, **plot_kwargs)
        elif var is not None and isinstance(result, xr.Dataset) and var in result.data_vars:
            ax = self.plot(result[var], **plot_kwargs)
        else:
            ax = self.plot(result, **plot_kwargs)

        fig = np.atleast_1d(ax).ravel()[0].figure
        fig.suptitle(self.title(var=var, **title_kwargs))

        if out_dir is not None:
            out_file = Path(out_dir) / self.filename(ext=ext, var=var, **filename_kwargs)
            out_file.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_file)

        return ax

    @staticmethod
    def save(result, path):
        """Save a diagnostic result to `path` - `.to_netcdf()` for an xarray Dataset/
        DataArray/DataTree, `.to_csv()` for a DataFrame (e.g. a table-shaped diagnostic's
        output, one whose `plotting_function` isn't meant to be used at all - only its raw
        result matters). Creates any missing parent folders. The same helper `run`'s own
        `save_result=` uses internally, exposed directly for a caller that wants to save a
        result WITHOUT also plotting it - `run`/`render` always plot; this doesn't.
        """
        _save_result(result, path)

class DataSetDiagnostic(Diagnostic):
    """A class representing a diagnostic that operates on the level of single datasets."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="single", short_name=None
    ):
        """
        Initialize the DataSetDiagnostic.

        Parameters
        ----------
        plot_type : str
            The type of plot to create. Options are "single" or "facetted".
            If "single", plot_dt will plot all the leaves of the DataTree on the same axis.
            If "facetted", plot_dt will plot all the leaves of the DataTree on different axes.
        short_name : str, optional
            See Diagnostic.__init__.
        """
        if plot_type not in ["single", "facetted"]:
            raise ValueError("Invalid plot_type provided. Options are 'single' or 'facetted'.")
        self.plot_type = plot_type
        super().__init__(diagnostic_function, plotting_function, name, description, short_name)
        

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
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type=None, short_name=None
    ):
        """Initialize the DataTreeDiagnostic.
        Parameters
        ----------
        plot_type : str, optional
            The type of plotting function to use. Default is None, which means the plotting function will be used as is.
            Options are "single" or "facetted".
            If "single", plot_dt will plot all the leaves of the DataTree on the same axis.
            If "facetted", plot_dt will plot all the leaves of the DataTree on different axes.
        short_name : str, optional
            See Diagnostic.__init__.

        """
        if plot_type not in [None, "single", "facetted"]:
            raise ValueError("Invalid plot_type provided. Options are None, 'single', or 'facetted'.")
        self.plot_type = plot_type
        super().__init__(diagnostic_function, plotting_function, name, description, short_name)
        
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
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="single", short_name=None
    ):
        """Initialize the Model2Self diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type, short_name)

    def apply(self, ds: xr.Dataset, *args, mask_to_reference=None, **kwargs):
        """Apply the diagnostic to a single dataset, optionally masking it first to
        another dataset's own coverage.

        Parameters
        ----------
        mask_to_reference : xr.Dataset, optional
            If given, `ds` is restricted first to wherever `mask_to_reference` has a
            valid value anywhere along time (see `mask_to_reference_coverage`'s own
            docstring) - NOT variable-specific, since this diagnostic itself isn't
            (e.g. AnnualCycle, applied to every variable a dataset has at once): the
            restriction is to whichever variables `mask_to_reference` itself has, not
            one named variable, so this works whether `ds`/`mask_to_reference` carry
            one variable or several. Must be a Dataset here, not a path - there is no
            DataTree to resolve a path against when calling on a single Dataset
            directly; see `apply_dt` for that.
        """
        if mask_to_reference is not None:
            if isinstance(mask_to_reference, str):
                raise TypeError(
                    "mask_to_reference as a path is only valid when calling with a "
                    "DataTree (see apply_dt) - there is no tree to resolve it against here."
                )
            ds = _mask_ds_to_reference_coverage(ds, mask_to_reference, dim="time")
        return super().apply(ds, *args, **kwargs)

    def apply_dt(self, dt: DataTree, *args, mask_to_reference=None, **kwargs):
        """Apply the diagnostic to a DataTree, optionally masking every leaf first to
        another dataset's own coverage - see `apply`'s own docstring for what masking
        does and why it isn't variable-specific.

        Parameters
        ----------
        mask_to_reference : xr.Dataset or str, optional
            A Dataset (see `apply`), or a path naming one of `dt`'s OWN branches (e.g.
            "observations/CLIMATE_GRID"), resolved via `dt[mask_to_reference].ds`
            before masking - convenient when the reference is itself one of `dt`'s own
            branches, rather than an externally supplied Dataset.
        """
        if mask_to_reference is not None:
            if isinstance(mask_to_reference, str):
                mask_to_reference = dt[mask_to_reference].ds
            dt = mask_to_reference_coverage(dt, mask_to_reference)
        return super().apply_dt(dt, *args, **kwargs)


class _ReferenceComparisonNaming:
    """filename/title structure for diagnostics comparing data to a reference (Model2Ref,
    Ensemble2Ref): folds `reference` into the sentence as "<name> of <var> compared to
    <reference>". `reference` is a plain label the caller supplies for naming purposes only -
    `ref` itself (a Dataset/DataTree) has no name of its own.
    """

    def filename(self, ext="png", var=None, reference=None, **kwargs):
        """See Diagnostic.filename. `reference` becomes its own folder level, right after `var`."""
        ordered_kwargs = {"reference": reference, **kwargs}
        return self._build_path(self._filename_sections(var=var, **ordered_kwargs), ext)

    def title(self, var=None, long_name=None, reference=None, **kwargs):
        """See Diagnostic.title. `reference`, if given, extends the sentence as "... compared
        to <reference>" instead of appearing in the "(key=value, ...)" parenthetical.
        """
        var_label = self._var_label(var, long_name)
        base = f"{self._title_name} of {self._format_param_value(var_label)}" if var_label is not None else self._title_name
        if reference is not None:
            base = f"{base} compared to {self._format_param_value(reference)}"
        return base + self._detail(**kwargs)

class Model2Ref(_ReferenceComparisonNaming, DataSetDiagnostic):
    """A class representing a diagnostic that compares a model to a reference."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type="facetted", short_name=None
    ):
        """Initialize the Model2Ref diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type, short_name)

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
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type=None, short_name=None
    ):
        """Initialize the Ensemble2Self diagnostic."""
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type, short_name)

#: Default dataset attributes identifying an ensemble member for Ensemble2Ref's
#: ref-to-data matching (see match_ref_to_data) - matches the attrs intake-esm
#: attaches to a catalog-loaded dataset (InputManager's own catalog machinery),
#: not specific to any one project's diagnostics.
DEFAULT_IDENTITY_ATTRS = (
    "intake_esm_attrs:source_id", "intake_esm_attrs:driving_source_id", "intake_esm_attrs:driving_variant_label",
)


def _member_identity(ds, identity_attrs):
    return tuple(ds.attrs.get(attr) for attr in identity_attrs)


def _identity_map(dt, identity_attrs):
    """{identity: dataset} for every real leaf of `dt`, keyed by each leaf's own
    `identity_attrs` values - lossy if `dt` has more than one leaf sharing the same
    identity (last one wins), so only ever used on `ref` below, where a reference
    is expected to have at most one dataset per identity. `dt`'s own side of the
    match walks every leaf directly instead (see match_ref_to_data), specifically
    so that doesn't apply to it.
    """
    result = {}
    for path, node in dt.subtree_with_keys:
        if not path or node.children or node.dataset is None:
            continue
        result[_member_identity(node.dataset, identity_attrs)] = node.dataset
    return result


def match_ref_to_data(dt: DataTree, ref: DataTree, identity_attrs=DEFAULT_IDENTITY_ATTRS) -> DataTree:
    """Re-key `ref` so its member paths exactly match `dt`'s, pairing members by
    identity (`identity_attrs`) rather than by tree path.

    `dt` and `ref` frequently don't share tree paths even when they cover the same
    ensemble members - e.g. a reference/historical branch and several future
    scenario branches for the same model don't share a path segment for the
    scenario/experiment, since that's exactly the thing distinguishing them. A
    member present under several different branches of `dt` (e.g. the same model
    run under more than one future scenario) legitimately needs the same `ref`
    member duplicated under each of `dt`'s paths, not merged into one entry -
    they're independent comparisons against the same baseline, which is why this
    walks every leaf of `dt` individually rather than going through `dt`'s own
    identity map (see _identity_map's docstring on why that would silently drop
    all but one of several same-identity paths). Any `dt` member whose identity
    isn't found in `ref` is simply absent from the result.
    """
    ref_by_identity = _identity_map(ref, identity_attrs)
    matched = {}
    for path, node in dt.subtree_with_keys:
        if not path or node.children or node.dataset is None:
            continue
        identity = _member_identity(node.dataset, identity_attrs)
        if identity in ref_by_identity:
            matched[path] = ref_by_identity[identity]
    return DataTree.from_dict(matched)


class Ensemble2Ref(_ReferenceComparisonNaming, DataTreeDiagnostic):
    """A class representing a diagnostic that compares an ensemble to a reference."""

    def __init__(
        self, diagnostic_function, plotting_function, name=None, description=None, plot_type=None,
        identity_attrs=DEFAULT_IDENTITY_ATTRS, short_name=None,
    ):
        """Initialize the Ensemble2Ref diagnostic.

        Parameters
        ----------
        identity_attrs : tuple of str, optional
            Dataset attribute names identifying an ensemble member, used by `apply`
            to pair `ref`'s members to `dt`'s when `ref` is itself a DataTree (see
            `match_ref_to_data`). Default matches intake-esm-cataloged data's own
            attrs - override for data cataloged differently.
        short_name : str, optional
            See Diagnostic.__init__.
        """
        super().__init__(diagnostic_function, plotting_function, name, description, plot_type, short_name)
        self.identity_attrs = identity_attrs

    def apply(self, dt: DataTree, ref, **kwargs):
        """Apply the diagnostic to the data.

        Parameters
        ----------
        dt : DataTree
            The data to apply the diagnostic to.
        ref : xr.DataSet or DataTree
            The reference data to compare the data to. If a DataTree, its members
            are first re-paired to `dt`'s own members by identity (see
            `match_ref_to_data`) - `ref` need not already share `dt`'s tree paths.

        Returns
        -------
        DataTree or dict
            The data after applying the diagnostic as a DataTree or a dictionary of results with the tree nodes as keys.
        """
        # TODO: Add some checks to make sure the reference is a DataTree or a Dataset and contain common variables with the data.
        if isinstance(ref, DataTree):
            ref = match_ref_to_data(dt, ref, identity_attrs=self.identity_attrs)
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

def _save_result(result, path):
    """Save a diagnostic result to `path` - `to_netcdf` for an xarray Dataset/DataArray/
    DataTree, `to_csv` for a DataFrame (e.g. a table-shaped diagnostic's output). Used by
    `Diagnostic.run`'s `save_result`.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(result, (xr.Dataset, xr.DataArray, DataTree)):
        result.to_netcdf(path)
    elif hasattr(result, "to_csv"):
        result.to_csv(path)
    else:
        raise TypeError(f"Don't know how to save a result of type {type(result).__name__}.")
