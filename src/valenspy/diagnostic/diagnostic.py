from xarray import DataTree
import xarray as xr
import matplotlib.pyplot as plt
from valenspy.processing.mask import add_prudence_regions
from valenspy.diagnostic.plot_utils import _augment_kwargs
from valenspy._utilities import generate_parameters_doc
import numpy as np
import inspect
import re
import textwrap

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
        #Check how to deal with shared_cbar (shared vmin and vmas - should this be named differently?) and should the cbar really be shared?

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
        """A short, filename-safe identifier for this diagnostic, derived from `short_name` if
        given, else `name` (e.g. "Spatial Bias" -> "spatial-bias"). Words are joined with "-".
        Used as the top-level folder of a generated output path - see `filename`.
        """
        return re.sub(r"[^0-9a-zA-Z]+", "-", (self.short_name or self.name).strip()).strip("-").lower()

    @staticmethod
    def _format_param_value(value):
        """Render one parameter value as a short, filename/title-safe token. A list/tuple
        (e.g. a set of quantiles or future periods) is joined with "-" rather than rendered
        with Python's own repr punctuation.
        """
        if isinstance(value, (list, tuple)):
            return "-".join(Diagnostic._format_param_value(v) for v in value)
        return str(value)

    def _filename_sections(self, var=None, **kwargs):
        """[id, "var_<value>" (if var given), "<key>_<value>" for each other kwarg given] -
        the section list shared by `filename` and overridden by subclasses (e.g.
        `_ReferenceComparisonNaming`) that want extra named sections between `var` and the
        rest. Each section becomes its own path segment in `filename` (see `_build_path`), so
        unlike an earlier version of this method, sections don't need to be mutually
        distinguishable by punctuation alone - "_" simply joins a key to its value (e.g.
        "future_periods_gwl2-gwl3"), and "-" is reserved purely for joining a list/tuple
        value's own items (see `_format_param_value`).
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
        """Turn a `_filename_sections`-style list into a relative path: every section but the
        last becomes a subfolder, and the last gets the extension - e.g.
        ["spatial-bias", "var_tas", "reference_ERA5"] -> "spatial-bias/var_tas/reference_ERA5.png".
        Real, nested folders rather than one long, densely-punctuated filename - each
        parameter is legible on its own, and outputs sharing an earlier parameter (e.g. every
        var="tas" output) naturally land in the same folder.
        """
        *dirs, last = sections
        return "/".join(dirs + [f"{last}.{ext}"])

    def _detail(self, **kwargs):
        """" (<key>=<value>, ...)" parenthetical from whichever kwargs are given (not None),
        or "" if none are - the optional trailing part shared by `title` and its subclass
        overrides.
        """
        given = {k: v for k, v in kwargs.items() if v is not None}
        if not given:
            return ""
        return " (" + ", ".join(f"{k}={self._format_param_value(v)}" for k, v in given.items()) + ")"

    def filename(self, ext="png", var=None, **kwargs):
        """Build a short, readable relative path for one output of this diagnostic.

        Each parameter becomes its own subfolder, in the order given, with the last one
        doubling as the file itself (extension appended) - e.g. `filename(var="tas",
        reference="ERA5")` on SpatialBias returns "spatial-bias/var_tas/reference_ERA5.png".
        Real folders keep each parameter legible without resorting to a punctuation-heavy
        single filename, and outputs sharing a parameter (e.g. every `var="tas"` output, from
        any diagnostic call) land under the same folder. `id` (see `Diagnostic.id`) is always
        the top-level folder; `filename()` alone (no var/kwargs) is just `f"{self.id}.{ext}"`
        with no subfolders at all.

        Only parameters actually given (not None) produce a folder level.

        Parameters
        ----------
        ext : str, optional
            The file extension, without a leading dot. Default "png".
        var : str, optional
            The variable this output is for, e.g. "tas". Virtually always given in practice
            (almost every diagnostic call is per-variable) - kept as its own parameter (rather
            than just another kwarg) so it always sits in the same position, right after `id`.
        **kwargs
            Any other parameter distinguishing this particular output from another produced
            by the same diagnostic, e.g. `region="belgium"`. A list/tuple value (e.g.
            `future_periods=["ssp245", "ssp585"]`) is joined with "-". Subclasses with a
            specific calling convention (e.g. Model2Ref's `reference`) may add their own named
            folder level - see the subclass's own `filename` if overridden.

        Returns
        -------
        str
            e.g. "spatial-bias/var_tas/region_belgium.png".
        """
        return self._build_path(self._filename_sections(var=var, **kwargs), ext)

    @property
    def _title_name(self):
        """`short_name` if given, else `name` - the name used to build `title` (and, via
        `id`, `filename`). `description`'s generated docstring always uses the full `name`
        regardless.
        """
        return self.short_name or self.name

    @staticmethod
    def _var_label(var, long_name):
        """The string to show for the variable in a title: `long_name` if given, else plain
        `var`, else None. `filename` has no equivalent - it always uses `var` itself (the
        short CF code), never `long_name`, to keep filenames compact.
        """
        return long_name if long_name is not None else var

    def title(self, var=None, long_name=None, **kwargs):
        """Build a readable title for one output of this diagnostic.

        The variable, if given, is folded directly into the sentence as "<name> of <var>"
        rather than shown as "var=<value>" - it's virtually always given (almost every
        diagnostic call is per-variable) and reads far more naturally inline than as a
        parameter. Every other parameter is optional: if given, it's appended as a lightweight
        "(key=value, ...)" parenthetical, just enough to disambiguate one output from another
        sharing the same name/var without cluttering the headline - nothing beyond the
        diagnostic's own name is ever required. Subclasses with a specific calling convention
        may fold a particular parameter into the sentence too, the same way this does for the
        variable - see e.g. `_ReferenceComparisonNaming`'s `reference`.

        Parameters
        ----------
        var : str, optional
            The variable this output is for, e.g. "tas" - shown in the sentence unless
            `long_name` is also given.
        long_name : str, optional
            A more readable label for the variable, e.g. "Near-Surface Air Temperature" -
            valenspy has no variable-name lookup of its own, so this is only ever what the
            caller passes; give it explicitly wherever a nicer label is wanted than the raw
            `var` code. When given, it's shown in the sentence in place of `var` - `var` itself
            is not otherwise referenced by `title` (unlike `filename`, which always uses `var`,
            never `long_name`, to keep filenames short).
        **kwargs
            Any other parameter worth noting, e.g. `region="belgium"`.

        Returns
        -------
        str
            e.g. "Spatial Bias of Near-Surface Air Temperature (region=belgium)", or plain
            "Spatial Bias" if neither the variable nor any kwarg is given. Uses `short_name` in
            place of `name` if one was set on this diagnostic (see `Diagnostic.__init__`).
        """
        var_label = self._var_label(var, long_name)
        base = f"{self._title_name} of {self._format_param_value(var_label)}" if var_label is not None else self._title_name
        return base + self._detail(**kwargs)

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


class _ReferenceComparisonNaming:
    """Shared filename/title structure for diagnostics comparing data to a reference
    (Model2Ref, Ensemble2Ref) - folds an optional `reference` kwarg into the same "<name> of
    <var>" sentence Diagnostic.title builds for `var`, as "<name> of <var> compared to
    <reference>", matching this apply(data, ref) calling convention. `reference` is a plain
    label (e.g. a dataset name) supplied by the caller for naming purposes only - `ref` itself
    is a Dataset/DataTree and has no name of its own to fall back on.
    """

    def filename(self, ext="png", var=None, reference=None, **kwargs):
        """See Diagnostic.filename. `reference`, if given, becomes its own "reference_<value>"
        folder level, positioned right after `var`'s.
        """
        ordered_kwargs = {"reference": reference, **kwargs}
        return self._build_path(self._filename_sections(var=var, **ordered_kwargs), ext)

    def title(self, var=None, long_name=None, reference=None, **kwargs):
        """See Diagnostic.title. `long_name`, if given, is shown in place of `var` in the
        sentence, same as Diagnostic.title. `reference`, if given, extends the sentence as
        "... compared to <reference>" rather than appearing in the "(key=value, ...)"
        parenthetical.
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
