from functools import wraps

def _merge_kwargs(def_kwargs, kwargs):
    return {**def_kwargs, **kwargs}

def set_wrapped_suptitle(fig, title, **kwargs):
    """Set `title` as `fig`'s suptitle, wrapped to actually fit the figure's own width -
    unlike `textwrap.wrap`, which needs a character-count guess that goes stale the moment
    figsize/fontsize/n_cols change independently of each other.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    title : str
    **kwargs
        Passed to `fig.suptitle`, e.g. `fontsize=11`.

    Returns
    -------
    matplotlib.text.Text
    """
    suptitle = fig.suptitle(title, **kwargs)
    renderer = fig.canvas.get_renderer()
    fig_width = fig.get_window_extent(renderer).width

    lines, line = [], ""
    for word in title.split():
        candidate = f"{line} {word}".strip()
        suptitle.set_text(candidate)
        if line and suptitle.get_window_extent(renderer).width > fig_width:
            lines.append(line)
            line = word
        else:
            line = candidate
    lines.append(line)
    suptitle.set_text("\n".join(lines))
    return suptitle

def cbar_scale(vmin, vmax, kind):
    """{"vmin": ..., "vmax": ...} for a shared colorbar - the vmin/vmax -> plot kwargs
    conversion behind every `shared_cbar` option in this package.

    Parameters
    ----------
    vmin, vmax : float
        The data's own range.
    kind : str
        "min_max": `vmin`/`vmax` as given. "abs": a symmetric range around 0, sized to
        whichever of `vmin`/`vmax` has the larger absolute value.

    Returns
    -------
    dict
        {"vmin": ..., "vmax": ...}
    """
    if kind == "min_max":
        return {"vmin": vmin, "vmax": vmax}
    elif kind == "abs":
        abs_max = max(abs(vmin), abs(vmax))
        return {"vmin": -abs_max, "vmax": abs_max}
    raise ValueError(f"Invalid cbar kind {kind!r}. Options are 'min_max' or 'abs'.")

def _augment_kwargs(def_kwargs, **kwargs):
    """
    Augment the user provided keyword arguments with the default plot keyword arguments, subplot keyword arguments and colorbar keyword arguments.

    Parameters
    ----------
    def_kwargs : dict
        Default plot keyword arguments for the plotting function. 
        subplot_kws and cbar_kwargs can also be set and will also be augmented to the user provided subplot_kws and cbar_kwargs.
    kwargs : dict
        User provided keyword arguments.

    Returns
    -------
    dict
        Augmented keyword arguments.
    """

    if 'subplot_kws' in def_kwargs:
        subplot_kws = _merge_kwargs(def_kwargs.pop('subplot_kws'), kwargs.pop('subplot_kws', {}))
        def_kwargs['subplot_kws'] = subplot_kws
    
    if 'cbar_kwargs' in def_kwargs:
        cbar_kwargs = _merge_kwargs(def_kwargs.pop('cbar_kwargs'), kwargs.pop('cbar_kwargs', {}))
        def_kwargs['cbar_kwargs'] = cbar_kwargs
    
    return _merge_kwargs(def_kwargs, kwargs)

######################################
############## Wrappers ##############
######################################

def default_plot_kwargs(kwargs):
    """
    Decorator to set the default keyword arguments for the plotting function. User will override and/or be augmented with the default keyword arguments.
    subplot_kws and cbar_kwargs can also be set as default keyword arguments for the plotting function.

    Parameters
    ----------
    kwargs : dict
        Default keyword arguments for the plotting function. Can also include subplot_kws and cbar_kwargs as dictionarys in the kwargs dictionary.
    
    Examples
    --------
    The following example sets the default colorbar orientation to horizontal for the plotting function. 
    
    >>> @plot_kwarg_defaults({'cbar_kwargs': {'orientation': 'horizontal'}})
    ... def plot_function(*args, **kwargs):
    ...     pass

    If unspecified by the user, the colorbar orientation will be horizontal.
    If the user specifies the colorbar orientation, it will override the default orientation.
    If the user passes cbar_kwargs={'label': 'Label'}, the default orientation will still be horizontal and the label will be 'Label'.
    """
    
    def decorator(plotting_function):
        """Decortor function to set the default keyword arguments for the plotting function."""

        @wraps(plotting_function)
        def wrapper(*args, **kwargs):
            return plotting_function(*args, **_augment_kwargs(def_kwargs=kwargs, **kwargs))

        return wrapper

    return decorator