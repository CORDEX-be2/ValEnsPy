from functools import wraps

def _merge_kwargs(def_kwargs, kwargs):
    return {**def_kwargs, **kwargs}

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