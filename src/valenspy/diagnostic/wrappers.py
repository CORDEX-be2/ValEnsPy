######################################
############## Wrappers ##############
######################################
from functools import wraps


def required_variables(variables):
    """
    A decorator that checks if the required variables are present in the dataset (and reference dataset if applicable) before applying the diagnostic.
    The required variables are specified as a list of strings. Only if all the required variables are present, the diagnostic is applied.
    Note that this is a minimum requirement; the dataset may contain other variables than the required ones.

    Parameters
    ----------
    variables : str or list of str
        The variable(s) required to apply the diagnostic.

    Examples
    --------
    >>> # The diagnostic function requires the variables 'tas' and 'pr' to be present in the dataset.
    >>> @required_variables(["tas", "pr"])
    >>> def my_diagnostic(ds: xr.Dataset):
    >>>    return ds.tas + ds.pr
    
    >>> # This also checks if the variables are present in both the data and the reference.
    >>> # An error is raised if the required variables are not present in the data or the reference.
    >>> @required_variables(["tas", "pr"])
    >>> def my_diagnostic(ds: xr.Dataset, ref: xr.Dataset):
    >>>    return ds.tas + ref.pr
    """

    def decorator(diagnostic_function):
        @wraps(diagnostic_function)
        def wrapper(ds, *args, **kwargs):
            required_vars = [variables] if isinstance(variables, str) else variables
            # Do the check for the ds
            if not all(var in ds.data_vars for var in required_vars):
                raise ValueError(
                    f"Variables {required_vars} are required to apply the diagnostic."
                )
            # Do the check for the reference if it is present, the reference is the second argument after the ds argument and should be a xr.Dataset.
            if len(args) > 0 and isinstance(args[0], xr.Dataset):
                ref = args[0]
                if not all(var in ref.data_vars for var in required_vars):
                    raise ValueError(
                        f"Variables {required_vars} are required to apply the diagnostic."
                    )
            return diagnostic_function(*args, **kwargs)

        return wrapper

    return decorator

def acceptable_variables(variables):
    """
    Decorator that checks if the dataset contains at least one of the acceptable variables before applying the diagnostic.
    The required variables are specified as a list of strings. If at least one of the acceptable variables is present the diagnostic is applied.
    If an acceptable variable is present all acceptable variables are passed to the diagnostic function - other variables will be dropped.

    Parameters
    ----------
    variables : str or list of str
        The variable(s) that are acceptable to apply the diagnostic.

    Examples
    -----
    >>> # The diagnostic function accepts the variables 'tas', 'tas_max' and 'tas_min'.
    >>> @acceptable_variables(["tas", "tas_max", "tas_min"])
    >>> def my_diagnostic(ds: xr.Dataset):
    >>>    #Function which is valid for tas, tas_max and tas_min
    >>>    return result
    """

    def decorator(diagnostic_function):
        @wraps(diagnostic_function)
        def wrapper(ds, *args, **kwargs):
            acceptable_vars = [variables] if isinstance(variables, str) else variables
            # Check if at least one of the acceptable variables is present in the dataset
            if not any(var in ds.data_vars for var in acceptable_vars):
                raise ValueError(
                    f"At least one of the variables {acceptable_vars} is required to apply the diagnostic."
                )
            # Drop all variables that are not in the acceptable variables
            ds = ds.drop_vars([var for var in ds.data_vars if var not in acceptable_vars])
            return diagnostic_function(ds, *args, **kwargs)

        return wrapper
    
    return decorator