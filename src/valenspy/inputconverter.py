from pathlib import Path
from typing import Callable, Union


class InputConverter:
    """A class that converts input files to netCDF files in CF convention."""

    def __init__(self, converter: Callable):
        """Initialize the InputProcessor with a converter function.

        Parameters
        ----------
        converter : function
            The function to convert the input file(s) to a netCDF file in CF convention.
            The function should take a string or list of strings as input and return a netCDF file.
        """
        self.converter = converter

    def convert_input(
        self, paths: Union[Path, list[Path]], metadata_info=None
    ) -> Union[Path, list[Path]]:
        """Convert the paths file to CF convention.

        Parameters
        ----------
        paths : Path or list(Path)
            The input file or list of input files to convert.
        """
        paths = self.converter(paths, metadata_info)

        return paths


# Idea is to extend the shared functionality here (with subclasses if required) while the inputconvertor_functions are model specific.

# Needed:
#  - Some helper functions to extend input to glob arguments, str arguments etc.
#  - Extend input so that already loaded datasets can also be input
#  - CF Checker functionality

# To be discussed -> Do we expect inputconvertor function to work at the file level? Or should they be able to manage a collection of files?
# If file level the InputConverter can seperately provide them to the function and handel concetanating/joining these.
