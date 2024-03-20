from valenspy.modeldata import Modeldata
from valenspy.ensmember import Ensmember
class InputProcessor:
    """Class to read an input file and convert it to CF convention."""
    def __init__(self, converter):
        """Initialize the InputProcessor with a converter function.
        
        Parameters
        ----------
        converter : function
            The function to convert the input file(s) to a netCDF file in CF convention. 
            The function should take a string or list of strings as input and return a netCDF file.
        """
        self.converter = converter

    def convert_input(self, input):
        """Convert the input file to CF convention.
        
        Parameters
        ----------
        input : str or list(str)
            The input file or list of input files to convert.
        """
        input = self.converter(input)

        if isinstance(input, str):
            return Modeldata(file_location=input)
        elif isinstance(input, list):
            return [Modeldata(file_location=file) for file in input]

#Maybe classes are not needed here rather objects of the InputProcessor class with different converter functions.
#To be discussed.
    
class FA_CORDEXTRACTOR(InputProcessor):
    """Class to convert FA files to netCDF CF convention files using the CORDEXTRACTOR."""
    from valenspy.inputprocessor_functions import fa_cordextractor

    def __init__(self):
        """Initialize the FA_CORDEXTRACTOR InputProcessor with the fa_cordextractor function."""
        super().__init__(self.fa_cordextractor)

class FA_PyFa(InputProcessor):
    """Class to convert FA files to netCDF CF convention files using PyFa."""
    from valenspy.inputprocessor_functions import fa_pyfa

    def __init__(self):
        """Initialize the FA_PyFa InputProcessor with the fa_pyfa function."""
        super().__init__(self.fa_pyfa)

