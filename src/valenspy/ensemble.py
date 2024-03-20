class Ensemble:
    """A class representing an ensemble which consists of multiple ensemble members."""

    def __init__(self, ensemble_members=[]):
        """Initialize an ensemble

        Parameters
        ----------
        ensemble_members : list of Ensmember objects
            The ensemble members that make up the ensemble
        """

        self.ensemble_members = ensemble_members
    