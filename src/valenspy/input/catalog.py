import pandas as pd
import intake_esm

from valenspy._utilities._formatting import parse_string_to_time_period

class ValenspyEsmDatastore(intake_esm.esm_datastore):
    """
    Subclass of intake_esm.ESMDataStore for ValEnsPy.
    
    This extends the ESMDataStore class with a adittional search functionality for time based searching using the time_period column.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def search(self, require_all_on: str | list[str] | None = None, **query):
        """
        Search for entries in the catalog.

        Standard search function of the intake_esm.esm_datastore class extended with time based searching based on the time_period column.

        Parameters
        ----------
        require_all_on : str, optional
             A dataframe column or a list of dataframe columns across
            which all entries must satisfy the query criteria.
            If None, return entries that fulfill any of the criteria specified
            in the query, by default None.
        **query:
            keyword arguments corresponding to user's query to execute against the dataframe.

        See Also
        --------
        :func:`intake_esm.esm_datastore.search`
        """

        time_query = query.pop("time_period", None)

        if len(query) == 0:
            cat = self.__class__({"esmcat": self.esmcat.dict(), "df": self.esmcat._df})
        else:
            cat = super().search(require_all_on=require_all_on, **query)

        if time_query:
            df = cat.esmcat.df
            if isinstance(time_query, str):
                start, end = parse_string_to_time_period(time_query)
            elif isinstance(time_query, list):
                start, _ = parse_string_to_time_period(time_query[0])
                _, end = parse_string_to_time_period(time_query[1])
            else:
                raise ValueError("time_period should be a string or a list of strings")

            #Filter keeping files which cover a period which overlaps with the time_period
            df = df[(pd.to_datetime(df["time_period_start"]) <= end) & (start <= pd.to_datetime(df["time_period_end"]))]
            cat.esmcat._df = df

        return cat

    def to_datatree(self, levels: list[str] = None, **kwargs):

        #Hack to avoid deepcopy the esmcat object which breaks the search function? Dont know why
        #Probably when updating the catalog the df is not really updated
        #Maybe overwrite the deepcopy function?
        if levels:
            self.esmcat.aggregation_control.groupby_attrs, old_agg = levels, self.esmcat.aggregation_control.groupby_attrs

        dt = super().to_datatree(**kwargs)

        if levels:
            self.esmcat.aggregation_control.groupby_attrs = old_agg

        return dt
        