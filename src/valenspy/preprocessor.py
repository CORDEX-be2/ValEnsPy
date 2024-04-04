from valenspy.preprocessing_tasks.task import PreprocessingTask
from datatree import DataTree

class Preprocessor:
    """A class that organizes the preprocessing of the input data to prepare it for a diagnostic."""

    def __init__(self):
        """Initialize the Preprocessor, it consists of an ordered list of preprocessing tasks."""

        self.preprocessing_tasks = []

    def add_preprocessing_task(self, task: PreprocessingTask):
        """Add a preprocessing task to the preprocessor."""
        self.preprocessing_tasks.append(task)

    def apply_preprocessing(self, dt: DataTree):
        """Apply all preprocessing tasks to the data in the DataTree."""
        for task in self.preprocessing_tasks:
            dt = dt.map_over_subtree(task.apply)
        return dt

    