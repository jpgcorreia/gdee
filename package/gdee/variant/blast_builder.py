"""
"""


class BlastBuilder:
    def __init__(self, parameters, database):
        self.parameters = parameters
        self.database = database
        self._task_list = list(range(20))

    def next_job(self):
        if self._task_list:
            return self._task_list.pop()

        return None

    def save_results(self, data):
        print("REMOVE ME! Pipeline: Saving results:", data)
