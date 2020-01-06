"""
"""


class BlastBuilder:
    def __init__(self, parameters):
        self.parameters = parameters
        self._task_list = list(range(20))
        print("Blast builder created")

    def next_job(self):
        if self._task_list:
            return self._task_list.pop()

        return None

