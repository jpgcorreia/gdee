"""
"""


class BaseVina:
    def __init__(self, parameters):
        self.parameters = parameters

    def run(self, job_data):
        return job_data


class VinaDocking(BaseVina):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class VinardoDocking(BaseVina):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
