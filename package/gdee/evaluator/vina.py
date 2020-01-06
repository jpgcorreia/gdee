"""
"""


class BaseVina:
    def __init__(self, parameters):
        self.parameters = parameters


class VinaDocking(BaseVina):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class VinardoDocking(BaseVina):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
