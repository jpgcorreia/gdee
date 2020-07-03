"""
"""


class DataContainer(dict):
    def __setattr__(self, name, value):
        self[name] = value

    def __getattr__(self, name):
        if name in self:
            return self[name]
        raise AttributeError(name)
