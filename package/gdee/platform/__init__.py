"""
"""


from .simple_platform import SimplePlatform
from .mpi_platform import MPIPlatform


__all__ = ["PlatformFactory"]


class PlatformFactory:
    def __init__(self):
        self.name = None
        self.pipeline = None

    def make(self):
        if self.name == "simple":
            return SimplePlatform(self.pipeline)

        elif self.name == "mpi":
            return MPIPlatform(self.pipeline)

        else:
            raise RuntimeError("Platform '{}' does not exists.".format(self.name))
