"""
"""


from .measurer import Measurer
from .distance import EuclideanDistance


__all__ = ["MeasurerFactory"]


class MeasurerFactory:
    def __init__(self):
        self.measurements = []
        self.metrics = {
            EuclideanDistance.name(): EuclideanDistance,
        }

    def make(self):
        measurer = Measurer()
        for name, prot_sel, lig_sel in self.measurements:
            if name in self.metrics:
                measurer.add(self.metrics[name](), prot_sel, lig_sel)

            else:
                raise RuntimeError("Measurement category '{}' does not exists.".format(name))
        
        return measurer
