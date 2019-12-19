"""
"""


from .blast_builder import BlastBuilder
from .mutation_builder import MutationBuilder


__all__ = ["VariantBuilderFactory"]


class VariantBuilderFactory:
    def __init__(self):
        self.name = None
        self.parameters = {}

    def make(self):
        print("Choosing a variant builder")
        if self.name == "blast":
            return BlastBuilder(self.parameters)

        elif self.name == "mutation":
            return MutationBuilder(self.parameters)

        else:
            raise RuntimeError("Variant builder '{}' does not exists.".format(self.name))
