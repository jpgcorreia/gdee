"""
"""


from ..database import Database
from .blast_builder import BlastBuilder
from .mutation_builder import MutationBuilder


__all__ = ["VariantBuilderFactory"]


class VariantBuilderFactory:
    def __init__(self):
        self.parameters = {}

    def make(self):
        if self.parameters["pdb_file"] is None:
            raise RuntimeError("Invalid template PDB file.")

        database = Database(self.parameters["db_name"])

        name = self.parameters["name"]
        if name == "blast":
            return BlastBuilder(self.parameters, database)

        elif name == "mutation":
            return MutationBuilder(self.parameters, database)

        else:
            raise RuntimeError("Variant builder '{}' does not exists.".format(name))
