"""
"""


from ..database import Database
from .blast_builder import BlastBuilder
from .mutation_builder import MutationBuilder


__all__ = ["VariantBuilderFactory"]


class VariantBuilderFactory:
    def __init__(self):
        self.name = None
        self.db_name = ""
        self.pdb_file = None
        self.parameters = {}

    def make(self):
        database = Database(self.db_name)
        if "pdb" not in self.parameters:
            self.parameters["pdb"] = self.pdb_file

        if self.name == "blast":
            return BlastBuilder(self.parameters, database)

        elif self.name == "mutation":
            return MutationBuilder(self.parameters, database)

        else:
            raise RuntimeError("Variant builder '{}' does not exists.".format(self.name))
