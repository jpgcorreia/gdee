"""
"""


from .modeller_builder import ModellerBuilder


__all__ = ["ModelBuilderFactory"]


class ModelBuilderFactory:
    def __init__(self):
        self.parameters = {}

    def make(self):
        if self.parameters["pdb_file"] is None:
            raise RuntimeError("Invalid template PDB file.")

        name = self.parameters["name"]
        if name == "modeller":
            return ModellerBuilder(self.parameters)

        # if name == "rosetta":
        #     return RosettaBuilder(self.parameters)

        else:
            raise RuntimeError("Modeler '{}' does not exists.".format(name))
