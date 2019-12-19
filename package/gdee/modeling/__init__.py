"""
"""


from .modeller_builder import ModellerBuilder


__all__ = ["ModelBuilderFactory"]


class ModelBuilderFactory:
    def __init__(self):
        self.name = None
        self.work_dir = None
        self.pdb_file = None
        self.parameters = {}

    def make(self):
        if self.pdb_file is None:
            raise RuntimeError("Invalid template PDB file.")

        params = ModelBuilderParameters(self.work_dir, self.pdb_file, self.parameters)

        if self.name == "modeller":
            return ModellerBuilder(params)

        # elif self.name == "rosetta":
        #     return RosettaBuilder(params)

        else:
            raise RuntimeError("Modeler '{}' does not exists.".format(self.name))


class ModelBuilderParameters:
    def __init__(self, work_dir, pdb_file, extra_params):
        self.work_dir = work_dir
        self.pdb_file = pdb_file
        self.extra_params = extra_params
