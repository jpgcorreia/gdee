"""
"""


from .vina import VinaDocking, VinardoDocking


__all__ = ["EvaluatorFactory"]


class EvaluatorFactory:
    def __init__(self):
        self.name = None
        self.work_dir = None
        self.pdb_file = None
        self.box = None
        self.box_center = None
        self.parameters = {}

    def make(self):
        if self.pdb_file is None:
            raise RuntimeError("Invalid ligand PDB file.")

        if self.box is None or len(self.box) != 3:
            raise RuntimeError("Invalid evaluator seraching box sizes.")

        if self.box_center is None or len(self.box_center) != 3:
            raise RuntimeError("Invalid searching box center.")

        params = EvaluatorParameters(
            self.work_dir,
            self.pdb_file,
            self.box,
            self.box_center,
            self.parameters
        )

        if self.name == "vina":
            return VinaDocking(params)

        elif self.name == "vinardo":
            return VinardoDocking(params)

        else:
            raise RuntimeError("Evaluator '{}' does not exists.".format(self.name))


class EvaluatorParameters:
    def __init__(self, work_dir, pdb_file, box, box_center, extra_params):
        self.work_dir = work_dir
        self.pdb_file = pdb_file
        self.sizes = box
        self.center = box_center
        self.extra_params = extra_params
