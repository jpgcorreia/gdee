"""
"""


from gdee.pipeline import PipelineFactory
from gdee.platform import PlatformFactory
import os


__all__ = ["ProteinEngineering"]


class ProteinEngineering:
    def __init__(self, protein_name, database):
        self.work_dir = os.getcwd()
        self.protein_name = protein_name
        self.db_file = database
        self.pdb = None
        self.platform = "simple"
        self.variant = {"name": "mutation", "matrix": "blosum62", "selection": "", "conservative": True, "max_iterations": 1000}
        self.model = {"name": "modeller"}
        self.evaluator = {"name": "vina", "exhaustiveness": 50}

    def run(self):
        pipeline_factory = PipelineFactory()
        pipeline_factory.protein_name = self.protein_name
        pipeline_factory.work_dir = self.work_dir
        pipeline_factory.pdb = self.pdb
        pipeline_factory.db_file = self.db_file
        pipeline_factory.variant_parameters = self.variant
        pipeline_factory.model_parameters = self.model
        pipeline_factory.evaluator_parameters = self.evaluator
        pipeline = pipeline_factory.make()

        platform_factory = PlatformFactory()
        platform_factory.name = self.platform
        platform_factory.pipeline = pipeline
        platform = platform_factory.make()

        platform.run()
