"""
"""


from .pipeline import PipelineFactory
from .platform import PlatformFactory
import os


__all__ = ["ProteinEngineering"]


class ProteinEngineering:
    def __init__(self, database):
        self.work_dir = os.getcwd()
        self.db_file = database

        self.platform = "simple"

        self.variant_name = "mutation"
        self.variant_parameters = {}

        self.model_name = "modeller"
        self.model_pdb = None
        self.model_parameters = {}

        self.evaluator_name = "vina"
        self.evaluator_parameters = {"exhaustiveness": 50}
        self.ligand_pdb = None
        self.ligand_box = None
        self.ligand_box_center = None

    def run(self):
        pipeline_factory = PipelineFactory()
        pipeline_factory.work_dir = self.work_dir
        pipeline_factory.db_file = self.db_file
        pipeline_factory.variant_name = self.variant_name
        pipeline_factory.variant_parameters = self.variant_parameters
        pipeline_factory.model_name = self.model_name
        pipeline_factory.model_pdb = self.model_pdb
        pipeline_factory.model_parameters = self.model_parameters
        pipeline_factory.evaluator_name = self.evaluator_name
        pipeline_factory.evaluator_parameters = self.evaluator_parameters
        pipeline_factory.ligand_pdb = self.ligand_pdb
        pipeline_factory.ligand_box = self.ligand_box
        pipeline_factory.ligand_box_center = self.ligand_box_center
        pipeline = pipeline_factory.make()

        platform_factory = PlatformFactory()
        platform_factory.name = self.platform
        platform_factory.pipeline = pipeline
        platform = platform_factory.make()

        platform.run()
