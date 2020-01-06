"""
"""


from .variant import VariantBuilderFactory
from .modeling import ModelBuilderFactory
from .evaluator import EvaluatorFactory
import os


__all__ = ["PipelineFactory"]


class PipelineFactory:
    def __init__(self):
        self.work_dir = None
        self.db_name = None
        self.variant_name = None
        self.variant_parameters = None
        self.model_name = None
        self.model_pdb = None
        self.model_parameters = None
        self.evaluator_name = None
        self.evaluator_parameters = None
        self.ligand_pdb = None
        self.ligand_box = None
        self.ligand_box_center = None

    def make(self):
        pipeline = Pipeline()

        variant_factory = VariantBuilderFactory()
        variant_factory.name = self.variant_name
        variant_factory.parameters = self.variant_parameters
        variant_factory.db_name = self.db_name
        pipeline.variant_builder = variant_factory.make()

        model_factory = ModelBuilderFactory()
        model_factory.name = self.model_name
        model_factory.work_dir = self.work_dir
        model_factory.pdb_file = self.model_pdb
        model_factory.parameters = self.model_parameters
        pipeline.add_task(model_factory.make())

        evaluator_factory = EvaluatorFactory()
        evaluator_factory.name = self.evaluator_name
        evaluator_factory.work_dir = self.work_dir
        evaluator_factory.pdb_file = self.ligand_pdb
        evaluator_factory.box = self.ligand_box
        evaluator_factory.box_center = self.ligand_box_center
        evaluator_factory.parameters = self.evaluator_parameters
        pipeline.add_task(evaluator_factory.make())

        return pipeline


class Pipeline:
    def __init__(self):
        self.database = None
        self._variant_builder = None
        self.task_list = []

    @property
    def variant_builder(self):
        return self._variant_builder

    @variant_builder.setter
    def variant_builder(self, obj):
        self._variant_builder = obj

    def add_task(self, task):
        self.task_list.append(task)

    def next_job(self):
        return self.variant_builder.next_job()

    def run_pipeline(self, job_data):
        return -job_data

    def save_results(self, results):
        self._variant_builder.save_results(results)
