"""
"""


from gdee.variant import VariantBuilderFactory
from gdee.modeling import ModelBuilderFactory
from gdee.evaluator import EvaluatorFactory
import os


__all__ = ["PipelineFactory"]


class PipelineFactory:
    def __init__(self):
        self.protein_name = None
        self.work_dir = None
        self.pdb = None
        self.db_file = None
        self.variant_parameters = {}
        self.model_parameters = {}
        self.evaluator_parameters = {}

    def make(self):
        pipeline = Pipeline()

        variant_factory = VariantBuilderFactory()
        self.variant_parameters["db_file"] = self.db_file
        self.variant_parameters["pdb_file"] = self.pdb
        self.variant_parameters["protein_name"] = self.protein_name
        variant_factory.parameters = self.variant_parameters
        pipeline.variant_builder = variant_factory.make()

        model_factory = ModelBuilderFactory()
        self.model_parameters["work_dir"] = self.work_dir
        self.model_parameters["pdb_file"] = self.pdb
        model_factory.parameters = self.model_parameters
        pipeline.add_task(model_factory.make())

        evaluator_factory = EvaluatorFactory()
        self.evaluator_parameters["work_dir"] = self.work_dir
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
        return job_data

    def save_results(self, results):
        self._variant_builder.save_results(results)
