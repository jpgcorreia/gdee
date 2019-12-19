"""
"""


from .database import Database
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
        print("Making a pipeline")
        pipeline = Pipeline()

        print("Making a database")
        pipeline.database = Database(self.db_name)

        print("Making a variant")
        variant_factory = VariantBuilderFactory()
        variant_factory.name = self.variant_name
        variant_factory.parameters = self.variant_parameters
        pipeline.add_initial_task(variant_factory.make())

        print("Making a model builder")
        model_factory = ModelBuilderFactory()
        model_factory.name = self.model_name
        model_factory.work_dir = self.work_dir
        model_factory.pdb_file = self.model_pdb
        model_factory.parameters = self.model_parameters
        pipeline.add_task(model_factory.make())

        print("Making an evaluator")
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
        self.initial_task_list = []
        self.task_list = []
        self._job_list = []

    def add_initial_task(self, task):
        print("Added task to initial list")
        self.initial_task_list.append(task)

    def add_task(self, task):
        print("Added task to list")
        self.task_list.append(task)

    def available(self):
        return bool(len(self._job_list))

    def next_job(self):
        return self._job_list.pop()

    def run_initial(self):
        print("Pipeline: Running initial")
        self._job_list += [1, 2, 3, 4]

    def run_pipeline(self, job_data):
        print("Pipeline: Running job:", job_data)
        return -job_data

    def save_results(self, results):
        print("Pipeline: Saving results:", results)
