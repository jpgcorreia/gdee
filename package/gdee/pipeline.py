"""
"""


from gdee.variant import VariantBuilderFactory
from gdee.modeling import ModelBuilderFactory
from gdee.evaluator import EvaluatorFactory
from path import Path
import os
import signal


__all__ = ["PipelineFactory"]


class PipelineFactory:
    def __init__(self):
        self.protein_name = None
        self.programs = {}
        self.work_dir = None
        self.pdb = None
        self.ligand_pdbqt = None
        self.db_file = None
        self.variant_parameters = {}
        self.model_parameters = {}
        self.evaluator_parameters = {}

    def make(self):
        self.pdb = Path(self.pdb).abspath()
        pipeline = Pipeline()
        pipeline.work_dir = Path(self.work_dir).abspath() / "files"
        pipeline.work_dir.makedirs_p()
        signal.signal(signal.SIGTERM, pipeline.catch_signals)

        variant_factory = VariantBuilderFactory()
        self.variant_parameters["db_file"] = self.db_file
        self.variant_parameters["pdb_file"] = self.pdb
        self.variant_parameters["protein_name"] = self.protein_name
        variant_factory.parameters = self.variant_parameters
        pipeline.variant_builder = variant_factory.make()

        model_factory = ModelBuilderFactory()
        self.model_parameters["pdb_file"] = self.pdb
        model_factory.parameters = self.model_parameters
        pipeline.add_task(model_factory.make())

        evaluator_factory = EvaluatorFactory()
        self.evaluator_parameters["ligand_pdbqt"] = self.ligand_pdbqt
        self.evaluator_parameters.update(self.programs)
        evaluator_factory.parameters = self.evaluator_parameters
        pipeline.add_task(evaluator_factory.make())

        return pipeline


class Pipeline:
    def __init__(self):
        self.database = None
        self.work_dir = Path()
        self._variant_builder = None
        self.task_list = []
        self._terminate = False

    @property
    def variant_builder(self):
        return self._variant_builder

    @variant_builder.setter
    def variant_builder(self, obj):
        self._variant_builder = obj

    def add_task(self, task):
        self.task_list.append(task)

    def next_job(self, size):
        if self._terminate:
            return None

        job_list = []
        for i in range(size):
            job = self.variant_builder.next_job()
            if job is None:
                break

            job_list.append(job)

        return job_list

    def run_pipeline(self, job_data):
        job_dir = self.work_dir / job_data["variant_dir"]
        job_dir.makedirs_p()
        job_data["job_dir"] = job_dir
        job_data["fatal_error"] = False

        for step in self.task_list:
            job_data = step.run(job_data)
            if job_data["fatal_error"]:
                return job_data

        return job_data

    def save_results(self, data):
        for result in data:
            if result["fatal_error"]:
                self._variant_builder.remove_variant(result)
                print("Error while processing variant: {}".format(result["variant"].name))

            self._variant_builder.save_results(result)

    def catch_signals(self, signal, frame):
        self._terminate = True
        print("\nCaught termination signal. Finalizing all workers. This may take some time\n")
