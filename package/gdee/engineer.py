"""
"""


from gdee.pipeline import PipelineFactory
from gdee.platform import PlatformFactory
import os
import signal


__all__ = ["ProteinEngineering"]


class ProteinEngineering:
    def __init__(self, protein_name, database):
        self.work_dir = os.getcwd()
        self.protein_name = protein_name
        self.db_file = database
        self.pdb = None
        self.ligand_pdbqt = None
        self.platform = {"name": "simple", "local_cpu": 1}
        self.programs = {"mgltools": "mgltools", "vina": "vina", "vinardo": "smina"}
        self.variant = {"name": "mutation", "matrix": "blosum62", "selection": "", "conservative": True, "max_iterations": 1000, "msa": ""}
        self.model = {"name": "modeller", "optimize_radius": 0, "num_models": 5, "optimize_level": 0}
        self.evaluator = {"name": "vina", "exhaustiveness": 50}
        self._measurements = []
        self._pipeline = None
        self._terminate = False
        signal.signal(signal.SIGTERM, self.catch_signals)

    def add_measurement(self, metric, protein_sel, ligand_sel):
        self._measurements.append((metric, protein_sel, ligand_sel))

    def run(self):
        pipeline_factory = PipelineFactory()
        pipeline_factory.protein_name = self.protein_name
        pipeline_factory.programs = self.programs
        pipeline_factory.work_dir = self.work_dir
        pipeline_factory.pdb = self.pdb
        pipeline_factory.ligand_pdbqt = self.ligand_pdbqt
        pipeline_factory.db_file = self.db_file
        pipeline_factory.variant_parameters = self.variant
        pipeline_factory.model_parameters = self.model
        pipeline_factory.evaluator_parameters = self.evaluator
        pipeline_factory.measurements = self._measurements
        self.pipeline = pipeline_factory.make()

        platform_factory = PlatformFactory()
        platform_factory.parameters = self.platform
        platform_factory.pipeline = self.pipeline
        platform = platform_factory.make()

        platform.run()

        if self._terminate:
            raise RuntimeError("Processing interrupted by a system signal. Everything should be fine")

    def catch_signals(self, signal, frame):
        self._terminate = True
        self.pipeline.terminate()
        print("\nCaught termination signal. Finalizing all workers. This may take some time\n")
