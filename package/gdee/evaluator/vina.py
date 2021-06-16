"""
"""


import numpy as np
from path import Path
import MDAnalysis as mda
import subprocess
from tempfile import TemporaryDirectory
from .pdbqt import PDBQT
from gdee.misc import DataContainer
import warnings

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


def external_command(arguments, name):
    proc = subprocess.run(
        arguments,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if proc.returncode:
        raise RuntimeError("Error processing job '{}':\n{}\n".format(name, proc.stderr.decode("UTF-8")))


class BaseVina:
    def __init__(self, parameters):
        self.parameters = parameters
        self.name = ""
        self.program = ""
        self.ligand = parameters["ligand"]
        self.extra_arguments = []
        self.prepare_receptor = Path(parameters["mgltools"]) / "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"

    def run(self, job_data):
        job_dir = job_data.job_dir
        temp_dir = TemporaryDirectory(prefix="gdee_docking")
        temp_path = Path(temp_dir.name)
        Path(self.ligand.filename).copy(temp_path / "ligand.pdbqt")

        with temp_path:
            for idx, model in enumerate(job_data.modeling.models):
                if "evals" not in model:
                    model.evals = {}

                if model.rejected:
                    continue

                protein = mda.Universe(str(job_dir / model.pdb))
                pos = protein.atoms.positions
                size = np.array(self.parameters["box_size"], np.float32) + 1
                center = np.array(self.parameters["box_center"], np.float32)
                upper = center + size
                lower = center - size
                atoms = np.all((pos <= upper) & (pos >= lower), axis=1)
                smaller = protein.atoms[atoms].residues.atoms
                smaller.write(str(temp_path / "model.pdb"))

                try:
                    self.run_docking(job_data)
                    pdbqt = PDBQT("results.pdbqt")

                except Exception as error:
                    print(error)

                else:
                    # Process and save results
                    if pdbqt.size():
                        results_pdb = "docking_{}_{:04d}.pdb".format(self.ligand.name,
                                                                     idx)
                        pdbqt.write_pdb(job_dir / results_pdb)

                        docking = DataContainer()
                        docking.ligand_name = self.ligand.name
                        docking.ligand_file = self.ligand.filename
                        docking.method = self.name
                        docking.pdb = results_pdb
                        docking.energies = [model.energy for model in pdbqt]

                        model.evals[self.ligand.name] = docking

        return job_data

    def run_docking(self, job_data):
        # Generate model's PDBQT
        command = [
            self.prepare_receptor,
            "-r", "model.pdb",
            "-o", "model.pdbqt",
            "-A", "checkhydrogens",
        ]

        external_command(command, job_data.variant.name)

        # Run docking
        box_center = "--center_x {:.2f} --center_y {:.2f} --center_z {:.2f}".format(*self.parameters["box_center"])
        box_size = "--size_x {:.2f} --size_y {:.2f} --size_z {:.2f}".format(*self.parameters["box_size"])

        command = [
            self.program,
            "--exhaustiveness", self.parameters["exhaustiveness"],
            "--cpu", "1",
            "--num_modes", "500",   # Exaggerated
            "--energy_range", "30", # numbers
            "--receptor", "model.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--out", "results.pdbqt"
        ] + box_center.split(" ") + box_size.split(" ")
        command += self.extra_arguments
        command = list(map(str, command))

        external_command(command, job_data.variant.name)


class VinaDocking(BaseVina):
    def __init__(self, parameters, *args, **kwargs):
        super().__init__(parameters, *args, **kwargs)
        self.name = "vina"
        self.program = parameters["vina"]


class VinardoDocking(BaseVina):
    def __init__(self, parameters, *args, **kwargs):
        super().__init__(parameters, *args, **kwargs)
        self.name = "vinardo"
        self.program = parameters["vinardo"]
        self.extra_arguments = ["--scoring", "vinardo"]
