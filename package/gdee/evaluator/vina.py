"""
"""


from .pdbqt import PDBQT
from path import Path
import subprocess
from tempfile import TemporaryDirectory, mkdtemp


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
        self.extra_arguments = []
        self.prepare_receptor = Path(parameters["mgltools"]) / "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"

    def run(self, job_data):
        job_dir = job_data["job_dir"]
        temp_dir = TemporaryDirectory(prefix="gdee_docking")
        temp_path = Path(temp_dir.name)
        ligand_file = self.parameters["ligand_pdbqt"]
        Path(ligand_file).copy(temp_path / "ligand.pdbqt")
        docking_data = []

        with temp_path:
            for idx, model_pdb in enumerate(job_data["models"]["pdbs"]):
                (job_dir / model_pdb).copy(temp_path / "model.pdb")
                self.run_docking(job_data)

                # Process and save results
                pdbqt = PDBQT("results.pdbqt")
                results_pdb = "docking_{:04d}.pdb".format(idx)
                pdbqt.write_pdb(job_dir / results_pdb)

                docking_data.append({
                    "ligand_file": ligand_file,
                    "method": self.name,
                    "pdb": results_pdb,
                    "energies": [model.energy for model in pdbqt]
                })

        job_data["evaluations"] = docking_data

        return job_data

    def run_docking(self, job_data):
        # Generate model's PDBQT
        command = [
            self.prepare_receptor,
            "-r", "model.pdb",
            "-o", "model.pdbqt",
        ]

        external_command(command, job_data["variant"].name)

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

        external_command(command, job_data["variant"].name)


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
