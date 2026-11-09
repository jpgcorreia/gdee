"""AutoDock Vina/Vinardo molecular docking implementation."""

import os
import time
import socket
import numpy as np
from path import Path
import MDAnalysis as mda
import subprocess
import re
from tempfile import TemporaryDirectory
from .pdbqt import PDBQT
from gdee.misc import DataContainer
import warnings
import fcntl
from contextlib import contextmanager

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


def external_command(arguments, name, env=None):
    """Execute external command and check for errors.

    Args:
        arguments: List of command arguments
        name: Job name for error reporting
        env: Environment variables for the subprocess

    Returns:
        subprocess.CompletedProcess: Command result

    Raises:
        RuntimeError: If command fails
    """
    proc = subprocess.run(
        arguments,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )

    if proc.returncode:
        raise RuntimeError("Error processing job '{}':\n{}\n".format(name, proc.stderr.decode("UTF-8")))

@contextmanager
def _gpu_slot(n_gpus, poll_s=0.2, lock_dir="/tmp", lock_prefix="gdee_vina_gpu"):
    if n_gpus <= 0:
        yield 0
        return

    host = socket.gethostname()
    fds = []
    try:
        for gid in range(int(n_gpus)):
            p = os.path.join(lock_dir, f"{lock_prefix}_{host}_{gid}.lock")
            fd = open(p, "a+")  # don't truncate
            fds.append((gid, fd, p))

        while True:
            for gid, fd, _p in fds:
                try:
                    fcntl.flock(fd.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    try:
                        yield gid
                    finally:
                        fcntl.flock(fd.fileno(), fcntl.LOCK_UN)
                    return
                except BlockingIOError:
                    continue
            time.sleep(poll_s)
    finally:
        for _gid, fd, _p in fds:
            try:
                fd.close()
            except Exception:
                pass

class BaseVina:
    """Base class for Vina-based docking engines.

    Handles protein preparation, search box validation, and result processing.
    """
    def __init__(self, parameters):
        """Initialize Vina docking base.

        Args:
            parameters: Docking configuration dictionary
        """
        self.parameters = parameters
        self.name = ""
        self.program = ""
        self.ligand = parameters["ligand"]
        self.extra_arguments = []
        self.prepare_receptor = Path(parameters["mgltools"]) / "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        self.atom_type = parameters['atom_type']

    def run(self, job_data):
        """Execute docking for all models in job.

        Args:
            job_data: Job data with models to dock

        Returns:
            DataContainer: Job data with docking results
        """
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
                if not len(smaller):
                    raise RuntimeError("No protein atoms inside the docking search box")
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
        """Prepare receptor PDBQT from PDB.
           Execute docking command.

        Args:
            job_data: Job data
        """
        # Generate model's PDBQT
        command = [
            self.prepare_receptor,
            "-r", "model.pdb",
            "-o", "model.pdbqt",
            "-A", "checkhydrogens",
            "-U", "nphs_lps_waters"
        ]

        external_command(command, job_data.variant.name)

        # Run atom type patch
        if self.atom_type:
            self.atom_type_patch('model.pdbqt')


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


    def atom_type_patch(self, pdbqt):
        """
        Patch PDBQT file to update atom types based on provided mapping.
        """
        new_pdbqt = 'model_patched.pdbqt'

        with open(pdbqt, 'r') as f, open(new_pdbqt, 'w') as nf:
            for line in f:
                write_line = True
                if line.startswith('ATOM'):
                    # "{name}:{chain}:{id}"
                    atom = "{}:{}:{}".format(line[12:16].strip(), line[21:22], line[22:26].strip())
                    if atom in self.atom_type:
                        if self.atom_type[atom] == "r":
                            write_line = False
                        else:
                            line = line[:77] + "{}\n".format(self.atom_type[atom])
                if write_line:
                    nf.write(line)

        os.remove(pdbqt)
        os.rename(new_pdbqt, pdbqt)


class VinaDocking(BaseVina):
    """AutoDock Vina docking implementation."""
    def __init__(self, parameters, *args, **kwargs):
        """Initialize Vina docking.

        Args:
            parameters: Docking configuration
            *args: Additional positional arguments
            **kwargs: Additional keyword arguments
        """
        super().__init__(parameters, *args, **kwargs)
        self.name = "vina"
        self.program = parameters["vina"]


class VinardoDocking(BaseVina):
    """Vinardo docking implementation.

    Uses Vinardo scoring function via Smina.
    """
    def __init__(self, parameters, *args, **kwargs):
        """Initialize Vinardo docking.

        Args:
            parameters: Docking configuration
            *args: Additional positional arguments
            **kwargs: Additional keyword arguments
        """
        super().__init__(parameters, *args, **kwargs)
        self.name = "vinardo"
        self.program = parameters["vinardo"]
        self.extra_arguments = ["--scoring", "vinardo"]

    def score_pose(self, receptor_pdb, ligand_pdbqt):
        """Score an existing receptor + ligand pose with Vinardo.

        Args:
            receptor_pdb: path-like to receptor PDB file
            ligand_pdbqt: path-like to ligand PDBQT file

        Returns:
            float: Vinardo score (kcal/mol) for the pose (first/only pose)
        """
        # Create temporary working directory and reuse run_docking machinery,
        # but run smina in score-only mode. Keep atom_type handling and prepare_receptor.
        temp_dir = TemporaryDirectory(prefix="gdee_vinardo_score")
        temp_path = Path(temp_dir.name)
        try:
            # Copy inputs to temp and convert PDB ligands for smina.
            Path(receptor_pdb).copy(temp_path / "model.pdb")
            ligand_path = Path(ligand_pdbqt)
            if ligand_path.suffix.lower() == ".pdbqt":
                ligand_path.copy(temp_path / "ligand.pdbqt")
            else:
                prepare_ligand = Path(self.prepare_receptor).parent / "prepare_ligand4.py"
                if not prepare_ligand.exists():
                    raise RuntimeError("A PDBQT ligand or MGLTools prepare_ligand4.py is required")
                subprocess.run(
                    [
                        str(prepare_ligand),
                        "-l", str(ligand_path.absolute()),
                        "-o", "ligand.pdbqt",
                    ],
                    check=True,
                    cwd=str(temp_path),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )

            # prepare receptor if needed (reuse current logic)
            command = [
                str(self.program),
                "--score_only",
                "--receptor", "model.pdbqt",
                "--ligand", "ligand.pdbqt"
            ]
            # if this class uses extra arguments for vinardo scoring, keep them
            command += self.extra_arguments
            command = list(map(str, command))

            # if prepare_receptor script exists, run it first (same as docking flow)
            if Path(self.prepare_receptor).exists():
                prep_cmd = [
                    str(self.prepare_receptor),
                    "-r", "model.pdb",
                    "-o", "model.pdbqt",
                    "-A", "checkhydrogens",
                    "-U", "nphs_lps_waters"
                ]
                subprocess.run(
                    list(map(str, prep_cmd)),
                    check=True,
                    cwd=str(temp_path),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )

            proc = subprocess.run(
                command,
                check=False,
                cwd=str(temp_path),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if proc.returncode != 0:
                raise RuntimeError("Vinardo scoring failed: {}".format(proc.stderr.decode("utf-8")))

            # parse stdout or output file depending on smina behavior; smina prints score lines to stdout
            out = proc.stdout.decode("utf-8").strip().splitlines()
            # find numeric score in output (simple heuristic)
            for line in reversed(out):
                match = re.search(r"(?:Affinity|score)\s*:?\s*(-?\d+(?:\.\d+)?)", line, re.IGNORECASE)
                if match:
                    return float(match.group(1))

            # fallback: raise
            raise RuntimeError("Could not parse vinardo score from smina output")

        finally:
            temp_dir.cleanup()

class VinaGPUDocking(BaseVina):
    """AutoDock Vina GPU docking implementation."""
    def __init__(self, parameters, *args, **kwargs):
        """Initialize Vina GPU docking.

        Args:
            parameters: Docking configuration
            *args: Additional positional arguments
            **kwargs: Additional keyword arguments
        """
        super().__init__(parameters, *args, **kwargs)
        self.name = "vina_gpu"
        self.program = parameters["vina_gpu"]

    def run_docking(self, job_data):

        command = [
            self.prepare_receptor,
            "-r", "model.pdb",
            "-o", "model.pdbqt",
            "-A", "checkhydrogens",
            "-U", "nphs_lps_waters"
        ]
        external_command(command, job_data.variant.name)

        if self.atom_type:
            self.atom_type_patch('model.pdbqt')

        box_center = "--center_x {:.2f} --center_y {:.2f} --center_z {:.2f}".format(*self.parameters["box_center"])
        box_size = "--size_x {:.2f} --size_y {:.2f} --size_z {:.2f}".format(*self.parameters["box_size"])

        command = [
            self.program,
            "--thread", self.parameters["gpu_thread"],
            "--opencl_binary_path", self.parameters["opencl_binary_path"],
            "--num_modes", "500",
            "--energy_range", "30",
            "--receptor", "model.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--out", "results.pdbqt"
        ] + box_center.split(" ") + box_size.split(" ")
        command += self.extra_arguments
        command = list(map(str, command))

        with _gpu_slot(self.parameters["n_gpu"]) as gpu_id:
            cmd = list(command)
            env = dict(os.environ)
            env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

            external_command(cmd, job_data.variant.name, env=env)
