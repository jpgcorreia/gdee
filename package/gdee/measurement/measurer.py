"""
"""


from path import Path
import MDAnalysis as mda


class Task:
    def __init__(self, metric, prot_text, lig_text):
        self.name = "{}|{}|{}".format(metric.name(), prot_text, lig_text)
        self.metric = metric
        self.prot_text = prot_text
        self.lig_text = lig_text
        self.data = []

    def clear(self):
        self.data = []

    def set_system(self, protein, ligand):
        self.prot_sel = protein.select_atoms(self.prot_text)
        self.lig_sel = ligand.select_atoms(self.lig_text)

        if not len(self.prot_sel) or not len(self.lig_sel):
            return False

        return True

    def compute(self):
        self.data.append(self.metric.compute(self.prot_sel.positions, self.lig_sel.positions))


class Measurer:
    def __init__(self):
        # Use a dict to avoid repeated measurements
        self.task_list = {}

    def add(self, metric, prot_sel, lig_sel):
        task = Task(metric, prot_sel, lig_sel)
        self.task_list[task.name] = task

    def run(self, job_data):
        if "models" not in job_data or "pdbs" not in job_data["models"] \
           or not job_data["models"]["pdbs"]:
            raise RuntimeError("No models PDBs found for measurement")

        values = {}
        job_data["measurements"] = values

        if "evaluations" not in job_data or not job_data["evaluations"] \
            or not self.task_list:
            return job_data

        job_dir = Path(job_data["job_dir"])
        protein = mda.Universe(str(job_dir / job_data["models"]["pdbs"][0]))
        ligand = mda.Universe(str(job_dir / job_data["evaluations"][0]["pdb"]))

        for name, task in self.task_list.items():
            if task.set_system(protein, ligand):
                values[name] = []

            else:
                print("Could not apply measurement '{}' to variant '{}'".format(
                          name, job_data["variant"].name))

        models = job_data["models"]["pdbs"]
        evaluations = job_data["evaluations"]
        for model_pdb, evalued in zip(models, evaluations):
            protein.load_new(str(job_dir / model_pdb))
            ligand.load_new(str(job_dir / evalued["pdb"]))

            for ts in ligand.trajectory:
                for name in values:
                    self.task_list[name].compute()

            for name in values:
                task = self.task_list[name]
                values[name].append(task.data)
                task.clear()

        return job_data
