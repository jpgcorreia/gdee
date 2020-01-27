"""
"""


from path import Path
import modeller as mdl
import MDAnalysis as mda
from modeller import automodel
import shutil
from tempfile import TemporaryDirectory


class ModellerBuilder:
    def __init__(self, parameters):
        self.parameters = parameters
        mdl.log.none()
        self.env = mdl.environ()
        mdl.log.level(0, 0, 0, 0, 0)

    def run(self, job_data):
        temp_dir = TemporaryDirectory(prefix="gdee_modeller")
        temp_path = Path(temp_dir.name)
        Path(self.parameters["pdb_file"]).copy(temp_path / "template.pdb")

        with temp_path:
            self.write_alignment(job_data)

            model_data = self.build_models(job_data)
            if not model_data:
                job_data["fatal_error"] = True
                return job_data

            model_data.sort(key=lambda x: x["molpdf"])
            top_models = model_data[:self.parameters["num_models"]]

            pdb_list = [model["name"] for model in top_models]
            structure = mda.Universe(pdb_list[0], pdb_list)
            self.rename_models(structure, job_data["variant"])

            pdb_list = []
            for i, ts in enumerate(structure.trajectory):
                filename = "model_{:04d}.pdb".format(i)
                pdb_list.append(filename)
                structure.atoms.write(str(job_data["job_dir"] / filename))

        job_data["models"] = {
            "method": "modeller",
            "scores": [model["molpdf"] for model in top_models],
            "pdbs": pdb_list
        }

        return job_data

    def write_alignment(self, job_data):
        template = ">P1;template\nstructureX:template.pdb:.:.:.:.::::\n{}*\n\n>P1;model\nsequence:model.pdb:.:.:.:.::::\n{}*\n"
        with open("alignment.ali", "w") as fd:
            fd.write(template.format(
                job_data["wildtype"].to_modeller(),
                job_data["variant"].to_modeller()
            ))

    def build_models(self, job_data):
        model = MutationModel(
            self.env,
            alnfile="alignment.ali",
            knowns="template",
            sequence="model",
        )

        model.select_opt_residues(
            job_data["mut_index"],
            self.parameters["optimize_radius"]
        )
        model.starting_model = 1
        model.ending_model = 3 * self.parameters["num_models"]

        opt_level = self.parameters["optimize_level"]
        if opt_level == 0:
            model.library_schedule = automodel.autosched.very_fast
            model.md_level = automodel.refine.fast

        elif opt_level == 1:
            model.library_schedule = automodel.autosched.normal
            model.md_level = automodel.refine.slow

        else:
            model.library_schedule = automodel.autosched.slow
            model.md_level = automodel.refine.very_slow

        model.make()
        model_data = [data for data in model.outputs if "molpdf" in data]

        return model_data

    def rename_models(self, structure, prot_seq):
        res_idx = 0
        for chain in prot_seq:
            segment = structure.add_Segment(segid=chain.code)

            for seq_pos in chain:
                if seq_pos.is_gap:
                    continue

                residue = structure.residues[res_idx]
                residue.resid = seq_pos.resid
                residue.resname = seq_pos.resname
                residue.segment = segment
                res_idx += 1


class MutationModel(automodel.automodel):
    def select_opt_residues(self, residues, coff):
        self._opt_residues = tuple(map(int, residues))
        self._opt_residues_coff = coff

    def select_atoms(self):
        if not self._opt_residues or self._opt_residues_coff == 0:
            return mdl.selection(self.residues)

        res = [self.residues[pos] for pos in self._opt_residues]
        sel = mdl.selection(res)

        return sel.select_sphere(self._opt_residues_coff).by_residue()
