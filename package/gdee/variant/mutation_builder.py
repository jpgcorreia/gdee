"""
"""


from .sequence import ProtSeq, MatrixMutation, Blosum62Mutation
from collections import defaultdict


class MutationBuilder:
    def __init__(self, parameters, database):
        self.parameters = parameters
        self.db = database
        self._initialized = False
        self.prot_id = None
        self.protein = None
        self.variant = None
        self.mut_sel = []
        self.invert_weights = self.parameters["conservative"]
        self.max_iter = self.parameters["max_iterations"]
        self.iterations = 0

        if self.parameters["matrix"] == "blosum62":
            self.matrix = Blosum62Mutation()

        else:
            self.matrix = MatrixMutation(
                self.parameters["matrix_aa"],
                self.parameters["matrix_weights"]
            )

    def initialize(self):
        self._initialized = True
        self.prot_id = self.db.register_protein(self.parameters["protein_name"])
        self.protein = ProtSeq(self.parameters["protein_name"], self.parameters["pdb_file"])
        self.variant = self.protein.copy()

        if not self.parameters["selection"]:
            raise RuntimeError("No residues were selected to be mutated")

        self.mut_sel = self.select(self.parameters["selection"])
        self.mut_index = [res.index for res in self.mut_sel]
        self.mut_index.sort()

    def select(self, selection):
        selected = []
        chain_sel = defaultdict(set)

        for residue in selection.split(" "):
            chain, resid = residue.split(":")
            chain_sel[chain].add(int(resid))

        for chain, sel_list in chain_sel.items():
            for residue in self.variant[chain]:
                if residue.resid in sel_list:
                    selected.append(residue)
                    sel_list.remove(residue.resid)

            if sel_list:
                not_found = ", ".join(map(str, sorted(sel_list)))
                raise RuntimeError("Residues not found in chain '{}': {}".format(chain, not_found))

        selected.sort(key=lambda x: x.index)

        return selected

    def mutations(self):
        return "|".join("{}:{}:{}".format(res.chain, res.resid, res.code) for res in self.mut_sel)

    def next_job(self):
        if self.iterations >= self.max_iter:
            return None

        if not self._initialized:
            self.initialize()  # Lazy initialization

        wildtype = True # First one is always wildtype
        mut_name = self.mutations()

        while self.db.variant_exists(self.prot_id, mut_name):
            for residue in self.mut_sel:
                residue.code = self.matrix.mutate(residue.code, self.invert_weights)

            mut_name = self.mutations()
            wildtype = False

        variant_dir = mut_name.replace("|", "_").replace(":", "")
        variant_id = self.db.register_variant(self.prot_id, mut_name, variant_dir, wildtype)
        variant = self.variant.copy()
        variant.name = mut_name
        self.iterations += 1

        job = {
            "variant_dir": variant_dir,
            "variant_id": variant_id,
            "wildtype": self.protein.copy(),
            "variant": variant,
            "mut_index": self.mut_index.copy()
        }
        return job

    def remove_variant(self, data):
        self.db.remove_variant(data["variant_id"])

    def save_results(self, data):
        variant_id = data["variant_id"]
        models = data["models"]

        for model_idx in range(len(models["pdbs"])):
            model_id = self.db.register_model(
                variant_id,
                models["method"],
                [models["scores"][model_idx]],
                models["pdbs"][model_idx]
            )

            eval_data = data["evaluations"][model_idx]

            eval_id = self.db.register_evaluation(
                variant_id,
                model_id,
                eval_data["ligand_file"],
                eval_data["method"],
                eval_data["energies"],
                [eval_data["pdb"]]
            )
