"""
"""


from .sequence import ProtSeq


class BaseBuilder:
    def __init__(self, parameters, database):
        self.parameters = parameters
        self.db = database
        self._initialized = False
        self.prot_id = None
        self.protein = None

    def initialize(self):
        self._initialized = True
        self.prot_id = self.db.register_protein(self.parameters["protein_name"])
        self.protein = ProtSeq(self.parameters["protein_name"], self.parameters["pdb_file"])
        self.special_initialize()

    def special_initialize(self):
        raise NotImplementedError("Child classes must implement this method")

    def next_job(self):
        if not self._initialized:
            self.initialize()  # Lazy initialization

        return self.fetch_next_job()

    def fetch_next_job(self):
        raise NotImplementedError("Child classes must implement this method")

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
                eval_data["pdb"]
            )
