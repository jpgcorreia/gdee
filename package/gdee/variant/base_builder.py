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
        self._excluded = parameters["excluded"]
        self._variants = set()

    def is_excluded(self, residue, code):
        key = "{}:{}".format(residue.chain, residue.resid)
        return code in self._excluded.get(key, "")

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

    def variant_exists(self, name):
        if name in self._variants:
            return True

        return self.db.variant_exists(self.prot_id, name)

    def new_variant(self, name):
        if name in self._variants:
            raise RuntimeError("Variant {} already exists".format(name))

        self._variants.add(name)

    def save_results(self, data):
        name = data.variant.name
        if data.fatal_error:
            self._variants.remove(name)
            print("Error while processing variant: {}".format(name))
            return

        variant_id = self.db.register_variant(self.prot_id, name,
                                              data.variant.to_modeller(),
                                              data.variant_dir,
                                              data.is_wildtype)

        modeling = data.modeling
        for model in modeling.models:
            model_id = self.db.register_model(variant_id, modeling.method,
                                              model.score, model.pdb)

            for ligand_name, evaluation in model.evals.items():
                eval_id = self.db.register_evaluation(variant_id, model_id,
                                                      evaluation)
                pose_id_list = self.db.register_poses(eval_id, evaluation.energies)

                self.db.register_measurements(eval_id, pose_id_list,
                                              evaluation.measurements)

        print("Ended with variant '{}'".format(data.variant.name))

        return True
