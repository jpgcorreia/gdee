"""
"""


from .base_builder import BaseBuilder
from .sequence import BLOSUM62_AA, ResidueIndex
import itertools


class ExhaustiveBuilder(BaseBuilder):
    def __init__(self, parameters, database):
        super().__init__(parameters, database)
        self.combinations = None

    def special_initialize(self):
        self.variant = self.protein.copy()

        if not self.parameters["selection"]:
            raise RuntimeError("No residues were selected to be mutated")

        selection = ResidueIndex(self.protein, self.parameters["selection"])
        self.wildtype_sel = selection.apply(self.protein)
        self.variant_sel = selection.apply(self.variant)
        for wt, mut in zip(self.wildtype_sel, self.variant_sel):
            assert wt == mut # Order check

        self.mut_index = [res.index for res in self.variant_sel]
        self.mut_index.sort()

    def mutations(self):
        mutations = []
        for wt_res, mut_res in zip(self.wildtype_sel, self.variant_sel):
            if wt_res.code != mut_res.code:
                mutations.append("{}:{}{}{}".format(wt_res.chain, wt_res.code, wt_res.resid, mut_res.code))

        return "|".join(mutations)

    def fetch_next_job(self):
        if self.combinations is None:
            self.combinations = itertools.product(BLOSUM62_AA[:-4], repeat=len(self.variant_sel))

        mut_name = self.mutations()

        while self.db.variant_exists(self.prot_id, mut_name):
            try:
                mutations = next(self.combinations)
            except StopIteration:
                return None

            for mut_res, new_code in zip(self.variant_sel, mutations):
                mut_res.code = new_code

            mut_name = self.mutations()

        variant_dir = mut_name.replace("|", "_").replace(":", "")
        variant_id = self.db.register_variant(
            self.prot_id,
            mut_name,
            self.variant.to_modeller(),
            variant_dir,
            False if mut_name else True
        )
        variant = self.variant.copy()
        variant.name = mut_name

        job = {
            "variant_dir": variant_dir,
            "variant_id": variant_id,
            "wildtype": self.protein.copy(),
            "variant": variant,
            "mut_index": self.mut_index.copy()
        }
        return job
