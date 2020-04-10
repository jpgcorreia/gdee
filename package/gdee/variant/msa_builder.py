"""
"""


from .sequence import ProtSeq
from .base_builder import BaseBuilder
import re
import warnings

from Bio import SeqIO, Align, BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonWarning)
    from Bio.Align import substitution_matrices


def get_valid_filename(name):
    name = str(name).strip().replace(' ', '_')
    return re.sub(r"(?u)[^-\w.]", "", name)


class MSABuilder(BaseBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.wt_seq = ""
        self.msa = tuple()
        self._iter = iter(self.msa)
        self._is_wildtype = True
        self._aligner = Align.PairwiseAligner()
        self._aligner.open_gap_score = -10
        self._aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    def special_initialize(self):
        self.wt_seq = self.protein.to_modeller().replace("/", "")
        self.msa = tuple(SeqIO.parse(self.parameters["msa"], "fasta"))
        self._iter = iter(self.msa)

    def variant_from_alignment(self, name, other_seq):
        alignment = self._aligner.align(self.wt_seq, other_seq)
        if not alignment:
            raise RuntimeError("No suitable alignment found for sequence '{}'".format(other_seq))

        variant = self.protein.copy()
        variant.name = name
        variant_iter = iter(variant.flatten())
        mut_index = []

        for query, match, target in zip(*alignment[0].format().split()):
            if query == "-":
                continue

            seq_pos = next(variant_iter)
            code = seq_pos.code
            assert code == query

            if target != "-" and target != code:
                seq_pos.code = target
                mut_index.append(seq_pos.index)

        return variant, mut_index

    def fetch_next_job(self):
        is_wildtype = False
        if self._is_wildtype:
            self._is_wildtype = False
            is_wildtype = True
            variant = self.protein.copy()
            mut_index = []

        else:
            try:
                next_seq = next(self._iter)
                variant, mut_index = self.variant_from_alignment(next_seq.name, str(next_seq.seq))

            except StopIteration:
                return None

        variant_dir = get_valid_filename(variant.name)
        variant_id = self.db.register_variant(
            self.prot_id,
            variant.name,
            variant.to_modeller(),
            variant_dir,
            is_wildtype
        )

        job = {
            "variant_dir": variant_dir,
            "variant_id": variant_id,
            "wildtype": self.protein.copy(),
            "variant": variant,
            "mut_index": mut_index
        }
        return job
