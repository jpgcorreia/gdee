import os
from oddt.virtualscreening import virtualscreening as vs
import math
from ..misc import DataContainer
from ..database import Ranked_Database
import sys
import traceback
from pathlib import Path
import joblib
import numpy as np
from gdee.evaluator.vina import VinardoDocking

import warnings
warnings.filterwarnings("ignore", category=UserWarning)


class Rescore:
    DEFAULT_FUNCTIONS = {
        "rfscore_v1": "RFScore_v1_pdbbind2016.pickle",
        "rfscore_v2": "RFScore_v2_pdbbind2016.pickle",
        "rfscore_v3": "RFScore_v3_pdbbind2016.pickle",
        "plecnn": "PLECnn_p5_l1_pdbbind2016_s65536.pickle",
        "plecrf": "PLECrf_p5_l1_pdbbind2016_s65536.pickle",
        "pleclinear": "PLEClinear_p5_l1_pdbbind2016_s65536.pickle",
        "nnscore": "NNScore_pdbbind2016.pickle",
    }
    METAMODEL_NAME = "metamodel"
    METAMODEL_FUNCTIONS = (
        "rfscore_v1",
        "rfscore_v2",
        "pleclinear",
        "plecrf",
        "vinardo",
    )

    def __init__(self, output_db, table, use_metamodel=True, functions=None, data_dir=None, mgltools="mgltools", smina="smina"):
        self.table = table
        self.path = None
        self.output_db = Ranked_Database(output_db)
        self.output_db.table = self.table

        self.functions = self._normalize_functions(functions)
        self.use_metamodel = self._should_use_metamodel(self.functions, use_metamodel)
        self.scoring_method = self._get_scoring_method()
        self.data_dir = Path(data_dir) if data_dir else Path(__file__).resolve().parents[1] / "data"
        self.mgltools = mgltools
        self.smina = smina
        self.function_paths = self._resolve_function_paths(self.functions)
        self.meta_model = None

        if self.use_metamodel:
            self.meta_model = joblib.load(self.data_dir / "StackingEnsemble.joblib")
            self.function_paths = self._default_function_paths()

    def _normalize_functions(self, functions):
        if functions is None:
            return []

        if isinstance(functions, (str, Path)):
            return [functions]

        return list(functions)

    def _should_use_metamodel(self, functions, use_metamodel):
        if not functions:
            return use_metamodel

        return any(str(function).strip().lower() == self.METAMODEL_NAME for function in functions)

    def _get_scoring_method(self):
        """Return a stable label describing the scoring method for output rows."""
        if self.use_metamodel:
            return self.METAMODEL_NAME

        if not self.functions:
            return "rfscore_v1"

        function_name = str(self.functions[0]).strip()
        normalized_name = function_name.lower()
        if normalized_name.endswith("vinardo"):
            return "vinardo"

        if normalized_name in self.DEFAULT_FUNCTIONS:
            return normalized_name

        return "custom:{}".format(Path(function_name).name)

    def _default_function_paths(self):
        """Return default scoring function paths."""
        return [
            self.data_dir / self.DEFAULT_FUNCTIONS[name]
            for name in self.METAMODEL_FUNCTIONS
            if name != "vinardo"
        ] + ["vinardo"]

    def _resolve_function_paths(self, functions):
        if not functions:
            return self._default_function_paths() if not self.use_metamodel else []

        if len(functions) > 1:
            non_meta_functions = [function for function in functions if str(function).strip().lower() != self.METAMODEL_NAME]
            if len(non_meta_functions) > 1:
                raise ValueError("Select 'metamodel' or provide one scoring function")

        function = functions[0]
        if str(function).strip().lower() == self.METAMODEL_NAME:
            return self._default_function_paths()

        return [self._resolve_function_path(function)]

    def _resolve_function_path(self, function):
        function_name = str(function).strip()
        normalized_name = function_name.lower()

        if normalized_name in self.DEFAULT_FUNCTIONS:
            return self.data_dir / self.DEFAULT_FUNCTIONS[normalized_name]

        return Path(function_name)

    def get_variants(self, db):
        variant_list = list(db.fetch_ranked_variants(self.table))
        self.variant_it = iter(variant_list)

    def fetch_next_job(self):
        try:
            variant = next(self.variant_it)
        except StopIteration:
            return None

        name = variant[1]
        variant_dir = variant[2]
        wt = variant[3]
        model_id = variant[4]
        variant_id = variant[5]
        eval_id = variant[6]
        pose = variant[7]
        docking_file = variant[8]

        job = DataContainer()
        job.variant_dir = variant_dir
        job.name = name
        job.pose = pose
        job.docking_file = docking_file
        job.wt = wt
        job.model_id = model_id
        job.variant_id = variant_id
        job.eval_id = eval_id
        return job

    def run(self, job_data):
        variant_dir = self.path / job_data.variant_dir
        docking_file = variant_dir / job_data.docking_file
        model_file = variant_dir / "model_{}".format(job_data.docking_file[-8:])
        pose = job_data.pose
        name = job_data.name

        results = self.run_rescore(docking_file, model_file, pose, name)

        job_data.results = results
        return job_data

    def kd_to_energy(self, pkd):
        R = 8.314
        T = 298.15
        kd = pow(10, -(float(pkd)))
        j = R * T * (math.log(float(kd)))
        kj = j / 1000
        kcal = kj / 4.18
        return kcal

    def _score_vinardo(self, model_file, docking_file):
        """Return vinardo score for the given model and docking pose file."""
        # instantiate a VinardoDocking with minimal parameters.
        params = {
            "ligand": None,
            "box_center": [0, 0, 0],
            "box_size": [1, 1, 1],
            "exhaustiveness": 1,
            "vina": "vina",
            "vinardo": self.smina,
            "mgltools": self.mgltools,
            "atom_type": None,
        }
        vin = VinardoDocking(params)
        # docking_file may be PDB with ligand models; Rescore should pass a ligand PDBQT path here.
        # Assume docking_file is a ligand pdbqt; if not, user-supplied files_path must include original ligand.pdbqt.
        score = vin.score_pose(str(model_file), str(docking_file))
        return score

    def run_ml_functions(self, docking_file, model_file, pose, name):
        function_paths = self.function_paths or self._default_function_paths()
        ml_scores = []

        for function in function_paths:
            func_name = str(function).strip().lower()
            if func_name == "vinardo" or func_name.endswith("vinardo"):
                try:
                    val = self._score_vinardo(model_file, docking_file)
                    ml_scores.append(val)
                except Exception as e:
                    print("Vinardo scoring failed for {}: {}".format(name, e))
                    ml_scores.append(float("nan"))
                continue

            vs_rescore = vs(n_cpu=1)
            vs_rescore.load_ligands("pdb", docking_file)
            vs_rescore.score(function=str(function), protein=model_file)

            # Save result
            results = []
            for mol in vs_rescore.fetch():
                data = mol.data.to_dict()

                if len(data) > 0:
                    data['name'] = mol.title
                else:
                    print('There is no data', file=sys.stderr)
                    return False

                for key in data:
                    if key.startswith('rf'):
                        results.append(data[key])
                    elif key.startswith('PLEC'):
                        results.append(data[key])

            try:
                ml_scores.append(self.kd_to_energy(results[pose]))

            except Exception as error:
                print("Exception caught for variant {}:".format(name), traceback.print_exc(), error, sep="\n")
        return ml_scores

    def run_rescore(self, docking_file, model_file, pose, name):
        ml_scores = self.run_ml_functions(docking_file, model_file, pose, name)

        if self.use_metamodel:
            features = np.asarray([ml_scores], dtype=float)
            return float(self.meta_model.predict(features)[0])

        # fallback: single scoring function
        return ml_scores[0]

    def save_results(self, data):
        name = data.name
        variant_dir = os.path.basename(data.variant_dir)

        self.output_db.register_variant(
            data.results,
            name,
            variant_dir,
            data.wt,
            data.model_id,
            data.variant_id,
            data.eval_id,
            data.pose,
            data.docking_file,
            scoring_method=self.scoring_method,
        )

        print("Ended with variant '{}'".format(name))
        return True
