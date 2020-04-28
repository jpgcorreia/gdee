"""
"""


import sqlite3 as sql
import os


def list_serialize(values):
    try:
        return "|".join(map(lambda x: "{:.4f}".format(x), values))
    except ValueError:
        pass

    return "|".join(map(lambda x: "{}".format(x), values))


class Database:
    def __init__(self, filename):
        if os.path.splitext(filename)[1] != ".sqlite3":
            filename += ".sqlite3"

        self.filename = filename
        self._conn = None
        self._metric_ids = {}

    def connect(self):
        # Connect to database
        self._conn = sql.connect(self.filename)
        self._conn.row_factory = sql.Row

    def create_tables(self):
        # Create tables
        with self.conn:
            self.conn.executescript(
                "CREATE TABLE IF NOT EXISTS"
                "    Proteins ("
                "        prot_id INTEGER PRIMARY KEY,"
                "        uniprot TEXT,"
                "        name TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    ProteinMetadata ("
                "        protmeta_id INTEGER PRIMARY KEY,"
                "        prot_id"
                "            REFERENCES Proteins(prot_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        data TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Variants ("
                "        variant_id INTEGER PRIMARY KEY,"
                "        name TEXT NOT NULL,"
                "        sequence TEXT NOT NULL,"
                "        prot_id"
                "            REFERENCES Proteins(prot_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        directory TEXT,"
                "        is_wildtype INT NOT NULL,"
                "        pdb_file TEXT,"
                "        pdb_code TEXT,"
                "        UNIQUE("
                "            name,"
                "            sequence,"
                "            prot_id"
                "        )"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Models ("
                "        model_id INTEGER PRIMARY KEY,"
                "        variant_id"
                "            REFERENCES Variants(variant_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        method TEXT NOT NULL,"
                "        scores TEXT NOT NULL,"
                "        pdb_file TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Evaluations ("
                "        eval_id INTEGER PRIMARY KEY,"
                "        variant_id"
                "            REFERENCES Variants(variant_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        model_id"
                "            REFERENCES Models(model_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        ligand_file TEXT NOT NULL,"
                "        method TEXT NOT NULL,"
                "        pdb_file TEXT NOT NULL,"
                "        UNIQUE("
                "            model_id,"
                "            ligand_file,"
                "            method"
                "        )"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Poses ("
                "        pose_id INTEGER PRIMARY KEY,"
                "        eval_id"
                "            REFERENCES Evaluations(eval_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        pdb_index INTEGER NOT NULL,"
                "        energy FLOAT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Metrics ("
                "        metric_id INTEGER PRIMARY KEY,"
                "        name TEXT UNIQUE NOT NULL,"
                "        identifier TEXT UNIQUE NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Measurements ("
                "        measurement_id INTEGER PRIMARY KEY,"
                "        metric_id"
                "            REFERENCES Metrics(metric_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        pose_id"
                "            REFERENCES Poses(pose_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        value REAL NOT NULL"
                "    );"
        )

    @property
    def conn(self):
        # Database is connected only when needed to allow
        # for multiple instantiations on MPI platform
        if self._conn is None:
            self.connect()
            self.create_tables()

        return self._conn

    def register_protein(self, name, uniprot=None):
        conn = self.conn
        cursor = conn.execute(
            "SELECT"
            "    prot_id "
            "FROM"
            "    Proteins "
            "WHERE"
            "    uniprot IS ?"
            "    AND"
            "    name IS ?;",
            (uniprot, name)
        )
        data = cursor.fetchone()

        if data:
            return data[0]

        cursor = conn.execute(
            "INSERT INTO"
            "    Proteins ("
            "        uniprot,"
            "        name"
            "    ) "
            "VALUES (?, ?);",
            (uniprot, name)
        )
        conn.commit()

        return cursor.lastrowid

    def variant_exists(self, prot_id, mutations):
        cursor = self.conn.execute(
            "SELECT EXISTS ("
            "    SELECT"
            "        1"
            "    FROM"
            "        Variants"
            "    WHERE"
            "        prot_id = ?"
            "        AND"
            "        name = ?"
            ");",
            (prot_id, mutations,)
        )
        return bool(cursor.fetchone()[0])

    def register_variant(self, prot_id, name, sequence, directory, wildtype, pdb_file=None, pdb_code=None):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Variants ("
            "        name,"
            "        sequence,"
            "        prot_id,"
            "        directory,"
            "        is_wildtype,"
            "        pdb_file,"
            "        pdb_code"
            "    ) "
            "VALUES (?, ?, ?, ?, ?, ?, ?);",
            (name, sequence, prot_id, directory, bool(wildtype), pdb_file, pdb_code)
        )
        conn.commit()

        return cursor.lastrowid

    def remove_variant(self, variant_id):
        conn = self.conn
        cursor = conn.execute(
            "DELETE FROM"
            "    Variants "
            "WHERE"
            "    variant_id = ?;",
            (variant_id,)
        )
        conn.commit()

        return cursor.lastrowid

    def register_model(self, variant_id, method, scores, pdb_file):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Models ("
            "        variant_id,"
            "        method,"
            "        scores,"
            "        pdb_file"
            "    ) "
            "VALUES (?, ?, ?, ?);",
            (variant_id, method, list_serialize(scores), pdb_file)
        )

        conn.commit()
        return cursor.lastrowid

    def register_evaluation(self, variant_id, model_id, ligand_file, method, pdb_file):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Evaluations ("
            "        variant_id,"
            "        model_id,"
            "        ligand_file,"
            "        method,"
            "        pdb_file"
            "    ) "
            "VALUES (?, ?, ?, ?, ?);",
            (variant_id, model_id, ligand_file, method, pdb_file)
        )
        return cursor.lastrowid

    def register_poses(self, eval_id, energy):
        pose_id = []
        conn = self.conn
        cursor = conn.cursor()
        for index, value in enumerate(energy):
            cursor.execute(
                "INSERT INTO"
                "    Poses ("
                "        eval_id,"
                "        pdb_index,"
                "        energy"
                "    ) "
                "VALUES (?, ?, ?);",
                (eval_id, index, value)
            )
            pose_id.append(cursor.lastrowid)

        conn.commit()
        return pose_id

    def register_metric(self, name, identifier):
        if name not in self._metric_ids:
            conn = self.conn
            cursor = conn.execute(
                "SELECT"
                "    metric_id "
                "FROM"
                "    Metrics "
                "WHERE"
                "    name = ?;",
                (name,)
            )
            data = cursor.fetchall()
            if data:
                self._metric_ids[name] = data[0][0]

            else:
                cursor = conn.execute(
                    "INSERT INTO"
                    "    Metrics ("
                    "        name,"
                    "        identifier"
                    "    ) "
                    "VALUES (?, ?);",
                    (name, identifier)
                )
                self._metric_ids[name] = cursor.lastrowid
                conn.commit()

        return self._metric_ids[name]

    def register_measurements(self, eval_id, metric_list, pose_id_list, measurements):
        conn = self.conn
        cursor = conn.cursor()
        for name, values_list in zip(metric_list, measurements):
            metric_id = self._metric_ids[name]

            for pose_id, value in zip(pose_id_list, values_list):
                cursor.execute(
                    "INSERT INTO"
                    "    Measurements ("
                    "        metric_id,"
                    "        pose_id,"
                    "        value"
                    "    ) "
                    "VALUES (?, ?, ?);",
                    (metric_id, pose_id, float(value))
                )

        conn.commit()
