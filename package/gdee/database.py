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
                "        energies TEXT NOT NULL,"
                "        pdb_files TEXT NOT NULL,"
                "        UNIQUE("
                "            model_id,"
                "            ligand_file,"
                "            method"
                "        )"
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

    def register_evaluation(self, variant_id, model_id, ligand_file, method, energies, pdb_list):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Evaluations ("
            "        variant_id,"
            "        model_id,"
            "        ligand_file,"
            "        method,"
            "        energies,"
            "        pdb_files"
            "    ) "
            "VALUES (?, ?, ?, ?, ?, ?);",
            (variant_id, model_id, ligand_file, method, list_serialize(energies), list_serialize(pdb_list))
        )

        conn.commit()
        return cursor.lastrowid
