"""
"""


import sqlite3 as sql
import os


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
                "        mutations TEXT UNIQUE NOT NULL,"
                "        prot_id"
                "            REFERENCES Proteins(prot_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        is_wildtype INT NOT NULL,"
                "        pdb_file TEXT,"
                "        pdb_code TEXT"
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
                "        pdb_files TEXT NOT NULL"
                "    );"
                ""
                "CREATE TABLE IF NOT EXISTS"
                "    Ligands ("
                "        ligand_id INTEGER PRIMARY KEY,"
                "        smiles TEXT UNIQUE NOT NULL,"
                "        pdb_file TEXT NOT NULL,"
                "        pdbqt_file TEXT NOT NULL,"
                "        atom_names TEXT NOT NULL,"
                "        atom_coords TEXT NOT NULL"
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
                "        ligand_id"
                "            REFERENCES Ligands(ligand_id)"
                "                ON DELETE CASCADE"
                "                ON UPDATE CASCADE,"
                "        method TEXT NOT NULL,"
                "        energy REAL NOT NULL,"
                "        pdb_file TEXT NOT NULL"
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

    def variant_exists(self, mutations):
        cursor = self.conn.execute(
            "SELECT EXISTS ("
            "    SELECT"
            "        1"
            "    FROM"
            "        Variants"
            "    WHERE"
            "        mutations = ?"
            ");",
            (mutations,)
        )
        return bool(cursor.fetchone()[0])

    def register_variant(self, prot_id, mutations, wildtype, pdb_file=None, pdb_code=None):
        conn = self.conn
        cursor = conn.execute(
            "INSERT INTO"
            "    Variants ("
            "        mutations,"
            "        prot_id,"
            "        is_wildtype,"
            "        pdb_file,"
            "        pdb_code"
            "    ) "
            "VALUES (?, ?, ?, ?, ?);",
            (mutations, prot_id, bool(wildtype), pdb_file, pdb_code)
        )
        conn.commit()

        return cursor.lastrowid
