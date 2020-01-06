"""
"""


class Database:
    def __init__(self, filename):
        self.filename = filename
        self._conn = None

    def connect(self):
        print("Database initialized on file", self.filename) # TODO

    @property
    def conn(self):
        # Database is connected only when needed to allow
        # for multiple instantiations on MPI platform
        if self._conn is None:
            self.connect()

        return self._conn

