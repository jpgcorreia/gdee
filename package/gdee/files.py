"""
"""
import tarfile
from path import Path
import multiprocessing as mp


class Archiver:
    def __init__(self, archive, buffer_size):
        self.archive = archive
        self.buffer_size = buffer_size
        self._buffer = []
        self._job = None

    def add(self, path, new_name):
        if len(self._buffer) > self.buffer_size:
            self.flush()

        self._buffer.append((path, new_name))

    def finalize(self):
        self.flush()
        self.join()

    def join(self):
        if self._job is not None:
            self._job.join()

    def flush(self):
        self.join()

        if not self._buffer:
            return

        self._job = mp.Process(target=self._save_results, args=(self._buffer,))
        self._job.start()
        self._buffer = []

    def _save_results(self, data):
        with tarfile.open(self.archive, "a") as tar:
            while data:
                path, new_name = data.pop()
                tar.add(path, new_name)
                Path(path).rmtree_p()
