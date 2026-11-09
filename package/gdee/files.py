"""
"""
import tarfile
import shutil
import tempfile
from glob import glob
from pathlib import Path as PathlibPath
from path import Path
import multiprocessing as mp


class InputFiles:
    """Resolve a rescoring directory or tar archive for one worker."""

    def __init__(self, source):
        self.source = str(source) if source is not None else None
        self.sources = None
        self.path = None
        self._temporary_path = None

    def prepare(self):
        if self.path is not None:
            return self.path
        if self.source is None:
            raise ValueError("An input directory or tar archive is required")
        source_paths = [PathlibPath(path) for path in sorted(glob(self.source))]
        if not source_paths:
            source_path = PathlibPath(self.source)
            if any(character in self.source for character in "*?["):
                raise FileNotFoundError(self.source)
            source_paths = [source_path]
        if len(source_paths) == 1 and source_paths[0].is_dir():
            self.sources = source_paths
            self.path = source_paths[0]
            return self.path
        if any(not path.exists() for path in source_paths):
            raise FileNotFoundError(self.source)
        if any(not tarfile.is_tarfile(path) for path in source_paths):
            raise ValueError("Input files must be directories or tar archives")

        extraction_path = PathlibPath(tempfile.mkdtemp(prefix="gdee_rescore_"))
        try:
            members = []
            member_names = set()
            member_files = set()
            for source_path in source_paths:
                with tarfile.open(source_path, "r:*") as archive:
                    archive_members = archive.getmembers()
                    self._validate_members(archive_members)
                    for member in archive_members:
                        if member.name in member_names and (member.name in member_files or not member.isdir()):
                            raise ValueError("Duplicate archive member: {}".format(member.name))
                        member_names.add(member.name)
                        if member.isfile():
                            member_files.add(member.name)
                    members.extend(archive_members)
                    archive.extractall(extraction_path)
            self.sources = source_paths
            self._temporary_path = extraction_path
            self.path = self._find_root(extraction_path, members)
            return self.path
        except Exception:
            shutil.rmtree(extraction_path, ignore_errors=True)
            raise

    @staticmethod
    def _validate_members(members):
        for member in members:
            member_path = PathlibPath(member.name)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise ValueError("Unsafe archive member: {}".format(member.name))
            if member.issym() or member.islnk() or not (member.isdir() or member.isfile()):
                raise ValueError("Unsupported archive member: {}".format(member.name))

    @staticmethod
    def _find_root(extraction_path, members):
        top_level = {
            PathlibPath(member.name).parts[0]
            for member in members
            if PathlibPath(member.name).parts
        }
        if len(top_level) == 1:
            candidate = extraction_path / next(iter(top_level))
            if candidate.is_dir():
                return candidate
        return extraction_path

    def cleanup(self):
        if self._temporary_path is not None:
            shutil.rmtree(self._temporary_path, ignore_errors=True)
            self._temporary_path = None
        self.path = None


class Archiver:
    def __init__(self, name, format, buffer_size):
        self.archive = Path(name)
        self.format = format + ".tar"
        self.buffer_size = buffer_size
        self._buffer = []
        self._job = None

        self._index = 0
        while self._next_filename().exists():
            pass
        self._index -= 1

    def _next_filename(self):
        new_file = self.archive + self.format.format(self._index)
        self._index += 1
        return new_file

    def add(self, path, new_name):
        if len(self._buffer) >= self.buffer_size:
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

        self._job = mp.Process(target=self._save_results, args=(self._next_filename(), self._buffer))
        self._job.start()
        self._buffer = []

    def _save_results(self, file_name, data):
        with tarfile.open(file_name, "a") as tar:
            while data:
                path, new_name = data.pop()
                tar.add(path, new_name)
                Path(path).rmtree_p()
