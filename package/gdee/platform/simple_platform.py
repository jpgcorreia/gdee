"""
"""
import socket
import glob
import os

def _cleanup_gpu_lock_files(lock_dir="/tmp", lock_prefix="gdee_vina_gpu"):
    host = socket.gethostname()
    pattern = os.path.join(lock_dir, f"{lock_prefix}_{host}_*.lock")
    for p in glob.glob(pattern):
        try:
            os.remove(p)
        except OSError:
            pass

class SimplePlatform:
    def __init__(self, parameters, pipeline):
        self.pipeline = pipeline

    def run(self):
        try:
            prepare = getattr(self.pipeline, "prepare", None)
            if callable(prepare):
                prepare()

            job_data = self.pipeline.next_job(1)

            while job_data:
                results = [self.pipeline.run_pipeline(data) for data in job_data]
                self.pipeline.save_results(results)
                job_data = self.pipeline.next_job(1)
        finally:
            try:
                self.pipeline.finalize()
            finally:
                _cleanup_gpu_lock_files(lock_dir="/tmp", lock_prefix="gdee_vina_gpu")
