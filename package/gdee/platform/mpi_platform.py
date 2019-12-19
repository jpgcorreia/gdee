"""
"""


class RunnerState:
    def __init__(self):
        pass


class Message:
    def __init__(self):
        pass


class MPIPlatform:
    def __init__(self, pipeline):
        self.pipeline = pipeline
        print("MPI platform created")

    def run(self):
        print("Running on MPI platform")
        self.pipeline.run_initial()

        while self.pipeline.available():
            job_data = self.pipeline.next_job()
            results = self.pipeline.run_pipeline(job_data)
            self.pipeline.save_results(results)


class MPIManager:
    def __init__(self):
        pass


class MPIRunner:
    def __init__(self):
        pass

