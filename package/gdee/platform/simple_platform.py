"""
"""


class SimplePlatform:
    def __init__(self, pipeline):
        self.pipeline = pipeline
        print("Simple platform created")

    def run(self):
        print("Running on Simple Platform")
        self.pipeline.run_initial()

        while self.pipeline.available():
            job_data = self.pipeline.next_job()
            results = self.pipeline.run_pipeline(job_data)
            self.pipeline.save_results(results)
