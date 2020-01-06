"""
"""
from mpi4py import MPI
from enum import IntEnum, auto


class RequestType(IntEnum):
    get_task = 0
    save = auto()
    run_task = auto()
    terminate = auto()


class Message:
    def __init__(self, request, data=None):
        self.request = request
        self.data = data


class MPIPlatform:
    def __init__(self, pipeline):
        self.pipeline = pipeline
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.root = 0
        self.manager = None
        self.runner = None

    def run(self):
        if self.rank == self.root:
            self.manager = MPIManager(self.comm, self.pipeline)
            self.manager.run()

        else:
            self.runner = MPIRunner(self.comm, self.root, self.pipeline)
            self.runner.run()


class MPIManager:
    def __init__(self, mpi_comm, pipeline):
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.status = MPI.Status()
        self.pipeline = pipeline
        self.alive = mpi_comm.Get_size() - 1 # Number of runners
        self.not_saved = 0

    def run(self):
        while self.alive or self.not_saved:
            message = self.comm.recv(source=MPI.ANY_SOURCE, status=self.status)

            if message.request == RequestType.get_task:
                self.send_task()

            if message.request == RequestType.save:
                self.pipeline.save_results(message.data)
                self.not_saved -= 1

    def send_task(self):
        task = self.pipeline.next_job()

        if task is None:
            message = Message(RequestType.terminate)
            self.alive -= 1

        else:
            message = Message(RequestType.run_task, task)
            self.not_saved += 1

        runner = self.status.Get_source()
        self.comm.send(message, dest=runner)


class MPIRunner:
    def __init__(self, mpi_comm, root, pipeline):
        self.comm = mpi_comm
        self.rank = mpi_comm.Get_rank()
        self.root = root
        self.pipeline = pipeline

    def run(self):
        message = self.request_task()

        while message.request != RequestType.terminate:
            results = self.pipeline.run_pipeline(message.data)
            self.comm.send(Message(RequestType.save, results), dest=self.root)

            message = self.request_task()

    def request_task(self):
        return self.comm.sendrecv(
            Message(RequestType.get_task),
            dest=self.root,
            source=self.root
        )
