
from twisted.spread.pb import PBClientFactory
from cassandra.labor.WorkerProxy import WorkerProxy

class Worker(PBClientFactory):

    def __init__(self):
        PBClientFactory.__init__(self, WorkerProxy)
