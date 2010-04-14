
from twisted.spread.pb import PBClientFactory
from cassandra.labor.WorkerProxy import WorkerProxy

class Worker(PBClientFactory):

    def __init__(self):
        PBClientFactory.__init__(self, WorkerProxy())
        
        deferred    = self.getRootObject()
        deferred.addCallback(self.init)


    def init(self, proxyobj):
        "Init taskmaster"
        self.taskmaster     = proxyobj	# TaskMasterProxy
        self.taskmaster.callRemote("workerInit")



