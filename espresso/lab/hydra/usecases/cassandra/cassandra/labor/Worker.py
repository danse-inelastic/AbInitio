
from collections import deque
from twisted.spread.pb import PBClientFactory
from cassandra.labor.WorkerProxy import WorkerProxy

class Worker(PBClientFactory):

    def __init__(self):
        PBClientFactory.__init__(self)

        # Attribute defaults:
        self.uuid               = ""
	self.computeChunkSize	= 10
	self.computeInterval	= 1 # sec
	self.generator		= None  # functor
	self.outputBuffer	= []
	self.reportInterval	= 5 # sec
        self.taskQueue		= deque()
	self.scheduler		= None  # set in WorkerDaemon
	self.taskmaster         = None  # set in WorkerDaemon
        
        deferred    = self.getRootObject()
        deferred.addCallback(self.initProxy)


    def initProxy(self, proxyobj):
        "Init taskmaster"
        self.taskmaster     = proxyobj      # TaskMasterProxy
        self.taskmaster.callRemote("initProxy", WorkerProxy()) #
        # Init worker
        self.taskmaster.callRemote("workerInit")



