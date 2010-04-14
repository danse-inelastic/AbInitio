
from collections import deque
from twisted.spread.pb import PBServerFactory

from cassandra.labor.TaskMasterProxy import TaskMasterProxy
from cassandra.labor.WorkerStatusRing import WorkerStatusRing

class TaskMaster(PBServerFactory):

    def __init__(self):
        PBServerFactory.__init__(self, TaskMasterProxy())

        # XXX: Some parameters move to config file
	self.bordomCounter      = 0
	self.bordomDelay 	= 20    # sec
	self.killDelay		= 10    # sec
	self.refuelingThreshold = 3
	self.scheduler          = None
	self.taskId             = 0     # unique ID for each task
	self.taskQueue          = deque()   # the global pile of work to do
	self.tickInterval 	= 1     # sec
	self.workerDB           = {}    # worker database is simply a dictionary of Employment records
	self.workerStatusRing   = WorkerStatusRing(self.killDelay)



    