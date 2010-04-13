
from twisted.spread.pb import PBServerFactory
from cassandra.labor.TaskMasterProxy import TaskMasterProxy

class TaskMaster(PBServerFactory):

    def __init__(self):
        PBServerFactory.__init__(self, TaskMasterProxy)

    