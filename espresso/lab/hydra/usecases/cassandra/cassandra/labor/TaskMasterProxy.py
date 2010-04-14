
from twisted.spread.pb import Root  # Referenceable

class TaskMasterProxy(Root):

    def remote_workerInit(self):
        print "workerInit"