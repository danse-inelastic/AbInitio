
from twisted.spread.pb import Root  # Referenceable

class TaskMasterProxy(Root):

    def remote_initProxy(self, worker):
        #self.factory.
        pass

    def remote_workerInit(self):
        print "workerInit"

    def remote_workerStatus(self, status):
        print 'workerStatus'

    def remote_taskCompleted(self, task, output):
        print 'taskCompleted'