
from twisted.spread.pb import Root, Referenceable

class WorkerProxy(Referenceable):

    def remote_perform(self, workload):
        print 'perform'


    def remote_createComponent(self, name, factory, registry):
        print 'createComponent'


    def remote_finalizeComponent(self):
        print 'finalizeComponent'


    def remote_reset(self):
        print 'reset'


#    def remote_workerInit(self):
#        print "workerInit"