
from twisted.application.service import MultiService

class TaskMaster(MultiService):
    debug = 0

    def __init__(self, basedir, configFile = "taskmaster.cfg"):
        MultiService.__init__(self)
        self._basedir       = basedir
        self._configFile    = configFile

    def startService(self):
        MultiService.startService(self)
