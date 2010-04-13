
import os
from twisted.python import log
from twisted.application.service import MultiService


class Daemon(MultiService):
    debug = 0

    def __init__(self, basedir, configFileName):
        MultiService.__init__(self)
        self._basedir           = basedir
        self._configFileName    = configFileName
        self._config            = None  # Configuration dictionary


    def startService(self):
        "Starts service"
        MultiService.startService(self)


    def stopService(self):
        "Stops service"
        MultiService.stopService(self)
        log.msg("stopping service")


    def _setConfig(self):
        configFile = os.path.join(self._basedir, self._configFileName)
        log.msg("loading configuration from %s" % configFile)

        try:
            f = open(configFile, "r")
        except IOError, e:
            log.msg("unable to open config file '%s'" % configFile)
            log.err(e)
            return

        try:
            self._config     = self._configDict(f)
        except:
            log.msg("error during parseConfig")
            log.err()
        f.close()


    def _configDict(self, file):
        "Parsing configuration file"
        localDict = {'basedir': os.path.expanduser(self._basedir)}
        try:
            exec file in localDict      # Ugly!
        except:
            log.msg("error while parsing config file")
            raise

        try:
            config = localDict['CassandraConfig']
            return config
        except KeyError:
            log.err("missing config dictionary")

        return None
