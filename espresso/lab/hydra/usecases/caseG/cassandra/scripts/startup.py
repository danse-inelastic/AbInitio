
import os
import sys
import time
from twisted.internet import reactor

def start(config):
    os.chdir(config['basedir'])
    # No setting.tac file, no daemon start
    if not os.path.exists("setting.tac"):
        print "This doesn't look like a cassandra base directory:"
        print "No setting.tac file. Giving up!"
        sys.exit(1)

    
    # fork a child to launch the daemon, while the parent process tails the logfile
    if os.fork():
        # this is the parent
        print "Starting reactor ..."
        reactor.run()

    # this is the child: give the logfile-watching parent a chance to start
    # watching it before we start the daemon
    time.sleep(0.2)
    launch(config)


    
#    if (not os.path.exists("buildbot.tac") and
#        not os.path.exists("Makefile.buildbot")):
#        print "This doesn't look like a buildbot base directory:"
#        print "No buildbot.tac or Makefile.buildbot file."
#        print "Giving up!"
#        sys.exit(1)
#    if config['quiet']:
#        return launch(config)

#    # we probably can't do this os.fork under windows
#    from twisted.python.runtime import platformType
#    if platformType == "win32":
#        return launch(config)

#    # fork a child to launch the daemon, while the parent process tails the
#    # logfile
#    if os.fork():
#        # this is the parent
#        print "Starting reactor ..."
#        reactor.run()
#
#    # this is the child: give the logfile-watching parent a chance to start
#    # watching it before we start the daemon
#    time.sleep(0.2)
#    launch(config)


def launch(config):
    sys.path.insert(0, os.path.abspath(os.getcwd()))
    
    argv = ["twistd",
            "--no_save",
            "--logfile=twistd.log", # windows doesn't use the same default
            "--python=setting.tac"]
    sys.argv = argv

    from twisted.scripts import twistd
    twistd.run()

#    run = twistd.run
#    run()