
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

    # One can fork process to trace logs. 
#    # Fork a child to launch the daemon, while the parent process tails the logfile
#    if os.fork():
#        # this is the parent
#        print "Starting reactor ..."
#        reactor.run()
#        reactor.stop()

    # This is the child: give the logfile-watching parent a chance to start
    # watching it before we start the daemon
    time.sleep(0.2)
    launch(config)


def launch(config):
    sys.path.insert(0, os.path.abspath(os.getcwd()))
    
    argv = ["twistd",
            "--no_save",
            "--logfile=twistd.log", # windows doesn't use the same default
            "--python=setting.tac"]
    sys.argv = argv

    from twisted.scripts import twistd
    twistd.run()
