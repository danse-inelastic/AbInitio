#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.applications.Script import Script as base

class HelloApp(base):

    def main(self):
        print 'Content-type: text/plain'
        print
        print "Hello World!"
        return

    def __init__(self, name='hello1'):
        super(HelloApp, self).__init__(name=name)
        return

# main
if __name__ == "__main__":
    app = HelloApp()
    app.run()

"""
print 'Content-type: text/plain'
print
print 'Hello, Alex!'
"""

__date__ = "$Jul 19, 2009 9:05:08 AM$"

# End of file
