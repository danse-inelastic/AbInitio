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

def main():


    from opal.applications.WebApplication import WebApplication


    class MainApp(WebApplication):


        def __init__(self):
            WebApplication.__init__(self, name='main')
            return


    app = MainApp()
    return app.run()

# main
if __name__ == '__main__':
    # invoke the application shell
    try:
        main()
    except:
        import traceback
        import time
        t = time.ctime()
        messages = traceback.format_exc().split('\n')

"""
# Simple web application that prints "Hello World!" using pyre

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

"""
# A very basic web application

print 'Content-type: text/plain'
print
print 'Hello, Alex!'
"""

__date__ = "$Jul 19, 2009 9:05:08 AM$"

# End of file
