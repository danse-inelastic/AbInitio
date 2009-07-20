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

from opal.applications.WebApplication import WebApplication as Base

class WebApplication(Base):
    def __init__(self, name):
        #Base.__init__(self, name)
        
        print 'Content-type: text/plain'
        print
        print 'Hello, WebApplication!'


    def run(self):
        return True


__date__ = "$Jul 19, 2009 11:07:10 PM$"


