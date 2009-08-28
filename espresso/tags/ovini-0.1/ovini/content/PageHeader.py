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

from opal.content.PageHeader import PageHeader as Base
#from opal.content.PageSection import PageSection as Base

class PageHeader(Base):
    def __init__(self): # No attributes can be passed!
        Base.__init__(self) #, id="head-wrapper"

    def identify(self, inspector):
        return inspector.onPageHeader(self)

        

__date__ = "$Jul 26, 2009 4:46:23 PM$"


