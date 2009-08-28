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

# Not used

from opal.content.Body import Body as Base

class Body(Base):
    def __init__(self, **kwds):
        Base.__init__(self, **kwds)
        #self._header = None

    def pageHeader(self, **kwds):
        from PageHeader import PageHeader
        self._header = PageHeader(**kwds)
        self.contents.append(self._header)
        return self._header

__date__ = "$Jul 26, 2009 4:54:44 PM$"


