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

from opal.content.Page import Page as Base

class Page(Base):

    def __init__(self, name, root, title=''):
        Base.__init__(self)

        head = self.head()
        head.base(url=root)
        head.title(title)

        head.stylesheet(rel="stylesheet", media="all", url="css/main.css")

        self.home = '%s/%s.html' % (root, name)

        return


# version
__id__ = "$Id$"

__date__ = "$Jul 20, 2009 7:48:10 AM$"

# End of file 