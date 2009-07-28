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

from opal.weaver.BodyMill import BodyMill as base

class BodyMill(base):
    def __init__(self, tagger):
        base.__init__(self, tagger)

    def onPageHeader(self, pageHeader):
        text = ['<h2><a href="/">Ovini</a> (Opal VNF Mini)</h2>']
        return text



__date__ = "$Jul 26, 2009 6:08:41 PM$"


