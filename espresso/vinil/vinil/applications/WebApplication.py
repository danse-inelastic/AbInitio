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

from luban.applications.WebApplication import WebApplication as base

class WebApplication(base):

    class Inventory(base.Inventory):
        import pyre.inventory

        blah    = pyre.inventory.str(name="blah", default="murmur")

    def _configure(self):
        super(WebApplication, self)._configure()

        self.blah = self.inventory.blah

    def _init(self):
        super(WebApplication, self)._init()


    def __init__(self, name):
        super(WebApplication, self).__init__(name=name)
        return

    def _getPrivateDepositoryLocations(self):
        return ['/tmp/luban-services', '../config', '../content']


__date__ = "$Nov 12, 2009 4:07:27 PM$"


