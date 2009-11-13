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

from luban.applications.SuperApp import SuperApp as base

class SuperApp(base):

    def runApp(self, config=None):

        if not config:
            config = ['/tmp/luban-services', '../config', '../content/components']

        from luban.applications.WebApplication import WebApplication

        class MainApp(WebApplication):

            # Adding test doesn't work
            class Inventory(WebApplication.Inventory):
                import pyre.inventory

                test    = pyre.inventory.str(name="test", default="murmur")

            def _configure(self):
                super(MainApp, self)._configure()

                self.test = self.inventory.test

            def _init(self):
                super(MainApp, self)._init()


            def __init__(self):
                WebApplication.__init__(self, name='main')
                return

            def _getPrivateDepositoryLocations(self):
                return config

        app = MainApp()
        return app.run()


__date__ = "$Nov 12, 2009 4:07:27 PM$"


