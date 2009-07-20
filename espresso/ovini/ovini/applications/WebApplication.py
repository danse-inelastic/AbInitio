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

    class Inventory(Base.Inventory):

        import opal.inventory
        import pyre.inventory

        # components
        actor = opal.inventory.actor(default='nyi')
        actor.meta['tip'] = "the component that defines the application behavior"


    def __init__(self, name):
        Base.__init__(self, name)
        


    def main(self, *args, **kwds):
        actor = self.actor
        if actor is None:
            inquiry = self.inventory._getTraitDescriptor('actor').inquiry
            actor = self.retrieveActor('nyi')
            actor.message = "Not implemented yet! actor=%s, routine=%s" % (
                inquiry, self.inventory.routine)
            self.actor = actor

        noErrors = True
        try:
            page = self.actor.perform(self, routine=self.inventory.routine, debug=self.debug)
            self.recordActivity()

            if isinstance(page, basestring):
                print page,
            else:
                self.render(page)
        except:
            noErrors = False
            try:
                self.fancyBugReport()
            except:
                # if we cannot generate a fancy report. we need a plain one
                self.plainBugReport()


    def retrievePage(self, name):
        page = super(WebApplication, self).retrievePage(name)
        if page:
            return page
        raise RuntimeError, "Unable to load page %s" % name

    def plainBugReport(self):
        print '<pre>'
        import traceback
        traceback.print_exc()
        print '</pre>'
        return

    def _configure(self):
        super(WebApplication, self)._configure()

        # custom weaver
        import os
        configurations = {
            'home': self.home,
            'cgihome':self.cgihome,
            }
        import vnf.weaver
        self.pageMill = vnf.weaver.pageMill(configurations)


        self.idd = self.inventory.idd
        self.clerk = self.inventory.clerk
        # this next line is a problem.  Technically, many of the components can be None at
        # this point....so trying to set an attribute of a None-type component throws an
        # exception....root of the problem may be in initializeConfiguration() in Application.py
        self.clerk.director = self
        self.dds = self.inventory.dds
        # same for this line
        self.dds.director = self
        self.scribe = self.inventory.scribe
        self.debug = self.inventory.debug
        self.csaccessor = self.inventory.csaccessor
        self.itaskmanager = self.inventory.itaskmanager

        from vnf.components import accesscontrol
        self.accesscontrol = accesscontrol()

        return



__date__ = "$Jul 19, 2009 11:07:10 PM$"


