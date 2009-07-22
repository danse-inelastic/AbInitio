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


        debug = pyre.inventory.bool(name="debug", default=True)
        debug.meta['tip'] = "suppress some html output for debugging purposes"

    def __init__(self, name):
        Base.__init__(self, name)

        # debugging mode
        self.debug = False


    def main(self, *args, **kwds):
        actor = self.actor
        if actor is None:
            inquiry = self.inventory._getTraitDescriptor('actor').inquiry
            self.actor = self.retrieveActor('nyi')
            print self.actor
            
            #actor.message = "Not implemented yet! actor=%s, routine=%s" % (
            #    inquiry, self.inventory.routine)

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

    def recordActivity(self):
        pass
        """
        from vnf.dom.Activity import Activity
        activity = Activity()

        from vnf.components.misc import new_id
        activity.id =  new_id(self)

        activity.actor = self.actor.name

        activity.username = self.sentry.username

        activity.routine = self.inventory.routine

        activity.remote_address = self._cgi_inputs.get('REMOTE_ADDR') or 'local'

        self.clerk.newRecord(activity)

        return
        """

    def _configure(self):
        super(WebApplication, self)._configure()

        # custom weaver
        import os
        configurations = {
            'home': self.home,
            'cgihome':self.cgihome,
            }
        import ovini.weaver
        self.pageMill = ovini.weaver.pageMill(configurations)

        self.debug = self.inventory.debug
        
        return

    def _init(self):
        super(WebApplication, self)._init()

    def _getPrivateDepositoryLocations(self):
        from os.path import join
        root = '..'
        content = join(root, 'content')
        config = join(root, 'config')

        #from ovini.depositories import depositories

        return [config, content] #depositories(content)+[config]

if __name__=='__main__':
    w = WebApplication(name = 'test')
    print w


__date__ = "$Jul 19, 2009 11:07:10 PM$"


