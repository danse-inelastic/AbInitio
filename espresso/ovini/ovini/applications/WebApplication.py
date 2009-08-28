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
        actor = opal.inventory.actor(default='greet')
        actor.meta['tip'] = "the component that defines the application behavior"

        import pyre.idd
        idd = pyre.inventory.facility('idd-session', factory=pyre.idd.session, args=['idd-session'])
        idd.meta['tip'] = "access to the token server"

        from ovini.components import clerk
        clerk = pyre.inventory.facility(name="clerk", factory=clerk)
        clerk.meta['tip'] = "the component that retrieves data from the various database tables"

        debug = pyre.inventory.bool(name="debug", default=True)
        debug.meta['tip'] = "suppress some html output for debugging purposes"

        imagepath   = pyre.inventory.str(name='imagepath', default = 'images' )

    def __init__(self, name):
        Base.__init__(self, name)

        # access to the token server
        self.idd = None

        # access to the data retriever
        self.clerk = None

        # debugging mode
        self.debug = False


    def main(self, *args, **kwds):
        actor = self.actor
        
        # If no actor is provided e.g.: example.com/main.cgi
        if actor is None:
            self.actor = self.retrieveActor('greet')

        try:
            page = self.actor.perform(self, routine=self.inventory.routine, debug=self.debug)

            if isinstance(page, basestring):
                print page
            else:
                self.render(page)
        except:
            self.plainBugReport()

    def retrievePage(self, name):
        page = super(WebApplication, self).retrievePage(name)
        if page:
            return page
        raise RuntimeError, "Unable to load page %s" % name

    def redirect(self, actor, routine, **kwds):
        self.inventory.routine = routine
        self.actor = self.retrieveActor(actor)

        if self.actor is not None:
            self.configureComponent(self.actor)
            for k,v in kwds.iteritems():
                setattr(self.actor.inventory, k, v)

        try:
            self.main()
        except:
            raise RuntimeError, "redirect to actor %r, routine %r, with kwds %r failed" % (
                actor, routine, kwds)
        return

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
            'imagepath':self.inventory.imagepath
            }
        import ovini.weaver
        self.pageMill = ovini.weaver.pageMill(configurations)

        self.idd = self.inventory.idd
        self.clerk = self.inventory.clerk
        # this next line is a problem.  Technically, many of the components can be None at
        # this point....so trying to set an attribute of a None-type component throws an
        # exception....root of the problem may be in initializeConfiguration() in Application.py
        self.clerk.director = self

        self.debug = self.inventory.debug
        
        return

    def _init(self):
        super(WebApplication, self)._init()

        """
        # Don't understand why I need this code
        # initialize table registry
        import vnf.dom
        vnf.dom.register_alltables()

        # set id generator for referenceset
        def _id():
            from vnf.components.misc import new_id
            return new_id(self)
        vnf.dom.set_idgenerator(_id)
        return
        """

    def _getPrivateDepositoryLocations(self):
        from os.path import join
        root = '..'
        content = join(root, 'content')
        config = join(root, 'config')

        return [config, content]

if __name__=='__main__':
    w = WebApplication(name = 'test')
    print w

__date__ = "$Jul 19, 2009 11:07:10 PM$"


