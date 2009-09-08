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

from luban.applications.WebApplication import WebApplication as Base

class WebApplication(Base):

    class Inventory(Base.Inventory):

        import opal.inventory
        import pyre.inventory

        # components
        actor = opal.inventory.actor(default='greet')
        actor.meta['tip'] = "the component that defines the application behavior"

        imagepath       = pyre.inventory.str(name='imagepath', default = 'images' )
        javascriptpath  = pyre.inventory.str(name='javascriptpath', default = 'javascripts' )

    def __init__(self, name):
        Base.__init__(self, name)

    def main(self, *args, **kwds):
        actor = self.actor

        # If no actor is provided e.g.: example.com/main.cgi
        if actor is None:
            self.actor = self.retrieveActor('greet')#render('greet') #

        try:
            page = self.actor.perform(self, routine=self.inventory.routine, debug=self.debug)

            if isinstance(page, basestring):
                print page
            else:
                self.render(page)
        except:
            self.plainBugReport()

    def retrievePage(self, name):
        page = super(WebApplication, self).retrievePage(name)#.retrieveVisual(name) #
        if page:
            return page
        raise RuntimeError, "Unable to load page %s" % name

    def plainBugReport(self):
        print '<pre>'
        import traceback
        traceback.print_exc()
        print '</pre>'
        return

    def render(self, page=None):
        self.weaver.resetRenderer()
        self.weaver.weave(document=page, stream=self.stream)
        return

    def _defaults(self):
        super(WebApplication, self)._defaults()
        from luban.components.weaver.web import weaver
        self.inventory.weaver = weaver()
        return

    def _configure(self):
        super(WebApplication, self)._configure()

        # I don't know if I need this line?
        #self.clerk.director = self

        return

    def _init(self):
        super(WebApplication, self)._init()
            
    def _getPrivateDepositoryLocations(self):
        from os.path import join
        root = '..'
        content = join(root, 'content')
        config = join(root, 'config')

        return [config, content]

if __name__=='__main__':
    w = WebApplication(name = 'test')
    print w


__date__ = "$Aug 5, 2009 4:23:32 PM$"


