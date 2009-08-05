#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2005  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

def main():
    from luban.applications.WebApplication import WebApplication

    class MainApp(WebApplication):

        def __init__(self):
            WebApplication.__init__(self, name='main') #, asCGI=True)
            print "WebApplication"
            return

        def _getPrivateDepositoryLocations(self):

            import os
            odbroot = os.path.abspath('../content')
            config = os.path.abspath('../config')

            depos = []
            depos.append(config)

            for entry in os.listdir(odbroot):
                path = os.path.join(odbroot, entry)
                if os.path.isdir(path): depos.append(path)
                continue

            return depos


    app = MainApp()
    return app.run()


# main
if __name__ == '__main__':
    # invoke the application shell
#    import journal,os
#    debug = journal.debug('main' )
#    debug.log(os.environ['PATH'] )
#    debug.log(os.environ['PYTHONPATH'] )
    main()


# version
__id__ = "$Id$"

# End of file


__date__ = "$Aug 4, 2009 8:56:40 PM$"


