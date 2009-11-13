#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                 Jiao Lin
#                      California Institute of Technology
#                        (C) 2008  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from luban.applications.utils import redirectWarningsToJournal
redirectWarningsToJournal('warning')

from vinil.application.SuperApp import SuperApp

def main():
    app = SuperApp('main-superapp')
    return app.run()


# main
if __name__ == '__main__':
    # invoke the application shell
    import journal
    debug = journal.debug('main' )
    debug.log(os.environ['PATH'] )
    debug.log(os.environ['PYTHONPATH'] )
    main()


# version
__id__ = "$Id$"

# End of file 
