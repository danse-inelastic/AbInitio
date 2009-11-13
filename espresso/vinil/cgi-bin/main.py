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

from vinil.applications.WebApplication import WebApplication

def main():
    app =WebApplication(name='main')
    return app.run()

# main
if __name__ == '__main__':
    # invoke the application shell
    import journal
    import os
    debug = journal.debug('main' )
    debug.log(os.environ['PATH'] )
    debug.log(os.environ['PYTHONPATH'] )
    main()


# version
__id__ = "$Id$"

# End of file 
