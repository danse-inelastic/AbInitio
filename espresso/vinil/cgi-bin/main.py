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

def main():

    from vinil.applications.WebApplication import WebApplication

    app = WebApplication(name='main')
    return app.run()

# main
if __name__ == '__main__':
    # invoke the application shell
    try:
        main()
    except:
        import traceback
        import time
        t = time.ctime()
        messages = traceback.format_exc()
        print messages



# version
__id__ = "$Id$"

# End of file 
