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

"""
Script that checks the status of job on the remote server.
Implemented for convenience
"""

def main():
    from jobmanager.applications.JobManager import JobManager

    app     = JobManager(name='main')
    return app.run()

if __name__ == "__main__":
    main()

__date__ = "$Oct 27, 2009 4:13:08 PM$"


