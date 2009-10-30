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
Script that periodically checks the status of job on the remote server.
Implemented for convenience
"""

def main():
    from jobmanager.applications.JobManager import JobManager

    app     = JobManager(name="jm", action="trace")
    return app.run()

if __name__ == "__main__":
    main()

__date__ = "$Oct 29, 2009 1:32:57 PM$"


