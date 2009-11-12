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
Job Manager Daemon that serves as a client (jmclient) for job manager server (jmserver)
"""

from vinil.applications.JMClient import JMClient
import sys

app = JMClient(name='jmd')
app.run()
#sys.exit(app.run())

__date__ = "$Nov 12, 2009 12:32:59 PM$"


