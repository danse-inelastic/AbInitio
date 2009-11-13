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

def main():
    from vinil.jmd.JMDaemon import JMDaemon

    app = JMDaemon(name='jmd')
    return app.run(spawn=True)

if __name__ == "__main__":
    main()


__date__ = "$Nov 12, 2009 12:32:59 PM$"


