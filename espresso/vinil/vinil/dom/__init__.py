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

maindom = 'vinil.dom'


def tables():
    # tables in this package
    tables = []

    for t in tablenames:
        m = '%s.%s' % (maindom, t)
        module = _import(m)
        tables.append( getattr(module, t) )

    return tables

def _import(package):
    return __import__(package, {}, {}, [''])

__date__ = "$Aug 5, 2009 1:13:27 PM$"


