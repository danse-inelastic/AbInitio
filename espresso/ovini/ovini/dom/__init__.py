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

maindom = 'ovini.dom'

def tables():
    # tables in this package
    tables = []

    # I have one table only :)
    m = '%s.%s' % (maindom, 'Job')
    module = _import(m)
    table = getattr(module, 'Job')
    tables.append( table )

    return tables

def _import(package):
    return __import__(package, {}, {}, [''])

__date__ = "$Jul 19, 2009 9:52:57 PM$"


