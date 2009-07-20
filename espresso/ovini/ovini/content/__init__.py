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

def page(**kwds):
    from Page import Page
    return Page(**kwds)

def action(*args, **kwds):
    from Action import Action
    return Action( *args, **kwds )

def image(*args,**kwds):
    from Image import Image
    return Image(*args,**kwds)

def button(*args, **kwds):
    from Button import Button
    return Button(*args, **kwds)


__date__ = "$Jul 19, 2009 9:52:50 PM$"


