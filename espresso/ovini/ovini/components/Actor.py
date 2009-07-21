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

from opal.components.Actor import Actor as base

class Actor(base):

    def nyi(self, director):
        """notify the user that the requested routine is not yet implemented"""
        page = director.retrievePage("nyi")
        main = page._body._content._main
        document = main.document(title = 'Under construction...')
        p = document.paragraph()
        p.text = [
            "Not implemented yet! actor=%s, routine=%s" % (
            self.name, director.inventory.routine),
            ]
        return page


    
__date__ = "$Jul 20, 2009 10:13:32 AM$"


