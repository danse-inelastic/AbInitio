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

from luban.content.Splitter import Splitter
from luban.content.Paragraph import Paragraph
from luban.content import load
from luban.content.Link import Link


class SimServer:

    def __init__(self, director):
        self._director  = director
        self._clerk     = director.clerk
        self._type      = "settings"

    def getServer(self, id):      # simulation id
        settings  = self._clerk.getConfigurations(where="simulationId='%s' AND type='%s'" % (id, self._type))

        if self._serverIsSet(settings):
            text    = Paragraph(text="Exists")
#                        Link(label=self._label(settings), Class="action-link",
#                        onclick=load(actor="espresso/settings-view",
#                        routine="link", id=id))
        else:
            text    = Paragraph(text="None")

        return text

    def _serverIsSet(self, settings):
        """Checks if server is set"""
#        if not settings:
#            return False

        import ConfigParser
        import StringIO

        config  = """
[server]
server-name = foxtrot.danse.us
"""
        if config:  # Implies that it has sections already
            fp  = StringIO.StringIO(config)

            parser  = ConfigParser.ConfigParser()
            parser.readfp(fp)
            if self._serverName(parser.get("server", "server-name")):
                return True

        return False
                
    def _serverName(self, name):
        if name == '':
            return False

        # Fix it
#        if self._director:
#            servers  = self._clerk.getServers()
#            for s in servers:
#                if name == s.sname:
#                    return True

        return True
#        return False

        
#        s   = settings[0]
#        if s:
#            config  = s.text
            
        



    def _label(self, settings):
        """Returns filename"""
        if settings:
            return settings[0].filename

        return "Add"


    def _getActor(self, settings):
        """Returns proper actor depending if 'input' exists"""
        if settings:   # View
            return "espresso/settings-view"

        return "espresso/settings-add" # Create New




if __name__ == "__main__":
    params   = SimServer(None)
    #params.serverIsSet(None)


__date__ = "$Nov 11, 2009 1:21:52 PM$"


