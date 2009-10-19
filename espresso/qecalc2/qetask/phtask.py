#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from qetask import QETask
from qeparser.qeoutput import QEOutput

class PHTask(QETask):
    def __init__(self, setting, cleanOutDir = False):
        QETask.__init__(self, setting, cleanOutDir)
        self.input = QEInput
        self.output = QEOutput
        self.cmdStr = self.setting.paraPrefix + " pw.x " +  \
                      self.setting.paraPostfix + " -inp " + \
                      self.setting.pwscfInput + " > " + \
                      self.setting.pwscfOutput + "< /dev/null"
        self.name = 'pw.x'