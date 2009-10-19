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

class PWPHTask(QETask):
    def __init__(self, pwTask, phTask, cleanOutDir = False):
        QETask.__init__(self, pwTask.setting, cleanOutDir)
        self.tasks = [pwTask, phTask]
        self.cmdStr = pwTask.cmdLine() + ' ; ' + phTask.cmdLine()
        self.name = 'pw.x -> ph.x'
    
    def launch(self):
        for task in self.tasks:
            self.tasks[task].input.parse()
        self._run()
        for task in self.tasks:
            self.tasks[task].output.parse()