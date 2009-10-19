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

class MultiTask(QETask):
    def __init__(self, setting, *tasks, cleanOutDir = False):
        QETask.__init__(self, setting, cleanOutDir)
        self.tasks = tasks
        self.cmdStr = tasks[0].cmdLine()
        self.name = tasks[0].name
        for task in tasks[1:]:
            self.cmdStr = self.cmdStr + ' ; ' + task.cmdLine()
            self.name = self.name + ' -> ' + self.cmdStr            
    
    def launch(self):
        for task in self.tasks:
            self.tasks[task].input.parse()
        self._run()
        for task in self.tasks:
            self.tasks[task].output.parse()