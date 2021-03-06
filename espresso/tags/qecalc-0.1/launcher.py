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
from setting import Setting
import qetorque
import os

class Launcher(Setting):
    def __init__(self, fname=None):
        Setting.__init__(self,fname)
        if self.useTorque:
            self.torque = qetorque.QETorque(self.configFileName)

    def _check(self, x):
        """Will check the exit status of QE program"""
        signal = x & 0xFF
        exitcode = (x >> 8) & 0xFF
        if exitcode != 0:
            raise Exception("Quantum Espresso crashed: check your settings and/or clean your 'outdir' directory")

    def _launch(self, cmdstr):
        if self.useTorque:
            self.torque.serial(cmdstr)
        else:
            self._check(os.system(cmdstr))

    def cleanOutDir(self):
        from parser.configParser import QEConfig
        import shutil
        qeConf = QEConfig(self.pwscfInput)
        qeConf.parse()
        outDir = qeConf.namelist('control').param('outdir')[1:-1]
        if self.useTorque:
            os.system('bpsh -a rm -r -f ' + outDir)
            os.system('bpsh -a mkdir ' + outDir)
        else:
            shutil.rmtree(outDir)
            os.mkdir(outDir)

    def pwscfLauncher(self):
        cmdstr = self.paraPrefix + " pw.x " +  self.paraPostfix + " -inp " + \
                 self.pwscfInput + " > " + self.pwscfOutput + "< /dev/null"
        print cmdstr         
        self._launch(cmdstr)

    def phLauncher(self):
        cmdstr_pw = self.paraPrefix + " pw.x " +  self.paraPostfix + " -inp " + \
                 self.pwscfInput + " > " + self.pwscfOutput + "< /dev/null"
        cmdstr_ph = self.paraPrefix + " ph.x " +  self.paraPostfix + " -inp " + \
                 self.phInput + " > " + self.phOutput + "< /dev/null"
        print cmdstr_ph        
        self._launch(cmdstr_pw + ' ; ' + cmdstr_ph)

    def singlePhononLauncher(self):
        #self.pwscfLauncher()
        self.phLauncher()

        cmdstr_dynmat = "dynmat.x < " + self.dynmatInput + " > " + self.dynmatOutput
        print cmdstr_dynmat
        self._check(os.system(cmdstr_dynmat))

    def matdynLauncher(self):
        """Execute matdyn.x after successful run of pw.x + ph.x + q2r.x"""
        cmdstr_matdyn = "matdyn.x -inp " + self.matdynInput + " > " + self.matdynOutput
        print cmdstr_matdyn
        self._check(os.system(cmdstr_matdyn))

    def multiPhononLauncher(self):
        """Runs complete sequence of programms needed to extract phonons except
        the last step: matdyn.x. Usecase: One then can regenerate matdyn.in
        for dispersions along different directions, phonon DOS etc"""
        #self.pwscfLauncher()
        self.phLauncher()
        cmdstr_q2r = "q2r.x < " + self.q2rInput + " > " + self.q2rOutput
        print cmdstr_q2r
        self._check(os.system(cmdstr_q2r))


    def multiPhononTaskLauncher(self):
        """Runs complete sequence of programms needed to extract phonons"""
        #self.pwscfLauncher()
        self.phLauncher()
        cmdstr_q2r = "q2r.x < " + self.q2rInput + " > " + self.q2rOutput
        print cmdstr_q2r
        self._check(os.system(cmdstr_q2r))
        cmdstr_matdyn = "matdyn.x < " + self.matdynInput + " > " + self.matdynOutput
        print cmdstr_matdyn
        self._check(os.system(cmdstr_matdyn))
