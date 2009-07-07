from config import Config

class Launcher(Config):
    def __init__(self, fname=None):
        Config.__init__(self,fname)

    def pwscfLauncher(self):
        import os
        cmdstr = self.paraPrefix + " pw.x " +  self.paraPostfix + " -inp " + \
                 self.pwscfInput + " > " + self.pwscfOutput + "< /dev/null"
        print cmdstr         
        os.system(cmdstr)

    def singlePhononLauncher(self):
        import os
        self.pwscfLauncher()
        cmdstr_ph = self.paraPrefix + " ph.x " +  self.paraPostfix + " -inp " + \
                 self.phInput + " > " + self.phOutput + "< /dev/null"
        print cmdstr_ph        
        os.system(cmdstr_ph)
        cmdstr_dynmat = "dynmat.x < " + self.dynmatInput + " > " + self.dynmatOutput
        print cmdstr_dynmat
        os.system(cmdstr_dynmat)
