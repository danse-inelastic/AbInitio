from task import Task

class EcutConverger(Task):
    def __init__(self, fname=None):
        # Default values, see explanations below:        
        configDic = {
        'isNormConserving': 'True',
        'ecutInit': '32',
        'ecutStep': '4'
        }    
        Task.__init__(self, fname)
        
        self.isNormConserving = self.config.getboolean('EcutConverger','isNormConserving')
        self.ecutInit = self.config.getfloat('EcutConverger','ecutInit')
        self.ecutStep = self.config.getfloat('EcutConverger','ecutStep')
        self.ecutConverger()

    def ecutConverger(self):
        import writetopwscf
        
        if self.isNormConserving:
            ecutrhoMult = 4.
        else:
            ecutrhoMult = 8.
        ecutwfc = self.ecutInit
        runHistory = []
        for iStep in range(self.nMaxSteps):
            ecutrho = ecutrhoMult*ecutwfc
            writetopwscf.varnameValue(self.pwscfInput,"ecutwfc", ecutwfc)
            writetopwscf.varnameValue(self.pwscfInput,"ecutrho", ecutrho)
            self.getLauncher()
            runHistory.append( self.getEstimator() )
            if iStep >= 2:
                if self.isConverged(runHistory): break
            ecutwfc = ecutwfc + self.ecutStep

        print " optimized ecut value : ", ecutwfc
        print runHistory

