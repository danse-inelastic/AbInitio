from task import Task
from string import split

class KConverger(Task):
    def __init__(self, fname=None):
        Task.__init__(self, fname)
        
        self.isMetallic =  self.config.getboolean('KConverger','isMetallic')
        self.kInit = [ int(k) for k in split(self.config.get('KConverger','kInit')) ]
        self.kStep = [ int(k) for k in split(self.config.get('KConverger','kStep')) ]
        
        self.kConverger()
        
    def kConverger(self):
        import writetopwscf   
        
        k_points = self.kInit
        runHistory = []
        for iK in range(self.nMaxSteps):
            writetopwscf.k_points(self.pwscfInput,k_points)
            self.getLauncher()
            runHistory.append( self.getEstimator() )
            if iK >= 2:
                if self.isConverged(runHistory): break            
            for i in range(len(self.kStep)):
                k_points[i] = k_points[i] + self.kStep[i]            
        print "optimized kpoints : ", k_points
        print runHistory
        return

