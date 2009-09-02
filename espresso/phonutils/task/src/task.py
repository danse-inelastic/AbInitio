from estimator import Estimator
from launcher import Launcher

class Task(Estimator, Launcher):
    def __init__(self, fname=None):
        # Default values, see explanations below:        
        taskDic = {
        'taskName': 'total energy',
        'tolerance': '1',
        'nMaxSteps': '10'
        }
        Estimator.__init__(self,fname)
        Launcher.__init__(self,fname)                
        # value to converge with respect to k-points or energy cutoffs
        # currently can be 'total energy', 'single phonon', or 'geometry':
        self.taskName = self.config.get('Task', 'taskName')

        # convergence criteria in percents:        
        self.tolerance = self.config.getfloat('Task','tolerance')
        
        # maximum number of optimization steps:
        self.nMaxSteps = self.config.getint('Task','nMaxSteps')  
        
        self.lookupTable = {
        'total energy' : (self.pwscfLauncher, self.getTotalEnergy),
        'single phonon': (self.singlePhononLauncher, self.getSinglePhonon),
        'geometry'     : (self.pwscfLauncher, self.getLatticeParameters),
        'multiple phonon': (self.multiPhononLauncher, self.getMultiPhonon)
        }        
        
        assert self.lookupTable.has_key(self.taskName), "Convergence \
        estimator's name is not known"
        
    def set(self,taskName):
        self.taskName = taskName
    
    def get(self):
        return self.taskName
    
    def isConverged(self,runHistory):
        import math
#	    'Check for convergence:  two last runs should be less than the tolerance value'
#        if there is a list of values, the code will choose one with maximum error
        tol1 = []
        tol2 = []
        valTol = 1e-7
        for i in range( len(runHistory[-1]) ):
            # check if the denominator is not zerro:
            if math.fabs( runHistory[-2][i] ) > valTol and \
               math.fabs( runHistory[-2][i] ) > valTol :
                tol1.append( math.fabs( runHistory[-1][i]/runHistory[-2][i] - 1.0) )
                tol2.append( math.fabs( runHistory[-2][i]/runHistory[-3][i] - 1.0) )
        if max(tol1) < self.tolerance/100. and max(tol2) < self.tolerance/100.:
            print "\nSuccess! ",self.taskName," estimator value in two \
            consecutive runs\ndiffers less than ",self.tolerance," percent: ", max(tol2)*100, max(tol1)*100
            return True
        else:
            print runHistory[-1]
            return False

    def getLauncher(self):          
        return self.lookupTable[self.taskName][0]()
    
    def getEstimator(self):
        return self.lookupTable[self.taskName][1]()
                   
