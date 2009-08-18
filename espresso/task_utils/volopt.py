from task.writetopwscf import varnameValue, atomic_positions
from task.task import Task
import numpy as np
from scipy.optimize import brent
import scipy

def getHexEnergy(c, *args ):
    """ Total energy launcher for scipy "brent" routine
        'args' is a tuple with volume and task name \
        Volume value should be directly related to lattice constants: \
        E.g.: Vhex = a^2*c ommiting all the constant factors \
        task should be set to 'total energy ' """
    volume = args[0]
    task = args[1]
    a = np.sqrt(volume/c)
    geometry = ['Al       0.000000000   0.0000000000000000   0.000000000',
                'B        0.500000000   0.2886751345948129   '+str(c/a/2.),
                'B        0.000000000   0.5773502691896257   '+str(c/a/2.)]
    varnameValue(task.pwscfInput, 'celldm(1)', a)
    varnameValue(task.pwscfInput, 'celldm(2)', a)
    varnameValue(task.pwscfInput, 'celldm(3)', c/a)
    atomic_positions(task.pwscfInput, geometry)
    
    task.getLauncher()
    return task.getEstimator()[0]

# will find optimal a and c of hexagonal lattice of fixed volume:
def hexVolOpt(a0, c0_a0, volumeExpansion):
    '''provide initial(equilibrium) lattice parameters a0, c0/a0, and \
    desired volume expansion in percents at which one should run \
    the optimization '''
    task = Task('config.ini')
    task.set('total energy')
    print task.get()
#   Initial(equilibrium) volume:    
    c0 = a0*c0_a0    
    volume = a0*a0*c0
            
    volume = volume + volume*volumeExpansion/100.
#   initial assumption: all latice parameters expand equally in percents    
    cExpansion = (1.+volumeExpansion/100.)**(1./3.)
    c = c0*cExpansion
    
    prcntC = 0.2 # percent of c we want to bracket around it for minima search(does not have to warantee the minima is inside)s

    brentOut = brent(getHexEnergy, (volume, task), (c-c*prcntC/100, c+c*prcntC/100), tol = 1.e-5, full_output = 1)    
    print brentOut
    c = brentOut[0]
    energy = brentOut[1]
    a = np.sqrt(volume/c)
    return a, c/a, energy

#    stepPrcntC = 0.1 # percent
#    cRange = np.r_[c-c*prcntC/100:c+c*prcntC/100:c*stepPrcntC/100]
#    aRange = np.zeros(len(cRange))
#    aRange = aRange + volume 
#    aRange = np.sqrt(aRange/cRange)
#    print cRange
#    print aRange
#    energies = []    
#    for a, c in zip(aRange, cRange):
#        varnameValue(task.pwscfInput, 'celldm(1)', a)
#        varnameValue(task.pwscfInput, 'celldm(2)', a)
#        varnameValue(task.pwscfInput, 'celldm(3)', c/a)
#        task.getLauncher()
#        energies.append(task.getEstimator())
#    file = open('energies'+strvolumeExpansion)+'.txt','w')
#    for a,c, e in zip(aRange, cRange, energies):
#        print a, c, e
        
if __name__ == '__main__':

    volPercRange = scipy.linspace(0.1, 3.0, 29)


    for volPrcnt in volPercRange:
        print hexVolOpt(5.6717525, 1.09041, volPrcnt)
