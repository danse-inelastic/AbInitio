import numpy
import os
from task.writetopwscf import varnameValue, atomic_positions
from task.task import Task

a_range = [ 5.67319, 5.67485, 5.67635, 5.67821, 5.67950, 5.68118, 5.68267, 5.68451, 5.68558, 5.68744, 5.68898, 5.69045, 5.69224, 5.69367, 5.69530, 5.69699, 5.69825, 5.69980, 5.70158, 5.70323, 5.70445, 5.70575, 5.70765, 5.70923, 5.71076, 5.71217, 5.71374, 5.71531, 5.71670]
c_a_range = [ 1.09067, 1.09084, 1.09110, 1.09116, 1.09154, 1.09169, 1.09196, 1.09202, 1.09253, 1.09257, 1.09281, 1.09307, 1.09316, 1.09345, 1.09363, 1.09378, 1.09416, 1.09438, 1.09447, 1.09463, 1.09504, 1.09540, 1.09542, 1.09561, 1.09584, 1.09613, 1.09634, 1.09653, 1.09684]

task = Task('config.ini')
task.set('multiple phonon')
indexRange = range(10,22,2)
print indexRange
for i in indexRange:
    varnameValue(task.pwscfInput, 'celldm(1)', a_range[i])
    varnameValue(task.pwscfInput, 'celldm(2)', a_range[i])
    varnameValue(task.pwscfInput, 'celldm(3)', c_a_range[i])
    geometry = ['Al       0.000000000   0.0000000000000000   0.000000000',
                'B        0.500000000   0.2886751345948129   '+str(c_a_range[i]/2.),
                'B        0.000000000   0.5773502691896257   '+str(c_a_range[i]/2.)]
    atomic_positions(task.pwscfInput, geometry)
    task.getLauncher()
    Pol, Omega, qPoints = task.getEstimator()
    polName = "polarizations"+str(i)+".txt"
    omegaName = "freqs"+str(i)+".txt"
    qPointsName = "qpoints"+str(i)+".txt"
    cpMatdynModesCmdStr = "cp "+"matdyn.modes "+ "matdyn" + str(i) + ".modes"
    os.system(cpMatdynModesCmdStr)
#   save pol vectors, frequencies and qpoints:  
    #numpy.savetxt(polName, Pol)
#    numpy.savetxt(omegaName, Omega)
#    numpy.savetxt(qPointsName, qPoints)
#   backup the current state
    backupName = 'backup'+str(i)+'.tgz'
    backupCmdStr = "tar -zcf" + backupName + " *  --exclude='*tgz' --no-recursion"
    os.system(backupCmdStr)
    
