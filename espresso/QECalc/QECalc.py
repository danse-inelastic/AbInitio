from task import Task
from configParser import *
import numpy

class QECalc(Task):
    def __init__(self, fname):
        
        Task.__init__(self, fname)

        self.__qeConfig = QEConfig(self.pwscfInput)
        self.__qeConfig.parse()

#        self._kpts = self.getkPoints()
        self._ekincutoff = self.getEkincutoff()
        print self._ekincutoff

    def getNamelistParameter(self, namelist, parameter):
        self.__qeConfig.parse()
        return self.__qeConfig.namelist(namelist).param(parameter)

    def setNamelistParameter(self, namelist, parameter, value):
        self.__qeConfig.namelists[namelist].params[parameter] = str(value)
        self.__qeConfig.save(self.pwscfInput)

    def getCard(self, name):
        self.__qeConfig.parse()
        return self.__qeConfig.cards[name].getLines()


    def getkPoints(self):
        kpoints = []
        for line in self.getCard('k_points'):
            kpoints.append([float(k) for k in line.split()])
        return numpy.array(kpoints)

    def setkPointsAutomatic(self, kpoints):
        self.__qeConfig.cards['k_points'].removeLines()
        self.__qeConfig.cards['k_points'].setArgument('AUTOMATIC')
        string = ""
        for k in kpoints:
            string = string + str(k) + " "
        self.__qeConfig.cards['k_points'].addLine(string)
        self.__qeConfig.save(self.pwscfInput)

    def getEkincutoff(self):
        self.__qeConfig.parse()
        return float(self.__qeConfig.namelist('system').param('ecutwfc'))


qe = QECalc('config.ini')
qe.setNamelistParameter('system', 'ecutwfc', 44)
print qe.getkPoints().shape
qe.setkPointsAutomatic(numpy.array([17, 33, 12, 0, 0, 0]))