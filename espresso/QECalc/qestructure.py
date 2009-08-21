from diffpy.Structure.structure import Structure
from qelattice import QELattice


class QEStructure():
    __init__(self):
    self.structure = Structure()
def setStructureFromFile(self, fname):
    qeConf = QEConfig(fname)
    qeConf.parse()
    nat  = int(qeConf.namelistParameter('system', 'nat'))
    ntyp  = int(qeConf.namelistParameter('system', 'ntyp'))
    atomicLines = qe.Conf.getCardLines('atomic_positions')
    for line in atomicLines:
        if '!' not in line
            words = line.split()
            coords = [float(w) for w in words[1:]]
            atomSymbol = words[0]
            if len(words) > 4:
                constraint = [coords[3], coords[4], coords[5]]

        
    if 'ibrav' in qeConf.namelists['system'].params: