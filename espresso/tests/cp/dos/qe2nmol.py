
from MMTK import *
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
import string


SPECIES     = [['N',1,[('n','N')],[]], ['H',3,[('h','H')],[]]]
VELFILE     = 'nh3_mol.vel'
POSFILE     = 'nh3_mol.pos'
NCFILE      = 'nh3_mol.nc'
N_VECTORS   = 2
PBC         = 3
BOX_X       = 12
BOX_Y       = 12
BOX_Z       = 12
DT          = 0.00120944
N_ATOMS     = 4

class QEData:

    def __init__(self):
        self.title      = "Generated from CP of Quantum Espresso"
        self.velfile    = open(VELFILE)
        self.posfile    = open(POSFILE)

        self.makeUniverse(PBC, SPECIES)
        self.universe.initializeVelocitiesToTemperature(0.)


    def makeUniverse(self, pbc, molecules):
        if pbc == 0:
            self.universe = InfiniteUniverse()
        else:
            self.universe = OrthorhombicPeriodicUniverse((0., 0., 0.))
        number = 0
        for mol_name, mol_count, atoms, constraints in molecules:
            for i in range(mol_count):
                atom_objects = []
                for element, name in atoms:
                    a = Atom(element, name = name)
                    a.number = number
                    number = number + 1
                    atom_objects.append(a)
                if len(atom_objects) == 1:
                    self.universe.addObject(atom_objects[0])
                else:
                    ac = AtomCluster(atom_objects, name = mol_name)
                    for i1, i2, d in constraints:
                        ac.addDistanceConstraint(atom_objects[i1],
                                                 atom_objects[i2],
                                                 d)
                    self.universe.addObject(ac)
        self.universe.configuration()


    def writeTrajectory(self, trajectory_name, block_size=1):
        trajectory  = Trajectory(self.universe, trajectory_name, 'w',
                                self.title, block_size=block_size)
        actions     = [TrajectoryOutput(trajectory, ["all"], 0, None, 1)]
        snapshot    = SnapshotGenerator(self.universe, actions=actions)
        conf        = self.universe.configuration()
        vel         = self.universe.velocities()
        grad        = ParticleVector(self.universe)
        nvectors    = N_VECTORS
        natoms      = N_ATOMS
        self._setSize()

        try:
            while True:
                vline   = self.velfile.readline()
                pline   = self.posfile.readline()
                if not vline:   
                    break

                vdata   = vline.split()
                pdata   = pline.split()

                if len(pdata) == 2: # Example: ["10", "0.00120944"]
                    step        = int(vdata[0])
                    step_data   = {'time': step*DT}

                for i in range(natoms):
                    conf.array[i]   = map(float, string.split(self.posfile.readline()))
                    vel.array[i]    = map(float, string.split(self.velfile.readline()))

                conf.array = Units.Ang*conf.array
                if nvectors > 1:
                    vel.array = Units.Ang/Units.ps * vel.array

                snapshot(data=step_data)
        finally:
            trajectory.close()


    def _setSize(self):
        "Set universe size"
        box_x = BOX_X*Units.Ang
        box_y = BOX_Y*Units.Ang
        box_z = BOX_Z*Units.Ang
        self.universe.setSize((box_x, box_y, box_z))



if __name__ == "__main__":
    data = QEData()
    data.writeTrajectory(NCFILE)
