###########################################################
#### Perform several calculations on random supercells
###########################################################

def generateV3SiDistorted222Supercell():
    zerodisp = np.real(pols[0][0]) * 0.0
    refuc = applyPolarizationToUnitCell(uc, zerodisp)
    refsupercell = generateSupercell(refuc,supersize=[2,2,2])
    totaldisp = getUnitCellDisplacementVector(refsupercell)
    temperature = 300
    scalefactor = 1
    for qindex in range(1,4):
        for modeindex in range(24):
            disp = np.real(pols[qindex][modeindex])
            qvec = qpts[qindex]
            recipvectors = uc.getRecipVectors()
            qvec = np.dot(qvec, recipvectors)
            energy = np.sqrt(om2s[qindex][modeindex]) * 4.1357 # in meV
            ucpol = applyPolarizationToUnitCell(uc, disp)
            supercell = generateSupercell(ucpol,supersize=[2,2,2])
            randphi = np.random.rand() * 2*pi
            phi0 = randphi
            superpol = getUnitCellDisplacementVector(supercell)
            superdisp = modulatePolarization(supercell, superpol, qvec, phi0)
            sudisp = applyMassTemperature(supercell, superdisp, energy, temperature)
            sudisp *= scalefactor
            totaldisp += sudisp       
    supercell = generateSupercell(refuc,supersize=[2,2,2])
    applyAtomDisplacementsToUnitCell(supercell, totaldisp)
    return supercell


###

from AbInitio.AbiCalc.VaspCalc import VaspCalc
import shutil

Ncalcs = 1
energies = []

for n in range(Ncalcs):
    scell = generateV3SiDistorted222Supercell()
    print "RMS displacement:", rmsDispCubicMatch(supercell, refsupercell)
    vc = VaspCalc(unitcell=scell, kpts=[2,2,2], ekincutoff=320, name='V3Si_'+str(n), xc='pawpbe', vaspcmd='mr vasp')
    vc._vasp.incar['ISMEAR'] = -1
    vc._vasp.incar['SIGMA'] = 0.1
    vc._vasp.incar['EDIFF']=1e-5
    vc._vasp.incar['LREAL'] = '.FALSE.'
    vc._vasp.incar['PREC'] = 'Accurate'
    vc._vasp.incar['ALGO'] = 'Normal'
    vc._vasp.incar['LPLANE'] = '.TRUE.'
    vc._vasp.incar['NPAR'] = 4
    vc._vasp.incar['LSCALU'] = '.FALSE.'
    vc._vasp.incar['NSIM'] = 4
    e = vc.getPotEnergy()
    energies.append(e)
    shutil.copy('OUTCAR', 'OUTCAR_'+str(n))
    shutil.copy('CHGCAR', 'CHGCAR_'+str(n))
    shutil.copy('vasprun.xml', 'vasprun.xml_'+str(n))

print energies


