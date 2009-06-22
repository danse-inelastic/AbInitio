#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                        (C) 2008  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

def intvec3(l):
    if isinstance(l, str):
        l = eval(l)
    l = list(l)
    assert len(l)==3
    for i in l: assert isinstance(i, int) or isinstance(i, long)
    return l


def floatvec2(l):
    if isinstance(l, str):
        l = eval(l)
    l = list(l)
    assert len(l)==2
    for i in l: assert isinstance(i, float) or isinstance(i, int)
    return l



from pyre.applications.Script import Script


class PhonApp(Script):


    class Inventory(Script.Inventory):

        import pyre.inventory

        name = pyre.inventory.str('name', default='')
        name.meta['tip'] = 'name of this computation'

        unitcell = pyre.inventory.str('unitcell', default='')

        dosaxis = pyre.inventory.str('dosaxis', default=(0,20), validator=floatvec2)
        
        amplitude = pyre.inventory.float('amplitude', default=0.1)

        qgridsize = pyre.inventory.str('qgridsize', default=[10,10,10], validator=intvec3)
        qgridsize.meta['tip'] = 'Q grid size'

        qgridhasgamma = pyre.inventory.bool('qgridhasgamma', default=True)
        fullqgrid = pyre.inventory.bool('fullqgrid', default=True)

        supersize = pyre.inventory.str('supersize', default=[2,2,2], validator=intvec3)
        supersize.meta['tip'] = 'super cell size'


    def main(self, *args, **kwds):
        #self._callVASP()
        #return
    
        pc = self._createPhonCalculator()
        #pc.genSupercell(self.supersize)
        pc.genPhonSupercell()
        pc._makePosFiles()

        pc.calcForces()
        pc.calcPhonons()

        # parse phon outputs and convert them to idf output files
        from AbInitio.AbiPhon.parsing.phonParsers import parsePhon2IDF, parseDOS_meV2IDF
        if self.qgridhasgamma and self.fullqgrid:
            # in this case we can generate Qgridinfo, and padding q grid
            padding = True
            qgridsize = self.qgridsize
            # write Qgridinfo
            self._writeQgridinfo()

        else:
            padding = False
            qgridsize = None
            
        parsePhon2IDF(padding=padding, qgridsize=qgridsize)

        # convert DOS.meV to idf format
        parseDOS_meV2IDF()

        return


    def _writeQgridinfo(self):
        f = open('Qgridinfo', 'w')
        unitcell = self.unitcell
        bs = unitcell.getRecipVectors()
        ns = self.qgridsize
        for i in range(3):
            b = bs[i]
            b = tuple(b)
            n = ns[i]+1
            f.write('n%d=%d; b%d=%s\n' % (i+1, n, i+1, b) )
            continue
        return
    
    
    def _callVASP(self):
        from crystal.crystalIO import converters
        from AbInitio.vasp.vasp import VASP
        v = VASP(pw=268, kpts=(2,2,2), xc='pawpbe', name=self._name, vaspcmd='vasp')

        uc = self.unitcell

        loa = converters.unitCell2ListOfAtom(uc)
        v._SetListOfAtoms(loa)
        v.atoms()
        v.incar['EDIFF']=1.0e-5
        v.makePoscarFile()
        v.makePotcarFile()

        v.GetPotentialEnergy()
        return


    def _createPhonCalculator(self):
        unitcell = self.unitcell
        name = self._name
        supersize = self.supersize
        qgridsize = self.qgridsize
        dosmin, dosmax = self.dosaxis
        amplitude = self.amplitude
        
        from AbInitio.AbiPhon.PhonCalc import PhonCalc
        return PhonCalc(
            unitcell,
            name=name,
            supersize=supersize,
            qgridsize=qgridsize,
            dosmin=dosmin, dosmax=dosmax,
            amplitude=amplitude,
            LGAMMA = '.'+str(self.qgridhasgamma).upper()+'.',
            LFULLQGRID = '.'+str(self.fullqgrid).upper()+'.',
            )


    def _readUnitCell(self, filepath):
        from sampleassembly.crystal.ioutils.xyzfile import read
        return read(filepath)


    def __init__(self):
        Script.__init__(self, 'phonapp')
        return


    def _configure(self):
        Script._configure(self)
        self._name = self.inventory.name
        self.unitcell_path = self.inventory.unitcell
        self.dosaxis = self.inventory.dosaxis
        self.qgridsize = self.inventory.qgridsize
        self.qgridhasgamma = self.inventory.qgridhasgamma
        self.fullqgrid = self.inventory.fullqgrid
        self.supersize = self.inventory.supersize
        self.amplitude = self.inventory.amplitude
        return


    def _init(self):
        Script._init(self)
        if not self._showHelpOnly:
            unitcell = self._readUnitCell(self.unitcell_path)
            from crystal.UnitCell import positive_volume_unitcell
            unitcell = positive_volume_unitcell(unitcell)
            self.unitcell = unitcell
        return



def main():
    app = PhonApp()
    return app.run()


# main
if __name__ == '__main__':
    # invoke the application shell
    main()


# version
__id__ = "$Id$"

# End of file 
