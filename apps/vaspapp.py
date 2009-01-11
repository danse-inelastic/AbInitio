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


from pyre.applications.Script import Script


class VaspApp(Script):


    class Inventory(Script.Inventory):

        import pyre.inventory

        name = pyre.inventory.str('name', default='')
        name.meta['tip'] = 'name of this computation'

        unitcell = pyre.inventory.str('unitcell', default='')

        ecutoff = pyre.inventory.float('ecutoff', default=0)

        xcf = pyre.inventory.str('xcf', default='PAW-PBE')
        xcf.meta['tip'] = 'exchange correlation functional'

        mpmesh = pyre.inventory.str('mpmesh', default=[1,1,1], validator=intvec3)
        mpmesh.meta['tip'] = 'Monkhorst-Pack mesh'

        generateInputsOnly = pyre.inventory.bool('generateInputsOnly', default=0)
        generateInputsOnly.meta['tip'] = 'If true, only generate vasp input files'

        cmd = pyre.inventory.str('cmd', default='vasp')
        

    def main(self, *args, **kwds):
        vc = self._createVaspCalculator()

        if self.generateInputsOnly:
            vc._vasp._makeVaspInputs()
        else:
            open('out.energy', 'w').write( '%s' % (vc.getPotEnergy(),) )
            open('out.forces', 'w').write( '%s' % (vc.getForces(),) )
            open('out.stress', 'w').write( '%s' % (vc.getStress(),) )
        return


    def _createVaspCalculator(self):
        unitcell = self.unitcell
        ecutoff = self.ecutoff
        xcf = self.xcf
        mpmesh = self.mpmesh
        name = self.name
        cmd = self.cmd
        
        from AbInitio.AbiCalc import VaspCalc as VaspCalc
        return VaspCalc.VaspCalc(
            unitcell = unitcell,
            kpts = mpmesh,
            ekincutoff=ecutoff,
            name=name,
            vaspcmd=cmd)


    def _readUnitCell(self, filepath):
        from sampleassembly.crystal.ioutils.xyzfile import read
        return read(filepath)


    def __init__(self):
        Script.__init__(self, 'vaspapp')
        return


    def _configure(self):
        Script._configure(self)
        self.name = self.inventory.name
        self.unitcell_path = self.inventory.unitcell
        self.ecutoff = self.inventory.ecutoff
        self.xcf  = self.inventory.xcf
        self.mpmesh = self.inventory.mpmesh
        self.generateInputsOnly = self.inventory.generateInputsOnly
        self.cmd = self.inventory.cmd
        return


    def _init(self):
        Script._init(self)
        if not self._showHelpOnly: self.unitcell = self._readUnitCell(self.unitcell_path)
        return



def main():
    app = VaspApp()
    return app.run()


# main
if __name__ == '__main__':
    # invoke the application shell
    main()


# version
__id__ = "$Id$"

# End of file 
