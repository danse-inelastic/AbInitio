# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2008  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
Runs simulation for Quantum Espresso (QE)
"""

from Simulator import Simulator

class QE(Simulator):

    def render(self, computation, db=None, dds=None):
        self.dds = dds

        self._files = []
        self._makeScript(computation)

        return self._files


    def _makeScript(self, computation):
        engine = computation.engine
        handler = '_make_%s_script' % engine
        handler = getattr(self, handler)
        handler(computation)
        return


#    def _make_vasp_script(self, computation):
#
#        self._make_vasp_script1(computation)
#
#        job = computation.getJob(self.db)
#        np = job.numprocessors
#        cmds = [
#            '#!/usr/bin/env sh',
#            '. ~/.abinitio-env',
#            'chmod +x %s' % self.shscript1name,
#            'mpirun -np %d ./%s' % (np, self.shscript1name),
#            '',
#            ]
#        script = self.shscriptname
#        path = self._path(script)
#        open(path, 'w').write('\n'.join(cmds))
#
#        self._files.append(self.shscriptname)
#        return script
#
#
#    shscript1name = 'run1.sh'
#    def _make_vasp_script1(self, computation):
#        matter = computation.matter.dereference(self.db)
#        xyzfilename = self._makeXYZfile(matter)
#        xcFunctional = computation.xcFunctional
#        xcf = _xcf[xcFunctional]
#
#        params = [
#            ('name', matter.chemical_formula or computation.short_description or computation.id),
#            ('ecutoff', computation.kineticEnergyCutoff),
#            ('xcf', xcf),
#            ('mpmesh', ','.join([str(i) for i in computation.monkhorstPackMesh]) ),
#            ('unitcell', xyzfilename),
#            ('generateInputsOnly', computation.generateInputsOnly),
#            ]
#        cmds = [
#            '#!/usr/bin/env sh',
#            'vaspapp.py ' + ' '.join(['-%s=%s' % (k,v) for k,v in params]),
#            ]
#        shscript1name = self.shscript1name
#        path = self._path(shscript1name)
#        open(path, 'w').write('\n'.join(cmds))
#        self._files.append(shscript1name)
#
#        if not computation.generateInputsOnly:
#            # copy vasp input files if they alreday exist
#            dds = self.dds
#            files = [
#                'INCAR',
#                'POSCAR',
#                'POTCAR',
#                'KPOINTS',
#                ]
#            job = computation.getJob(self.db)
#            import shutil, os
#            for f in files:
#                src = dds.abspath(computation, f)
#                dest = dds.abspath(job, f)
#                debug.log('copying file %s to %s' % (src, dest))
#                if os.path.exists(src):
#                    shutil.copyfile(src, dest)
#                    self._files.append(f)
#                else:
#                    debug.log('file %s does not exist' % src)
#                continue
#
#        return
#
#
#    def _makeXYZfile(self, crystal):
#        from vnfb.components.matter_utils import crystal2xyz
#        contents = crystal2xyz(crystal)
#        filename = 'matter.xyz'
#        path = self._path(filename)
#        open(path, 'w').write('\n'.join(contents))
#        self._files.append(filename)
#        return filename
#
#_xcf = {
#    'PAW-PBE': 'pawpbe',
#    'PAW-GGA': 'pawgga',
#    'PAW-LDA': 'pawlda',
#    'USPP-GGA': 'usppgga',
#    'USPP-LDA': 'uspplda',
#    }





# End of file 
