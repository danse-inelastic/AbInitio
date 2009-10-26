# -*- Python -*-
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



"""
a job is something that can be delivered to computation server and be launched.

It needs to contain:

  - a bash script to be launched on the computation server
  - all necessary files that are needed for the bash script to run, such as
    * configuration files
    * data files

History: Adapted from JobBuilder.py
"""

import os
from pyre.components.Component import Component

class Simulator(Component):

    shscriptname = 'run.sh'

    def __init__(self, name, path):
        super(Simulator, self).__init__(name)
        
        self.path = path
        self.initRootDir()
        return


    def build(self, computation, db=None, dds=None):
        self.dds = dds
        return self.render(computation, db=db, dds=dds)
    

    def render(self, computation, db=None, dds=None):
        raise NotImplementedError


    def _path(self, filename):
        return os.path.join(self.path, filename)


    def _initRootDir(self):
        if not os.path.exists(self.path):
            try:
                os.makedirs(self.path)
            except Exception, err:
                raise RuntimeError, "unable to create directory %r. %s: %s" % (
                    self.path, err.__class__.__name__, err)
        return


#    def registerDependency(self, dependency):
#        '''register a dependency to my dependency list
#
#        A dependency is a db record. For example, a job requires a density of
#        states curve, and that curve is stored as db record DOS(id=ABCDE).
#        The data files for the db record DOS(id=ABCDE) need to be made
#        available at the computation server.
#        '''
#
#
#        computation = self.computation
#        job = computation.getJob(db)
#        dependencies = job.dependencies
#
#        dependencies.add(dependency, db=db)
#        return


#    def getDependencies(self):
#        db = self.db
#
#        computation = self.computation
#        job = computation.getJob(db)
#        dependencies = job.dependencies
#
#        return [value for key, value in dependencies.dereference(db)]
    

    
    pass # end of JobBuilder


# version
__id__ = "$Id$"

# End of file 
