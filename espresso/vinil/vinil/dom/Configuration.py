#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

configPW = """ &control
    calculation='scf'
    restart_mode='from_scratch',
    tprnfor = .true.
    prefix='ni',
    pseudo_dir = '',
    outdir=''
 /
 &system
    ibrav=2,
    celldm(1) =6.65,
    nat=  1,
    ntyp= 1,
    nspin=2,
    starting_magnetization(1)=0.5,
    degauss=0.02,
    smearing='gauss',
    occupations='smearing',
    ecutwfc =27.0
    ecutrho =300.0
 /
 &electrons
    conv_thr =  1.0d-8
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
 Ni  26.98  Ni.pbe-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1"""

configPH = """phonons of Ni at gamma
&inputph
  tr2_ph=1.0d-16,
  prefix='ni',
  ldisp=.true.,
  nq1=2,
  nq2=2,
  nq3=2,
  amass(1)=58.6934,
  outdir='',
  fildyn='ni.dyn',
/"""

configDOS = """&inputpp
   prefix='ni',
   outdir='',
   fildos='',
   Emin=5.0,
   Emax=25.0,
   DeltaE=0.1
/"""

#outdir='/home/dexity/exports/vinil/content/temp/',
#fildos='/home/dexity/exports/vinil/output/ni.scf.dos.out',


from vinil.utils.utils import timestamp, newid, setname, ifelse
# Electron DOS (Ni_E_DOS): simulationId = 5

defaults    = (
               {"id": 1, "simulationId": 5, "type": "PW",
               "filename": "ni.scf.in", "text": configPW},
               {"id": 2, "simulationId": 6, "type": "PH",
                "filename": "ni.ph.in", "text": configPH},
               {"id": 3, "simulationId": 4, "type": "PP",
                "filename": "ni.pp.in", "text": configDOS},
                {"id": 4, "simulationId": 5, "type": "DOS",
                "filename": "ni.scf.dos.in", "text": configDOS},
              )


from pyre.db.Table import Table

class Configuration(Table):

    name = "configuration"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    simulationId    = pyre.db.varchar(name="simulationId", length=8)
    #simulationId.constraints = 'REFERENCES simulation (id)'    # Important
    simulationId.meta['tip'] = "simulationId"

    type        = pyre.db.varchar(name="type", length=1024, default='')
    type.meta['tip'] = "Type of configuration. Example: PW, PP"

    # Later on can be tranformed to a separate File table
    filename    = pyre.db.varchar(name="filename", length=1024, default='')
    filename.meta['tip'] = "Filename assiciated with this configuration"

    # To separate table?
    parser    = pyre.db.varchar(name="parser", length=1024, default='')
    parser.meta['tip'] = "Parser for configuration"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"

    timeCreated = pyre.db.varchar(name="timeCreated", length=16, default='')
    timeCreated.meta['tip'] = "timeCreated"

    timeModified = pyre.db.varchar(name="timeModified", length=16, default='')
    timeModified.meta['tip'] = "timeModified"

    text = pyre.db.varchar(name="text", length=8192, default='')
    text.meta['tip'] = "text"

    def setClerk(self, clerk):

    def updateRecord(self, director, params):
        """
        Updates configuration row (even if key in params is not present).
        'id' ans 'timeCreated' cannot be updated!
        """
        self.simulationId  = setname(params, self, 'simulationId')
        self.type          = setname(params, self, 'type')
        self.filename      = setname(params, self, 'filename')
        self.description   = setname(params, self, 'description')
        self.timeModified  = timestamp()
        self.text          = setname(params, self, 'text')
        
        director.clerk.updateRecord(self)   # Update record

    def createRecord(self, director, params):
        """Inserts configuration row """
        self.id            = ifelse(params.get('id'), params.get('id'), newid(director))
        self.simulationId  = setname(params, self, 'simulationId')
        self.type          = setname(params, self, 'type')
        self.filename      = setname(params, self, 'filename')
        self.description   = setname(params, self, 'description')
        self.timeCreated   = timestamp()
        self.timeModified  = timestamp()
        self.text          = setname(params, self, 'text')
        
        director.clerk.insertRecord(self)


    def deleteRecord(self, director):
        """Deletes record"""
        director.clerk.deleteRecord(self, id=self.id)

# For debugging
def inittable(db):

    for params in defaults:
        c   = Configuration()
        c.createRecord(db, params)

#    def configuration(params):
#        r               = Configuration()
#        r.id            = params['id']
#        r.simulationId  = params['simulationId']
#        r.type          = params['type']
#        r.filename      = params['filename']
#        r.timeCreated   = timestamp()
#        r.timeModified  = timestamp()
#        r.text          = params['text']
#        return r
#
#    records = []
#    for e in defaults:
#        records.append(configuration(e))
#
#    for r in records: db.insertRow( r )
    return

def test():
    for e in defaults:
        s = ""
        for v in e.keys():
            s += "%s: %s " % (v, e[v])
        print s

if __name__ == "__main__":
    test()

__date__ = "$Oct 5, 2009 8:58:32 AM$"


