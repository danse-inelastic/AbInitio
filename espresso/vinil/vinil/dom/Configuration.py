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

configElectons = """ &control
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

configPhonons = """phonons of Ni at gamma
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

configPP = """&inputpp
   prefix='ni',
   outdir='',
   fildos='',
   Emin=5.0,
   Emax=25.0,
   DeltaE=0.1
/"""


from pyre.db.Table import Table

class Configuration(Table):

    name = "configuration"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    simulationId    = pyre.db.varchar(name="simulationId", length=8)
    #simulationId.constraints = 'REFERENCES simulation (id)'
    simulationId.meta['tip'] = "simulationId"

    cname = pyre.db.varchar(name="cname", length=1024, default='')
    cname.meta['tip'] = "cname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"

    date = pyre.db.varchar(name="date", length=16, default='')
    date.meta['tip'] = "date"

    text = pyre.db.varchar(name="text", length=8192, default='')
    text.meta['tip'] = "text"


def timestamp():
    import time
    return int(time.time())

# For debugging
def inittable(db):
    
    def configuration(params):
        r           = Configuration()
        r.id        = params['id']
        r.simulationId  = params['simulationId']
        r.cname     = params['cname']
        r.description   = params['description']
        r.date      = params['date']
        r.text      = params['text']
        return r

    records = [
        configuration( {"id": 1, "simulationId": 1, "cname": "Electron Simulation",
                        "description": "", "date": timestamp(),
                        "text": configElectons} ),
        configuration( {"id": 2, "simulationId": 2, "cname": "Phonon Simulation",
                        "description": "", "date": timestamp(),
                        "text": configPhonons} ),
        configuration( {"id": 3, "simulationId": 3, "cname": "Post Processing",
                        "description": "", "date": timestamp(),
                        "text": configPP} )        
        ]
    for r in records: db.insertRow( r )
    return

__date__ = "$Oct 5, 2009 8:58:32 AM$"


