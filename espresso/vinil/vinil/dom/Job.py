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

from pyre.db.Table import Table

class Job(Table):
    # 'name' attribute should be present in every class table.
    name = "job"
    import pyre.db
    
    id          = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    userId      = pyre.db.varchar(name="userId", length=8)
    #userId.constraints = 'REFERENCES user (id)'
    userId.meta['tip'] = ""

    simulationId = pyre.db.varchar(name="simulationId", length=8)
    #simulationId.constraints = 'REFERENCES simulation (id)'
    simulationId.meta['tip'] = ""

    serverId    = pyre.db.varchar(name="serverId", length=8)
    #serverId.constraints = 'REFERENCES server (id)'
    serverId.meta['tip'] = ""

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = ""

    status = pyre.db.varchar(name="status", length=64, default='')
    status.meta['tip'] = ""

    timeSubmitted = pyre.db.varchar(name="timeSubmitted", length=16, default='')
    timeSubmitted.meta['tip'] = "timeSubmitted"

    timeStarted = pyre.db.varchar(name="timeStarted", length=16, default='')
    timeStarted.meta['tip'] = "timeStarted"

    timeRestarted = pyre.db.varchar(name="timeRestarted", length=16, default='')
    timeRestarted.meta['tip'] = "timeRestarted"

    timeCompleted = pyre.db.varchar(name="timeCompleted", length=16, default='')
    timeCompleted.meta['tip'] = "timeCompleted"

    exitCode = pyre.db.integer(name="exitCode", default=-1)
    exitCode.meta['tip'] = "exitCode"

    numberProcessors = pyre.db.integer(name="numberProcessors", default=0)
    numberProcessors.meta['tip'] = "numberProcessors"

    errorFilename = pyre.db.varchar(name="errorFilename", length=256, default='')
    errorFilename.meta['tip'] = "errorFilename"

    outputFilename = pyre.db.varchar(name="outputFilename", length=256, default='')
    outputFilename.meta['tip'] = "outputFilename"

    statusMessage = pyre.db.varchar(name="statusMessage", length=256, default='')
    statusMessage.meta['tip'] = "statusMessage"


# Generates only id's.
# Do I need generate default jobs? For debugging, yes!
"""
def inittable(db):
    def job(params):
        r           = Job()
        r.id        = params['id']
        r.userId        = params['userId']
        r.simulationId  = params['simulationId']
        r.serverId      = params['serverId']

        return r

    records = [
        job( {"id": 1, "userId": 1, "simulationId": 1, "serverId": 1}),
        job( {"id": 2, "userId": 1, "simulationId": 1, "serverId": 1}),
        job( {"id": 3, "userId": 1, "simulationId": 1, "serverId": 1} )
        ]
    for r in records: db.insertRow( r )
    return
"""


__date__ = "$Jul 29, 2009 8:31:54 PM$"

config_electons = """ &control
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

config_phonons = """phonons of Ni at gamma
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

config_pp = """&inputpp
   prefix='ni',
   outdir='',
   fildos='',
   Emin=5.0,
   Emax=25.0,
   DeltaE=0.1
/"""

