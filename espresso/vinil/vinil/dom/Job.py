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

class Job(Table):

    name = "job"
    import pyre.db
    
    id = pyre.db.varchar(name="id", length=10)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    type = pyre.db.varchar(name="type", length=100)
    type.meta['tip'] = "Type of the simulation (electron/phonon)"

    status = pyre.db.varchar(name="status", length=100)
    status.meta['tip'] = "Status of the job (not started/running/finished)"

    created = pyre.db.varchar(name="created", length=100)
    created.meta['tip'] = "Date and time when the simulation started"

    config = pyre.db.varchar(name="config", length=5000)
    config.meta['tip'] = "Configuration text"


def inittable(db):
    def job(id, type, status, created, config):
        r           = Job()
        r.id        = id
        r.type      = type
        r.status    = status
        r.created   = created
        r.config    = config
        return r
    records = [
        job( 1, 'pw', 'finished', '2009-06-11 13:35:34.081567', config_electons),
        job( 2, 'ph', 'finished', '2009-06-11 16:47:52.790506', config_phonons),
        job( 3, 'dos', 'finished', '2009-06-11 16:44:32.345214', config_pp)
        ]
    for r in records: db.insertRow( r )
    return



__date__ = "$Jul 29, 2009 8:31:54 PM$"


