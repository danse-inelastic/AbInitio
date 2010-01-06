#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from qecalc.multiphononcalc import MultiPhononCalc

configString = """
# all the relevant input files must be preconfiguered for specific tasks
# before using this class

[Launcher]
# parallelization parameters
# if this section is empty - serial mode is used
paraPrefix:   mpiexec -n 8
paraPostfix: -npool 8

serialPrefix: mpiexec -n 1
serialPostfix:

#useTorque: True # default: False
#paraPrefix: mpirun --mca btl openib,sm,self
#paraPostfix: -npool 900
#
#serialPrefix: mpiexec
#serialPostfix:

#Name of a script to execute a command on multiple nodes
#relevant if outdir is not located on Parallel/Network File system.
#Default value is empty
#paraRemoteShell: bpsh -a

# this string will be passed to qsub, -d workingDir -V are already there:
paraTorqueParams: -l nodes=2:ppn=12 -N myjob -j oe
serialTorqueParams: -l nodes=1:ppn=1 -N myjob -j oe

outdir: /scratch/markovsk/d3paratest


[pw.x]
# pwscf input/output files
pwInput:  si.scf.in
pwOutput: si.scf.out


[ph.x]
#ph.x input/ouput, relevant to all phonon calculations:
phInput:  si.ph.in
phOutput: si.ph.out

[q2r.x]
# input/output files relevant to multiple phonon calculation
q2rInput:      q2r.in
q2rOutput:     q2r.out


[matdyn.x]
# input/output files relevant to multiple phonon calculation
matdynInput:   matdyn.in
matdynOutput:  matdyn.out

[d3.x]
d3Input:  si.d3.in
d3Output: si.d3.out

"""

qgrid = [2,2,2]

mphon = MultiPhononCalc(configString = configString)

mphon.ph.input.parse()
mphon.ph.input.namelist('inputph').remove('fildrho')
mphon.ph.input.qpoints.setAutomatic(qgrid)
mphon.ph.input.save()

mphon.launch([mphon.pwph, mphon.q2r])

#irred_qpoints = mphon.ph.output.property('qpoints')

grid, qpoints_indep, qpoints_full = mphon.ph.output.property('qpoints')

print qpoints_indep
print qpoints_full


#generate qpoint grid:


#[nq1,nq2,nq3], qpoints_indep, qpoints_full

#qpoints = kmesh.kMeshCart([2,2,2],mphon.pw.input.structure.lattice.reciprocalBase())

mphon.matdyn.input.parse()
mphon.matdyn.input.qpoints.set(qpoints_full)
mphon.matdyn.input.save()
mphon.matdyn.launch()

modes, freqs, qpoints =  mphon.matdyn.output.property('multi phonon')

print modes, freqs, qpoints