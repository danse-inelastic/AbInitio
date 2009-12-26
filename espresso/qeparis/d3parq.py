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
from mpi4py import MPI
import sys
from qecalc.d3calc import D3Calc
from qeutils import kmesh

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

[ph.x multi]
#ph.x input/ouput, relevant to all phonon calculations:
phInput:  ph_disp.in
phOutput: ph_disp.out

[dynmat.x]
#dynmat.x input/output files relevant to single phonon calculation
dynmatInput:  dynmat.in
dynmatOutput: dyn.out


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

#def D3ParQ(qpoints, d3):
#    outDir = d3.setting.outDir
#    myrank = MPI.COMM_WORLD.Get_rank()
#
#
#
#if __name__ == "__main__":
#
#    if len(sys.argv) != 2:
#        raise
#
#    configName = sys.argv[1]
#    d3 = D3Task(configName)
#    outDir = d3.setting.outDir
#    myrank = MPI.COMM_WORLD.Get_rank()
#    outDir = outDir + '/d3task_' + str(myrank)
#
#    d3.setting.outdir = outdir
#


if __name__ == "__main__":

    myrank = MPI.COMM_WORLD.Get_rank()
    qpoints = None
    if myrank == 0:
        d3calc = D3Calc(configString = configString)
        print d3calc.pw.setting.get('pwInput')
        d3calc.pw.input.parse()
        d3calc.pw.input.save()
        qpointGrid = [2,2,2]
        qpoints = kmesh.kMeshCart( qpointGrid, \
                            d3calc.pw.input.structure.lattice.reciprocalBase() )
        nprocs = MPI.COMM_WORLD.Get_size()
        if nprocs > qpoints.shape[0]:
            nprocs = qpoints.shape[0]
        pointsPerProc = int(qpoints.shape[0]/nprocs)
        nprocs = pointsPerProc*qpoints.shape[0]
        print qpoints.reshape(nprocs, pointsPerProc, qpoints.shape[1])
        #for i in range(nprocs):
        #    qpt = qpoints[i:(i+pointsPerProc-1), :]
        #    MPI.COMM_WORLD.Send([qpt, MPI.FLOAT], dest=i, tag=66)
    #if myrank != 0:
    #    MPI.COMM_WORLD.Recv([qpt, MPI.FLOAT], source=0, tag=66)

    #   if pointsPerProc == 0:
    #        raise NameError('too many processors')





__author__="Nikolay Markovskiy"
__date__ ="$Dec 11, 2009 2:20:08 PM$"
