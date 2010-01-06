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
from qecalc.multiphononcalc import MultiPhononCalc

import sys
from qecalc.d3calc import D3Calc
from qeutils import kmesh

import subprocess
import shutil
import os
import numpy

configString = """
# all the relevant input files must be preconfiguered for specific tasks
# before using this class

[Launcher]
# parallelization parameters
# if this section is empty - serial mode is used
#paraPrefix:   mpiexec -n 8
#paraPostfix: -npool 8

#serialPrefix: mpiexec -n 1
#serialPostfix:

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
#paraTorqueParams: -l nodes=2:ppn=12 -N myjob -j oe
#serialTorqueParams: -l nodes=1:ppn=1 -N myjob -j oe

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

def calcIndices(npoints, nproc):

        idx = []
        pointsPerProc = int(npoints/nproc)
        remainder = npoints % nproc

        uninproc = nproc
        if remainder != 0:
            uninproc = nproc - 1
            remainder = npoints % uninproc
            pointsPerProc = int(npoints/ uninproc)

        for i in range(uninproc):
            idx.append( range(i*pointsPerProc, (i*pointsPerProc+pointsPerProc) ) )
        if remainder != 0:
            idx.append( range( pointsPerProc*uninproc, (pointsPerProc*uninproc+remainder) ) )

        return idx


def main():

    myrank = MPI.COMM_WORLD.Get_rank()
    qpoints = None
    d3calc = D3Calc(configString = configString)

    #if myrank == 0:
    mphon = MultiPhononCalc(configString = configString)
    mphon.ph.syncSetting()
    mphon.ph.output.parse()
    grid, qpoints_indep, qpoints_full = mphon.ph.output.property('qpoints')
    qpoints = numpy.array(qpoints_indep)
    nprocs = MPI.COMM_WORLD.Get_size()
    assert nprocs <= qpoints.shape[0], 'Number of processors exceeds number of q-points !'
#        if nprocs > qpoints.shape[0]:
#            print 'Warning ! Number of processors exceeds number of q-points '
#            nprocs = qpoints.shape[0]
#            raise
    if myrank == 0:
        print qpoints
        for i, irange in enumerate(calcIndices(qpoints.shape[0], nprocs)):
            #print irange
            #print qpoints[irange, :]
            if i != 0:
                MPI.COMM_WORLD.send(qpoints[irange, :], dest = i, tag = 66)
            else:
                qpnts = qpoints[irange, :]



    if myrank != 0:
        qpnts = MPI.COMM_WORLD.recv(source=0, tag=66)



    #print 'irank = ', myrank, '\tqpoints = ', qpnts


    taskOutDir = d3calc.pw.setting.get('outdir') + '/d3calc' + str(myrank)
    #print taskOutDir

    # Create irank dir

    # clean if exists:
    shutil.rmtree(taskOutDir, ignore_errors = True)
    # create dir and all subdirs from taskOutDir
    os.makedirs(taskOutDir)

    for task in d3calc.getAllTasks():
        task.syncSetting()
        task.setting.set('outdir', 'tmp/')
        taskFileList = task.setting.getExistingFiles()
        for fileName in taskFileList:
            shutil.copyfile(fileName, taskOutDir + '/' + \
                                                    os.path.basename(fileName))
        #print task.name(), '\t', taskFileList
    d3calc.pw.setting.set('pseudo_dir', './')
    os.chdir(taskOutDir)

    print  d3calc.pw.cmdLine()
    d3calc.pw.launch()

    d3data = []

    fildrho = d3calc.d3.setting.get('fildrho').strip('_G')
    #print d3calc.d3.setting.get('fildrho')
    #print d3calc.d3.input.namelist('inputph').param('fildrho')
    #print fildrho
    # calculation at Gamma:
    fildrhoG = fildrho +  '_G'
    d3calc.ph.setting.set('fildrho', fildrhoG)
    d3calc.d3.setting.set('fildrho', fildrhoG)
    d3calc.d3.setting.set('fild0rho', fildrhoG)
    for task in d3calc.getAllTasks()[1:]:
        task.input.qpoints.set([0.0, 0.0, 0.0])
        task.input.save()
        task.launch()
    #d3calc.launch()

    # non Gamma points
    d3calc.ph.setting.set('fildrho', fildrho)
    d3calc.d3.setting.set('fildrho', fildrho)
    for qpoint in qpnts:
        for task in d3calc.getAllTasks()[1:]:
            print 'launching task', task.name(), 'at', qpoint
            task.input.qpoints.set(qpoint)
            task.input.save()
            task.launch()
        d3data.append(d3calc.d3.output.property('d3 tensor'))

    d3data = MPI.COMM_WORLD.gather(d3data, root = 0)

    if myrank == 0:
        #flatten list of data:
        d3results = []
        for idata in d3data:
            for qdata in idata:
                d3results.append(qdata)
        #d3results = d3data
        print d3results
        #print numpy.array(d3results).shape


if __name__ == "__main__":
    main()







__author__="Nikolay Markovskiy"
__date__ ="$Dec 11, 2009 2:20:08 PM$"
