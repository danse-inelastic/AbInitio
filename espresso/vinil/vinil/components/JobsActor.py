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

from vinil.utils.jobstable import tableJobs
from vinil.utils.table import tableController
from vinil.utils.utils import stamp2date

import os
from luban.content import select
from luban.content.Paragraph import Paragraph
from luban.content.Document import Document
from luban.content import load
from luban.content.Link import Link
from luban.content.Splitter import Splitter
from luban.content.FormSelectorField import FormSelectorField
from luban.content.Plot2D import Plot2D

# Variable that should be set in envs.sh
EXPORT_ROOT = os.environ.get('EXPORT_ROOT')


from vinil.components.Actor import Actor as base
class JobsActor(base):

    class Inventory(base.Inventory):
        import pyre.inventory
        id      = pyre.inventory.str('id')


    # TODO: Refactor
    def getTableData(self, director, where=None):
        """Gets table data from Simulations database table to display in table"""
        jobs     = director.clerk.getJobs(where=where)
        # Take data from User, Simulations
        jobids      = []
        simids      = []
        simnames    = []
        columns     = []

        # Traversing jobs
        for i in range(len(jobs)):
            # Get simulation based on jobid
            (simid, simname)    = self.simdata(director, jobs[i])

            jobids.append(jobs[i].id)
            simids.append(simid)
            simnames.append(simname)
            row     = ["foxtrot.danse.us", "dexity", stamp2date(jobs[i].timeSubmitted), jobs[i].status ]
            columns.append(row)

        return (jobids, simids, simnames, columns)


    # TODO: Refactor
    def addTable(self, document, title, data):
        """
        Adds table to the document
        data = (jobids, simnames, columns)
        """
        #[<Job ID> | <Simulation> | Server | User | Submitted | Status | <Delete> | <Check>]

        headers = ("Job ID", "Simulation", "Server", "User", "Submitted", "Status", " ", " ")

        if len(data[0]) != 0:
            document.add(Paragraph(text=title, Class="header-h2"))
            document.add(tableController(headers, 'jobs/index'))    # temp
            document.add(tableJobs(headers, data[0], data[1], data[2], data[3]))


    def simdata(self, director, job):
        """Returns simulation data based on job"""
        simid   = None
        simname = ""

        if job:
            sim     = director.clerk.getSimulations(id=job.simulationId)
            if sim:
                simid   = sim.id
                simname = sim.sname

        return (simid, simname)

    def __init__(self, *args, **kwds):
        super(JobsActor, self).__init__(*args, **kwds)

        return


    def _configure(self):
        super(JobsActor, self)._configure()
        self.id = self.inventory.id
        return


    def _init(self):
        super(JobsActor, self)._init()
        return


__date__ = "$Nov 11, 2009 11:01:56 AM$"


