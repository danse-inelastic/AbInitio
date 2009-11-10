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

# Temporary solution to create jobs tables for all jobs actors
# TODO: Refactor so that the table creation was more general (useful for actors other than jobs)

from luban.content import load
from luban.content.Link import Link

# [<Job ID> | <Simulation> | Server | User | Submitted | Status | <Delete> | <Check>]
# 8 columns

def tableJobs(headers, jobids, simids, simnames, columns):
    from luban.content.table import Table, Model, View

    # create a model class
    class model(Model):
        jobid       = Model.descriptors.link(name='jobid')
        simulation  = Model.descriptors.link(name='simulation')
        server      = Model.descriptors.str(name='server')
        user        = Model.descriptors.str(name='user')
        submitted   = Model.descriptors.str(name='submitted')
        status      = Model.descriptors.str(name='status')
        delete      = Model.descriptors.link(name='delete')
        check       = Model.descriptors.link(name='check')

    # create a view
    view = View( columns =  [ View.Column(label=headers[0], measure='jobid'),
                              View.Column(label=headers[1], measure='simulation'),
                              View.Column(label=headers[2], measure='server'),
                              View.Column(label=headers[3], measure='user'),
                              View.Column(label=headers[4], measure='submitted'),
                              View.Column(label=headers[5], measure='status'),
                              View.Column(label=headers[6], measure='delete'),
                              View.Column(label=headers[7], measure='check'),]
                              )
    # Populate the data list
    def jobid(i):
        link = Link(label=jobids[i], onclick = load(actor='jobs/view', routine='link', id=jobids[i]))
        return link

    def sim(i):
        link = Link(label=simnames[i], onclick = load(actor='espresso/sim-view', routine='link', id=simids[i])) 
        return link

    def delete(i):
        link = Link(label="Delete", onclick = load(actor='jobs/delete', routine='link', id=jobids[i]))
        return link

    def check(i):
        link = Link(label="Check", onclick = load(actor='jobs/index', routine='link', id=jobids[i]))       #ids[i]
        return link

    data    = []
    for i in range(len(jobids)):    # just need size of list (e.g. from jobids)
        n           = [jobid(i)]
        data.append(n)
        data[i]     += [sim(i)]
        data[i]     += columns[i]
        data[i]     += [check(i)]
        data[i]     += [delete(i)]

    # create the table
    table = Table(model=model, data=data, view=view)

    return table

__date__ = "$Nov 4, 2009 9:52:45 PM$"


