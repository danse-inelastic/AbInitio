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

# Temporary solution to create espresso tables 
# TODO: Refactor so that the table creation was more general

from luban.content import load
from luban.content.Link import Link

# 8 columns

def tableSimulations(headers, names, columns, ids):
    from luban.content.table import Table, Model, View

    # create a model class
    class model(Model):
        name        = Model.descriptors.link(name='name')
        type        = Model.descriptors.str(name='type')
        description = Model.descriptors.str(name='description')
        configs     = Model.descriptors.str(name='configs')
        created     = Model.descriptors.str(name='created')
        edit        = Model.descriptors.link(name='edit')
        remove      = Model.descriptors.link(name='remove')
        use         = Model.descriptors.link(name='use')

    # create a view
    view = View( columns =  [ View.Column(label=headers[0], measure='name'),
                              View.Column(label=headers[1], measure='type'),
                              View.Column(label=headers[2], measure='description'),
                              View.Column(label=headers[3], measure='configs'),
                              View.Column(label=headers[4], measure='created'),
                              View.Column(label=headers[5], measure='edit'),
                              View.Column(label=headers[6], measure='remove'),
                              View.Column(label=headers[7], measure='use'),]
                              )
    # Populate the data list
    def name(i):
        link = Link(label=names[i], onclick = load(actor='espresso/sim-view', routine='link', id=ids[i]))
        return link

    def edit(i):
        link = Link(label="Edit", onclick = load(actor='espresso/sim-edit', routine='link', id=ids[i]))
        return link

    def delete(i):
        link = Link(label="Delete", onclick = load(actor='espresso/sim-delete', routine='link', id=ids[i]))
        return link

    def use(i):
        link = Link(label="Clone", onclick = load(actor='espresso/sim-use', routine='link', id=ids[i]))
        return link

    data    = []
    for i in range(len(names)):
        n           = [name(i)]
        data.append(n)
        data[i]     += columns[i]
        data[i]     += [edit(i)]
        data[i]     += [delete(i)]
        data[i]     += [use(i)]

    # create the table
    table = Table(model=model, data=data, view=view)

    return table

__date__ = "$Nov 4, 2009 9:52:45 PM$"


