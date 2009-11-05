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

from luban.content.Document import Document
from luban.content import load
from luban.content.Link import Link
from luban.content.Splitter import Splitter
from luban.content.FormSelectorField import FormSelectorField
from luban.content.FormTextField import FormTextField


def tableController(headers, actor):
    s               = Splitter(Class="table-controller", orientation='horizontal')
    s_filter        = s.section()
    s_filter.add(FormTextField(label='Filter', Class="table-filter"))
    s_formselect    = s.section()
    s_formselect.add(FormSelectorField(entries=enumerate(headers)) )

    s_pagination    = s.section()
    p               = Document(Class="table-pagination")
    p.add(Link(label="Prev", Class="pagination-link", onclick = load(actor=actor, routine='link')) )
    p.add(Link(label="Next", Class="pagination-link", onclick = load(actor=actor, routine='link')) )

    s_pagination.add(p)

    return s


__date__ = "$Nov 4, 2009 10:05:44 PM$"


