#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                  Jiao Lin
#                     California Institute of Technology
#                       (C) 2008  All Rights Reserved
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

external_depositories = []

# I nowhere use 'depositories' function. Do I actually need it?
def depositories(contentroot):
    return standard_depositories(contentroot) + external_depositories

def standard_depositories(root, excludes = ['data', '.svn']):
    entries = [os.path.join(root, e) for e in os.listdir(root)]
    directories = filter(
        lambda entry: entry not in excludes and os.path.isdir(entry),
        entries)
    return directories

import os

__date__ = "$Jul 21, 2009 11:45:20 AM$"


