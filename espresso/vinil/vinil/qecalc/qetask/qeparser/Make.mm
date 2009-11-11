# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = vinil
PACKAGE = qecalc/qetask/qeparser


BUILD_DIRS = \
    inputs \
    outputs \

OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)


#--------------------------------------------------------------------------
#

all: export

tidy::
	BLD_ACTION="tidy" $(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse

distclean::
	BLD_ACTION="distclean" $(MM) recurse


#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    card.py \
    matdyninput.py \
    matdynqpoints.py \
    namelist.py \
    orderedDict.py \
    pwinput.py \
    pwkpoints.py \
    qeinput.py \
    qelattice.py \
    qeoutput.py \
    qeparser.py \
    qestructure.py \
    __init__.py \



export:: export-package-python-modules
	BLD_ACTION="export" $(MM) recurse


# version
# $Id: Make.mm,v 1.1.1.1 2006-11-27 00:09:19 aivazis Exp $

# End of file
