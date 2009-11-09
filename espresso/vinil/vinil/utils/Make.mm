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
PACKAGE = utils


#--------------------------------------------------------------------------
#

EXPORT_PYTHON_MODULES = \
        const.py \
        dos.py \
        espressotable.py \
        jobstable.py \
        orderedDict.py \
        pathbuilder.py \
        plot.py \
        pw.py \
        simchain.py \
        simulationsteps.py \
        stepper.py \
        table.py \
        utils.py \
        __init__.py \

BUILD_DIRS = \
        parser \

OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)

all: export-package-python-modules #export-python-modules
	BLD_ACTION="all" $(MM) recurse

tidy::
	BLD_ACTION="tidy" $(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse

distclean::
	BLD_ACTION="distclean" $(MM) recurse

#export:: export-package-python-modules
#	BLD_ACTION="export" $(MM) recurse


# version
# $Id$

# End of file