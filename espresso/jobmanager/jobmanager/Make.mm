# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2004  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = jobmanager
PACKAGE = jobmanager


# directory structure

BUILD_DIRS = \
    applications \
    components \
    simulations \
    utils \

OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)


#--------------------------------------------------------------------------
#

EXPORT_PYTHON_MODULES = \
	__init__.py \


all: export-python-modules
	BLD_ACTION="all" $(MM) recurse



# version
# $Id: Make.mm 1213 2006-11-18 16:17:08Z linjiao $

# End of file
