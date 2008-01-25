# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = AbInitio
PACKAGE = AbiPhon/pyphon

BUILD_DIRS = \
    libphon \
    phonmodule \
    python \

RECURSE_DIRS = $(BUILD_DIRS)

all: 
	$(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse

tidy::
	BLD_ACTION="tidy" $(MM) recurse

distclean::


# version
# $Id: Make.mm,v 1.1.1.1 2006/08/01 20:47:15 patrickh Exp $

# End of file
