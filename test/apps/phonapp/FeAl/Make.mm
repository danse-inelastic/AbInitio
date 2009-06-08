# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = AbInitio


PROJ_TIDY += \
	CHG \
	CHGCAR* \
	CONTCAR \
	DISP* \
	DOSCAR \
	DOS* \
	EIGENVAL \
	FORCES \
	IBZKPT \
	INCAR \
	INPHON \
	KPOINTS \
	OSZICAR \
	OUTCAR \
	PCDAT \
	POSCAR* \
	QPOINTS* \
	SPOSCAR \
	POTCAR \
	XDATCAR \
	WAVECAR \
	*.out \
	out.* \
	out_* \
	runhf* \
	vapsrun.xml* \


#--------------------------------------------------------------------------
#

BUILD_DIRS = \

OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)


#--------------------------------------------------------------------------
#


all: 
	BLD_ACTION="all" $(MM) recurse

tidy::
	BLD_ACTION="tidy" $(MM) recurse


# version
# $Id$

# End of file
