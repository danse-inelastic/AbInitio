# -*- Makefile -*-

# Patrick Hung
# Calif. Inst. of Tech.


PROJECT = pyphon
PACKAGE = _pyphon
MODULE = _pyphon


include std-pythonmodule.def
include local.def

PROJ_CXX_SRCLIB = -llapack -lblas -lgfortran -lphon
#-lg2c 
#-l_pyphon 
PROJ_LIBRARIES += $(DEV_LCXX_LIBRARIES) $(COMPILER_LCXX_FORTRAN)


PROJ_SRCS = \
	dummy.cc \

#	_pyphon.cc \

export:: 

# version
# $Id: Make.mm,v 1.2 2006/09/06 21:19:14 cummings Exp $

# End of file
