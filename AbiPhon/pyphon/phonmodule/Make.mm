# -*- Makefile -*-

# Patrick Hung
# Calif. Inst. of Tech.


PROJECT = AbInitio
PACKAGE = phonmodule
MODULE = _pyphon


include std-pythonmodule.def
include local.def

PROJ_CXX_SRCLIB = -l_pyphon -llapack -lblas -lg2c -lgfortran -lphon
PROJ_LIBRARIES += $(DEV_LCXX_LIBRARIES) $(COMPILER_LCXX_FORTRAN)


PROJ_SRCS = \
	dummy.cc
	_pyphon.cc

export:: 

# version
# $Id: Make.mm,v 1.2 2006/09/06 21:19:14 cummings Exp $

# End of file
