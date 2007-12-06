# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = pyphon
PACKAGE = libphon

PROJ_SAR = $(BLD_LIBDIR)/$(PACKAGE).$(EXT_SAR)
PROJ_DLL = $(BLD_LIBDIR)/$(PACKAGE).$(EXT_SO)
PROJ_TMPDIR = $(BLD_TMPDIR)/$(PROJECT)/$(PACKAGE)
PROJ_CLEAN += $(PROJ_SAR) $(PROJ_DLL)

EXT_F77 = f90

include local.def
include Listfile90.mk
include std-f90.def


PROJ_SRCS = \
    $(D90SOURCE) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build the library

all: $(PROJ_SAR) export

show:
	@echo " "
	@echo "  Environment"
	@echo "  -----------"
	@echo "             source tree : $(BLD_ROOT)"
	@echo "        operating system : $(PLATFORM_ID)"
	@echo "            build target : $(TARGET_ID)"
	@echo "               C compiler: $(CC_ID) as $(CC)"
	@echo "             C++ compiler: $(CXX_ID) as $(CXX)"
	@echo "         FORTRAN compiler: $(F77_ID) as $(F77)"
	@echo "         FORTRAN compiler: $(COMPILER_F90_NAME)"
	@echo " "

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# export

export:: export-libraries

EXPORT_LIBS = $(PROJ_SAR)
EXPORT_BINS = $(PROJ_DLL)


# version
# $Id: Make.mm,v 1.2 2006/08/30 02:08:14 cummings Exp $

#
# End of file
