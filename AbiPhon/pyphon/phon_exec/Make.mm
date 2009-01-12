# -*- Makefile -*-
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                              Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005 All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def
include std-f90.def

PROJECT = phon

PROJ_BIN = $(BLD_BINDIR)/$(PROJECT)
PHON = phon.f90

EXPORT_BINS = $(PROJ_BIN)
PROJ_BINDIR = $(EXPORT_BINDIR)

LIBRARIES = $(LCXX_FORTRAN) -lphon -llapack -lblas
# -lg2c
#LIBRARIES = $(EXTERNAL_LIBS) $(LCXX_FORTRAN) 

EXT_F77 = f90

F77_FLAGS += -I../libphon

PROJ_CLEAN += $(PROJ_BIN)

all: $(PROJ_BIN) export

export:: export-binaries

install: $(PROJ_BIN)
	$(CP_F) $(PROJ_BIN) $(TOOLS_DIR)/bin

$(PROJ_BIN): $(PHON)
	$(F77) -I../libphon $(F77_FLAGS) $(PROJ_LCXX_FLAGS) -o $@ $< $(LF77_FLAGS) $(LF77_LIBPATH) $(LIBRARIES) 


# version
# $Id: Make.mm,v 1.1.1.1 2005/03/08 16:13:30 aivazis Exp $

# End of file
