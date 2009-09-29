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

PROJECT = vinil
PACKAGE = html

EXPORT_DATADIRS = \
	css \
	images \
	javascripts \

EXPORT_DATAFILES = \
        main-superapp.pml \ 

OTHERS = \

#--------------------------------------------------------------------------
#

all: export-package-data
	BLD_ACTION="all" $(MM) recurse

tidy::
	BLD_ACTION="tidy" $(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse

distclean::
	BLD_ACTION="distclean" $(MM) recurse



RSYNC_A = rsync -a
EXPORT_DATA_PATH = $(EXPORT_ROOT)/$(PACKAGE)


export-package-data: export-package-data-dirs export-package-data-files


export-package-data-dirs:: $(EXPORT_DATADIRS) 
	mkdir -p $(EXPORT_DATA_PATH); \
	for x in $(EXPORT_DATADIRS); do { \
            if [ -d $$x ]; then { \
	        $(RSYNC_A) $$x/ $(EXPORT_DATA_PATH)/$$x/ ; \
            } fi; \
        } done

export-package-data-files:: $(EXPORT_DATAFILES)
	mkdir -p $(EXPORT_DATA_PATH); \
	for x in $(EXPORT_DATAFILES); do { \
            if [ -f $$x ]; then { \
	        $(RSYNC_A) $$x $(EXPORT_DATA_PATH)/$$x ; \
            } fi; \
        } done

# version
# $Id: Make.mm,v 1.1.1.1 2006-11-27 00:09:14 aivazis Exp $

# End of file
