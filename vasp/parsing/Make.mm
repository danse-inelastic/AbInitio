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
PACKAGE = vasp/parsing

#--------------------------------------------------------------------------
#

BUILD_DIRS = \

OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)


#--------------------------------------------------------------------------
#


all: export
	BLD_ACTION="all" $(MM) recurse


#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
	__init__.py \
	Array.py \
	cmatrix.py \
	config.py \
	cp4vasp.py \
	cStructure.py \
	db.py \
	Dictionary.py \
	eigenvectorParser.py \
	FArray.py \
	graph.py \
	matrix.py \
	message.py \
	ODPdom.py \
	OrderedDict.py \
	parser.py \
	parser2.py \
	Property.py \
	repository.py \
	schedule.py \
	Selection.py \
	sellang.py \
	setupstore.py \
	setutils.py \
	SQLSystemPM.py \
	store.py \
	Structure.py \
	SystemPM.py \
	util.py \

BINDING_MODULES = \
	_cp4vasp.so \

<<<<<<< .mine
export-binding-modules::
	for i in $(BINDING_MODULES); do \
	cp $$i $(EXPORT_ROOT)/modules/$(PROJECT)/$(PACKAGE) ; \
	done


EXPORT_BINS = \
=======
EXPORT_BINS = _cp4vasp.so \
>>>>>>> .r30



export:: export-binaries release-binaries export-package-python-modules export-binding-modules


# version
# $Id$

# End of file
