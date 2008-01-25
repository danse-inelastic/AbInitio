# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

#import AbInitio.vasp.parsing._cp4vasp
import _cp4vasp

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class ClassInterface(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ClassInterface, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ClassInterface, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ClassInterface instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ClassInterface, 'this', _cp4vasp.new_ClassInterface(*args))
        _swig_setattr(self, ClassInterface, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.ClassInterface_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_ClassInterface):
        try:
            if self.thisown: destroy(self)
        except: pass


class ClassInterfacePtr(ClassInterface):
    def __init__(self, this):
        _swig_setattr(self, ClassInterface, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ClassInterface, 'thisown', 0)
        _swig_setattr(self, ClassInterface,self.__class__,ClassInterface)
_cp4vasp.ClassInterface_swigregister(ClassInterfacePtr)


getAtomtypesRecordHash = _cp4vasp.getAtomtypesRecordHash
class AtomtypesRecord(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, AtomtypesRecord, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, AtomtypesRecord, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ AtomtypesRecord instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getClassName(*args): return _cp4vasp.AtomtypesRecord_getClassName(*args)
    def __init__(self, *args):
        _swig_setattr(self, AtomtypesRecord, 'this', _cp4vasp.new_AtomtypesRecord(*args))
        _swig_setattr(self, AtomtypesRecord, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_AtomtypesRecord):
        try:
            if self.thisown: destroy(self)
        except: pass

    __swig_setmethods__["hash"] = _cp4vasp.AtomtypesRecord_hash_set
    __swig_getmethods__["hash"] = _cp4vasp.AtomtypesRecord_hash_get
    if _newclass:hash = property(_cp4vasp.AtomtypesRecord_hash_get, _cp4vasp.AtomtypesRecord_hash_set)
    __swig_setmethods__["atomspertype"] = _cp4vasp.AtomtypesRecord_atomspertype_set
    __swig_getmethods__["atomspertype"] = _cp4vasp.AtomtypesRecord_atomspertype_get
    if _newclass:atomspertype = property(_cp4vasp.AtomtypesRecord_atomspertype_get, _cp4vasp.AtomtypesRecord_atomspertype_set)
    __swig_setmethods__["mass"] = _cp4vasp.AtomtypesRecord_mass_set
    __swig_getmethods__["mass"] = _cp4vasp.AtomtypesRecord_mass_get
    if _newclass:mass = property(_cp4vasp.AtomtypesRecord_mass_get, _cp4vasp.AtomtypesRecord_mass_set)
    __swig_setmethods__["valence"] = _cp4vasp.AtomtypesRecord_valence_set
    __swig_getmethods__["valence"] = _cp4vasp.AtomtypesRecord_valence_get
    if _newclass:valence = property(_cp4vasp.AtomtypesRecord_valence_get, _cp4vasp.AtomtypesRecord_valence_set)
    __swig_setmethods__["radius"] = _cp4vasp.AtomtypesRecord_radius_set
    __swig_getmethods__["radius"] = _cp4vasp.AtomtypesRecord_radius_get
    if _newclass:radius = property(_cp4vasp.AtomtypesRecord_radius_get, _cp4vasp.AtomtypesRecord_radius_set)
    __swig_setmethods__["covalent"] = _cp4vasp.AtomtypesRecord_covalent_set
    __swig_getmethods__["covalent"] = _cp4vasp.AtomtypesRecord_covalent_get
    if _newclass:covalent = property(_cp4vasp.AtomtypesRecord_covalent_get, _cp4vasp.AtomtypesRecord_covalent_set)
    __swig_setmethods__["n"] = _cp4vasp.AtomtypesRecord_n_set
    __swig_getmethods__["n"] = _cp4vasp.AtomtypesRecord_n_get
    if _newclass:n = property(_cp4vasp.AtomtypesRecord_n_get, _cp4vasp.AtomtypesRecord_n_set)
    __swig_setmethods__["red"] = _cp4vasp.AtomtypesRecord_red_set
    __swig_getmethods__["red"] = _cp4vasp.AtomtypesRecord_red_get
    if _newclass:red = property(_cp4vasp.AtomtypesRecord_red_get, _cp4vasp.AtomtypesRecord_red_set)
    __swig_setmethods__["green"] = _cp4vasp.AtomtypesRecord_green_set
    __swig_getmethods__["green"] = _cp4vasp.AtomtypesRecord_green_get
    if _newclass:green = property(_cp4vasp.AtomtypesRecord_green_get, _cp4vasp.AtomtypesRecord_green_set)
    __swig_setmethods__["blue"] = _cp4vasp.AtomtypesRecord_blue_set
    __swig_getmethods__["blue"] = _cp4vasp.AtomtypesRecord_blue_get
    if _newclass:blue = property(_cp4vasp.AtomtypesRecord_blue_get, _cp4vasp.AtomtypesRecord_blue_set)
    __swig_setmethods__["hidden"] = _cp4vasp.AtomtypesRecord_hidden_set
    __swig_getmethods__["hidden"] = _cp4vasp.AtomtypesRecord_hidden_get
    if _newclass:hidden = property(_cp4vasp.AtomtypesRecord_hidden_get, _cp4vasp.AtomtypesRecord_hidden_set)
    __swig_setmethods__["selected"] = _cp4vasp.AtomtypesRecord_selected_set
    __swig_getmethods__["selected"] = _cp4vasp.AtomtypesRecord_selected_get
    if _newclass:selected = property(_cp4vasp.AtomtypesRecord_selected_get, _cp4vasp.AtomtypesRecord_selected_set)
    def getElement(*args): return _cp4vasp.AtomtypesRecord_getElement(*args)
    def setElement(*args): return _cp4vasp.AtomtypesRecord_setElement(*args)
    def getPPType(*args): return _cp4vasp.AtomtypesRecord_getPPType(*args)
    def getPPSpecie(*args): return _cp4vasp.AtomtypesRecord_getPPSpecie(*args)
    def getPPVersion(*args): return _cp4vasp.AtomtypesRecord_getPPVersion(*args)
    def getPseudopotential(*args): return _cp4vasp.AtomtypesRecord_getPseudopotential(*args)
    def setPPType(*args): return _cp4vasp.AtomtypesRecord_setPPType(*args)
    def setPPSpecie(*args): return _cp4vasp.AtomtypesRecord_setPPSpecie(*args)
    def setPPVersion(*args): return _cp4vasp.AtomtypesRecord_setPPVersion(*args)
    def setPseudopotential(*args): return _cp4vasp.AtomtypesRecord_setPseudopotential(*args)
    def clean(*args): return _cp4vasp.AtomtypesRecord_clean(*args)
    def setAtomtypesRecord(*args): return _cp4vasp.AtomtypesRecord_setAtomtypesRecord(*args)
    def clone(*args): return _cp4vasp.AtomtypesRecord_clone(*args)

class AtomtypesRecordPtr(AtomtypesRecord):
    def __init__(self, this):
        _swig_setattr(self, AtomtypesRecord, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, AtomtypesRecord, 'thisown', 0)
        _swig_setattr(self, AtomtypesRecord,self.__class__,AtomtypesRecord)
_cp4vasp.AtomtypesRecord_swigregister(AtomtypesRecordPtr)

class AtomInfo(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, AtomInfo, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, AtomInfo, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ AtomInfo instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getClassName(*args): return _cp4vasp.AtomInfo_getClassName(*args)
    __swig_getmethods__["types"] = _cp4vasp.AtomInfo_types_get
    if _newclass:types = property(_cp4vasp.AtomInfo_types_get)
    __swig_getmethods__["allocation_step"] = _cp4vasp.AtomInfo_allocation_step_get
    if _newclass:allocation_step = property(_cp4vasp.AtomInfo_allocation_step_get)
    def __init__(self, *args):
        _swig_setattr(self, AtomInfo, 'this', _cp4vasp.new_AtomInfo(*args))
        _swig_setattr(self, AtomInfo, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_AtomInfo):
        try:
            if self.thisown: destroy(self)
        except: pass

    def append(*args): return _cp4vasp.AtomInfo_append(*args)
    def realloc(*args): return _cp4vasp.AtomInfo_realloc(*args)
    def allocate(*args): return _cp4vasp.AtomInfo_allocate(*args)
    def getRecord(*args): return _cp4vasp.AtomInfo_getRecord(*args)
    def setRecord(*args): return _cp4vasp.AtomInfo_setRecord(*args)
    def getRecordForAtom(*args): return _cp4vasp.AtomInfo_getRecordForAtom(*args)
    def getRecordForElement(*args): return _cp4vasp.AtomInfo_getRecordForElement(*args)
    def getRecordForElementSafe(*args): return _cp4vasp.AtomInfo_getRecordForElementSafe(*args)
    def speciesIndex(*args): return _cp4vasp.AtomInfo_speciesIndex(*args)
    def getNatoms(*args): return _cp4vasp.AtomInfo_getNatoms(*args)
    def setAtomInfo(*args): return _cp4vasp.AtomInfo_setAtomInfo(*args)
    def clean(*args): return _cp4vasp.AtomInfo_clean(*args)
    def clone(*args): return _cp4vasp.AtomInfo_clone(*args)
    def delitem(*args): return _cp4vasp.AtomInfo_delitem(*args)
    def len(*args): return _cp4vasp.AtomInfo_len(*args)
    def fillAttributesWithTable(*args): return _cp4vasp.AtomInfo_fillAttributesWithTable(*args)

class AtomInfoPtr(AtomInfo):
    def __init__(self, this):
        _swig_setattr(self, AtomInfo, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, AtomInfo, 'thisown', 0)
        _swig_setattr(self, AtomInfo,self.__class__,AtomInfo)
_cp4vasp.AtomInfo_swigregister(AtomInfoPtr)


setvec3d = _cp4vasp.setvec3d

deletevec3d = _cp4vasp.deletevec3d

createmat3d = _cp4vasp.createmat3d

setmat3d = _cp4vasp.setmat3d

deletemat3d = _cp4vasp.deletemat3d

getVecElement3d = _cp4vasp.getVecElement3d

setVecElement3d = _cp4vasp.setVecElement3d

getMatVecElement3d = _cp4vasp.getMatVecElement3d

setMatVecElement3d = _cp4vasp.setMatVecElement3d

getMatElement3d = _cp4vasp.getMatElement3d

setMatElement3d = _cp4vasp.setMatElement3d

add3d = _cp4vasp.add3d

plus3d = _cp4vasp.plus3d

createplus3d = _cp4vasp.createplus3d

createplusmat3d = _cp4vasp.createplusmat3d

sub3d = _cp4vasp.sub3d

minus3d = _cp4vasp.minus3d

createminus3d = _cp4vasp.createminus3d

createminusmat3d = _cp4vasp.createminusmat3d

neg3d = _cp4vasp.neg3d

createneg3d = _cp4vasp.createneg3d

createnegmat3d = _cp4vasp.createnegmat3d

scalmul3d = _cp4vasp.scalmul3d

createscalmultiply3d = _cp4vasp.createscalmultiply3d

scaldiv3d = _cp4vasp.scaldiv3d

createscaldivide3d = _cp4vasp.createscaldivide3d

copy3d = _cp4vasp.copy3d

clone3d = _cp4vasp.clone3d

copymat3d = _cp4vasp.copymat3d

clonemat3d = _cp4vasp.clonemat3d

veclength3d = _cp4vasp.veclength3d

normalize3d = _cp4vasp.normalize3d

scalprod3d = _cp4vasp.scalprod3d

crossprod3d = _cp4vasp.crossprod3d

createcrossprod3d = _cp4vasp.createcrossprod3d

createmultiplymatscal3d = _cp4vasp.createmultiplymatscal3d

createmultiplymatvec3d = _cp4vasp.createmultiplymatvec3d

multiplymatvec3d = _cp4vasp.multiplymatvec3d

mulmatvec3d = _cp4vasp.mulmatvec3d

multiplymatmat3d = _cp4vasp.multiplymatmat3d

createmultiplymatmat3d = _cp4vasp.createmultiplymatmat3d

mulmatmat3d = _cp4vasp.mulmatmat3d

createrotmat3d = _cp4vasp.createrotmat3d

createrotmat3da = _cp4vasp.createrotmat3da

identitymat3d = _cp4vasp.identitymat3d

createidentitymat3d = _cp4vasp.createidentitymat3d

zeromat3d = _cp4vasp.zeromat3d

createzeromat3d = _cp4vasp.createzeromat3d

detmat3d = _cp4vasp.detmat3d

transmat3d = _cp4vasp.transmat3d
class FArray1D(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FArray1D, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FArray1D, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ FArray1D instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, FArray1D, 'this', _cp4vasp.new_FArray1D(*args))
        _swig_setattr(self, FArray1D, 'thisown', 1)
    def clear(*args): return _cp4vasp.FArray1D_clear(*args)
    def size(*args): return _cp4vasp.FArray1D_size(*args)
    def get(*args): return _cp4vasp.FArray1D_get(*args)
    def set(*args): return _cp4vasp.FArray1D_set(*args)
    def printrepr(*args): return _cp4vasp.FArray1D_printrepr(*args)
    def parseString(*args): return _cp4vasp.FArray1D_parseString(*args)
    def cloneBuff(*args): return _cp4vasp.FArray1D_cloneBuff(*args)
    def clone(*args): return _cp4vasp.FArray1D_clone(*args)
    def getMinimum(*args): return _cp4vasp.FArray1D_getMinimum(*args)
    def getMaximum(*args): return _cp4vasp.FArray1D_getMaximum(*args)
    def getAverage(*args): return _cp4vasp.FArray1D_getAverage(*args)
    def getVariance(*args): return _cp4vasp.FArray1D_getVariance(*args)
    def getSigma(*args): return _cp4vasp.FArray1D_getSigma(*args)
    def getClassName(*args): return _cp4vasp.FArray1D_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_FArray1D):
        try:
            if self.thisown: destroy(self)
        except: pass


class FArray1DPtr(FArray1D):
    def __init__(self, this):
        _swig_setattr(self, FArray1D, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FArray1D, 'thisown', 0)
        _swig_setattr(self, FArray1D,self.__class__,FArray1D)
_cp4vasp.FArray1D_swigregister(FArray1DPtr)

class FArray2D(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FArray2D, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FArray2D, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ FArray2D instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, FArray2D, 'this', _cp4vasp.new_FArray2D(*args))
        _swig_setattr(self, FArray2D, 'thisown', 1)
    def getMinimum(*args): return _cp4vasp.FArray2D_getMinimum(*args)
    def getMaximum(*args): return _cp4vasp.FArray2D_getMaximum(*args)
    def getAverage(*args): return _cp4vasp.FArray2D_getAverage(*args)
    def getVariance(*args): return _cp4vasp.FArray2D_getVariance(*args)
    def getSigma(*args): return _cp4vasp.FArray2D_getSigma(*args)
    def sizeX(*args): return _cp4vasp.FArray2D_sizeX(*args)
    def sizeY(*args): return _cp4vasp.FArray2D_sizeY(*args)
    def clear(*args): return _cp4vasp.FArray2D_clear(*args)
    def get(*args): return _cp4vasp.FArray2D_get(*args)
    def set(*args): return _cp4vasp.FArray2D_set(*args)
    def printrepr(*args): return _cp4vasp.FArray2D_printrepr(*args)
    def parseString(*args): return _cp4vasp.FArray2D_parseString(*args)
    def getArray(*args): return _cp4vasp.FArray2D_getArray(*args)
    def cloneBuff(*args): return _cp4vasp.FArray2D_cloneBuff(*args)
    def clone(*args): return _cp4vasp.FArray2D_clone(*args)
    def cubicInterpolation(*args): return _cp4vasp.FArray2D_cubicInterpolation(*args)
    def smear(*args): return _cp4vasp.FArray2D_smear(*args)
    def cloneVector(*args): return _cp4vasp.FArray2D_cloneVector(*args)
    def getClassName(*args): return _cp4vasp.FArray2D_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_FArray2D):
        try:
            if self.thisown: destroy(self)
        except: pass


class FArray2DPtr(FArray2D):
    def __init__(self, this):
        _swig_setattr(self, FArray2D, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FArray2D, 'thisown', 0)
        _swig_setattr(self, FArray2D,self.__class__,FArray2D)
_cp4vasp.FArray2D_swigregister(FArray2DPtr)


createFArray1Dsimple = _cp4vasp.createFArray1Dsimple

createStructure = _cp4vasp.createStructure

createStructureN = _cp4vasp.createStructureN
class Structure(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Structure, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Structure, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ Structure instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getClassName(*args): return _cp4vasp.Structure_getClassName(*args)
    __swig_setmethods__["scaling_flag"] = _cp4vasp.Structure_scaling_flag_set
    __swig_getmethods__["scaling_flag"] = _cp4vasp.Structure_scaling_flag_get
    if _newclass:scaling_flag = property(_cp4vasp.Structure_scaling_flag_get, _cp4vasp.Structure_scaling_flag_set)
    __swig_setmethods__["allocation_step"] = _cp4vasp.Structure_allocation_step_set
    __swig_getmethods__["allocation_step"] = _cp4vasp.Structure_allocation_step_get
    if _newclass:allocation_step = property(_cp4vasp.Structure_allocation_step_get, _cp4vasp.Structure_allocation_step_set)
    __swig_getmethods__["scaling"] = _cp4vasp.Structure_scaling_get
    if _newclass:scaling = property(_cp4vasp.Structure_scaling_get)
    __swig_getmethods__["basis"] = _cp4vasp.Structure_basis_get
    if _newclass:basis = property(_cp4vasp.Structure_basis_get)
    __swig_getmethods__["rbasis"] = _cp4vasp.Structure_rbasis_get
    if _newclass:rbasis = property(_cp4vasp.Structure_rbasis_get)
    __swig_getmethods__["total_number_of_atoms"] = _cp4vasp.Structure_total_number_of_atoms_get
    if _newclass:total_number_of_atoms = property(_cp4vasp.Structure_total_number_of_atoms_get)
    __swig_getmethods__["allocated"] = _cp4vasp.Structure_allocated_get
    if _newclass:allocated = property(_cp4vasp.Structure_allocated_get)
    __swig_getmethods__["info"] = _cp4vasp.Structure_info_get
    if _newclass:info = property(_cp4vasp.Structure_info_get)
    __swig_setmethods__["comment"] = _cp4vasp.Structure_comment_set
    __swig_getmethods__["comment"] = _cp4vasp.Structure_comment_get
    if _newclass:comment = property(_cp4vasp.Structure_comment_get, _cp4vasp.Structure_comment_set)
    __swig_setmethods__["coordinates"] = _cp4vasp.Structure_coordinates_set
    __swig_getmethods__["coordinates"] = _cp4vasp.Structure_coordinates_get
    if _newclass:coordinates = property(_cp4vasp.Structure_coordinates_get, _cp4vasp.Structure_coordinates_set)
    def __init__(self, *args):
        _swig_setattr(self, Structure, 'this', _cp4vasp.new_Structure(*args))
        _swig_setattr(self, Structure, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_Structure):
        try:
            if self.thisown: destroy(self)
        except: pass

    def correctScaling(*args): return _cp4vasp.Structure_correctScaling(*args)
    def parse(*args): return _cp4vasp.Structure_parse(*args)
    def read(*args): return _cp4vasp.Structure_read(*args)
    def write(*args): return _cp4vasp.Structure_write(*args)
    def toString(*args): return _cp4vasp.Structure_toString(*args)
    def isSelective(*args): return _cp4vasp.Structure_isSelective(*args)
    def isCarthesian(*args): return _cp4vasp.Structure_isCarthesian(*args)
    def isDirect(*args): return _cp4vasp.Structure_isDirect(*args)
    def getSelectiveDOF(*args): return _cp4vasp.Structure_getSelectiveDOF(*args)
    def setSelectiveDOF(*args): return _cp4vasp.Structure_setSelectiveDOF(*args)
    def setSelective(*args): return _cp4vasp.Structure_setSelective(*args)
    def updateRecipBasis(*args): return _cp4vasp.Structure_updateRecipBasis(*args)
    def getRecipBasis(*args): return _cp4vasp.Structure_getRecipBasis(*args)
    def forceConvertToCarthesian(*args): return _cp4vasp.Structure_forceConvertToCarthesian(*args)
    def forceConvertToDirect(*args): return _cp4vasp.Structure_forceConvertToDirect(*args)
    def setCarthesian(*args): return _cp4vasp.Structure_setCarthesian(*args)
    def setDirect(*args): return _cp4vasp.Structure_setDirect(*args)
    def dir2cart(*args): return _cp4vasp.Structure_dir2cart(*args)
    def cart2dir(*args): return _cp4vasp.Structure_cart2dir(*args)
    def dirVectorToUnitCell(*args): return _cp4vasp.Structure_dirVectorToUnitCell(*args)
    def dirVectorToCenteredUnitCell(*args): return _cp4vasp.Structure_dirVectorToCenteredUnitCell(*args)
    def cartVectorToUnitCell(*args): return _cp4vasp.Structure_cartVectorToUnitCell(*args)
    def cartVectorToCenteredUnitCell(*args): return _cp4vasp.Structure_cartVectorToCenteredUnitCell(*args)
    def vectorToUnitCell(*args): return _cp4vasp.Structure_vectorToUnitCell(*args)
    def vectorToCenteredUnitCell(*args): return _cp4vasp.Structure_vectorToCenteredUnitCell(*args)
    def toUnitCell(*args): return _cp4vasp.Structure_toUnitCell(*args)
    def toCenteredUnitCell(*args): return _cp4vasp.Structure_toCenteredUnitCell(*args)
    def mindistCartVectors(*args): return _cp4vasp.Structure_mindistCartVectors(*args)
    def mindistDirVectors(*args): return _cp4vasp.Structure_mindistDirVectors(*args)
    def createMindistMatrix(*args): return _cp4vasp.Structure_createMindistMatrix(*args)
    def deleteMindistMatrix(*args): return _cp4vasp.Structure_deleteMindistMatrix(*args)
    def getMindist(*args): return _cp4vasp.Structure_getMindist(*args)
    def clean(*args): return _cp4vasp.Structure_clean(*args)
    def setStructure(*args): return _cp4vasp.Structure_setStructure(*args)
    def clone(*args): return _cp4vasp.Structure_clone(*args)
    def len(*args): return _cp4vasp.Structure_len(*args)
    def getNumberOfSpecies(*args): return _cp4vasp.Structure_getNumberOfSpecies(*args)
    def getRecord(*args): return _cp4vasp.Structure_getRecord(*args)
    def get(*args): return _cp4vasp.Structure_get(*args)
    def set(*args): return _cp4vasp.Structure_set(*args)
    def allocate(*args): return _cp4vasp.Structure_allocate(*args)
    def realloc(*args): return _cp4vasp.Structure_realloc(*args)
    def append(*args): return _cp4vasp.Structure_append(*args)
    def delitem(*args): return _cp4vasp.Structure_delitem(*args)
    def setScaling(*args): return _cp4vasp.Structure_setScaling(*args)

class StructurePtr(Structure):
    def __init__(self, this):
        _swig_setattr(self, Structure, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Structure, 'thisown', 0)
        _swig_setattr(self, Structure,self.__class__,Structure)
_cp4vasp.Structure_swigregister(StructurePtr)

createvec3d = _cp4vasp.createvec3d

createFArray2Dsimple = _cp4vasp.createFArray2Dsimple

createFArray2DsimpleN = _cp4vasp.createFArray2DsimpleN

class Process(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Process, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Process, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ Process instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def total(*args): return _cp4vasp.Process_total(*args)
    def step(*args): return _cp4vasp.Process_step(*args)
    def status(*args): return _cp4vasp.Process_status(*args)
    def error(*args): return _cp4vasp.Process_error(*args)
    def next(*args): return _cp4vasp.Process_next(*args)
    def __del__(self, destroy=_cp4vasp.delete_Process):
        try:
            if self.thisown: destroy(self)
        except: pass

    def __init__(self, *args):
        _swig_setattr(self, Process, 'this', _cp4vasp.new_Process(*args))
        _swig_setattr(self, Process, 'thisown', 1)

class ProcessPtr(Process):
    def __init__(self, this):
        _swig_setattr(self, Process, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Process, 'thisown', 0)
        _swig_setattr(self, Process,self.__class__,Process)
_cp4vasp.Process_swigregister(ProcessPtr)

class ReadChgcarProcess(Process):
    __swig_setmethods__ = {}
    for _s in [Process]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ReadChgcarProcess, name, value)
    __swig_getmethods__ = {}
    for _s in [Process]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ReadChgcarProcess, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ ReadChgcarProcess instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def next(*args): return _cp4vasp.ReadChgcarProcess_next(*args)
    def __del__(self, destroy=_cp4vasp.delete_ReadChgcarProcess):
        try:
            if self.thisown: destroy(self)
        except: pass


class ReadChgcarProcessPtr(ReadChgcarProcess):
    def __init__(self, this):
        _swig_setattr(self, ReadChgcarProcess, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ReadChgcarProcess, 'thisown', 0)
        _swig_setattr(self, ReadChgcarProcess,self.__class__,ReadChgcarProcess)
_cp4vasp.ReadChgcarProcess_swigregister(ReadChgcarProcessPtr)

class ChgcarPlaneProcess(Process):
    __swig_setmethods__ = {}
    for _s in [Process]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ChgcarPlaneProcess, name, value)
    __swig_getmethods__ = {}
    for _s in [Process]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ChgcarPlaneProcess, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ ChgcarPlaneProcess instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def next(*args): return _cp4vasp.ChgcarPlaneProcess_next(*args)
    def __del__(self, destroy=_cp4vasp.delete_ChgcarPlaneProcess):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getPlane(*args): return _cp4vasp.ChgcarPlaneProcess_getPlane(*args)

class ChgcarPlaneProcessPtr(ChgcarPlaneProcess):
    def __init__(self, this):
        _swig_setattr(self, ChgcarPlaneProcess, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ChgcarPlaneProcess, 'thisown', 0)
        _swig_setattr(self, ChgcarPlaneProcess,self.__class__,ChgcarPlaneProcess)
_cp4vasp.ChgcarPlaneProcess_swigregister(ChgcarPlaneProcessPtr)

class Chgcar(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Chgcar, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Chgcar, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ Chgcar instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["plane_minimum"] = _cp4vasp.Chgcar_plane_minimum_set
    __swig_getmethods__["plane_minimum"] = _cp4vasp.Chgcar_plane_minimum_get
    if _newclass:plane_minimum = property(_cp4vasp.Chgcar_plane_minimum_get, _cp4vasp.Chgcar_plane_minimum_set)
    __swig_setmethods__["plane_maximum"] = _cp4vasp.Chgcar_plane_maximum_set
    __swig_getmethods__["plane_maximum"] = _cp4vasp.Chgcar_plane_maximum_get
    if _newclass:plane_maximum = property(_cp4vasp.Chgcar_plane_maximum_get, _cp4vasp.Chgcar_plane_maximum_set)
    __swig_setmethods__["plane_average"] = _cp4vasp.Chgcar_plane_average_set
    __swig_getmethods__["plane_average"] = _cp4vasp.Chgcar_plane_average_get
    if _newclass:plane_average = property(_cp4vasp.Chgcar_plane_average_get, _cp4vasp.Chgcar_plane_average_set)
    __swig_setmethods__["plane_variance"] = _cp4vasp.Chgcar_plane_variance_set
    __swig_getmethods__["plane_variance"] = _cp4vasp.Chgcar_plane_variance_get
    if _newclass:plane_variance = property(_cp4vasp.Chgcar_plane_variance_get, _cp4vasp.Chgcar_plane_variance_set)
    def getClassName(*args): return _cp4vasp.Chgcar_getClassName(*args)
    __swig_setmethods__["structure"] = _cp4vasp.Chgcar_structure_set
    __swig_getmethods__["structure"] = _cp4vasp.Chgcar_structure_get
    if _newclass:structure = property(_cp4vasp.Chgcar_structure_get, _cp4vasp.Chgcar_structure_set)
    __swig_setmethods__["nx"] = _cp4vasp.Chgcar_nx_set
    __swig_getmethods__["nx"] = _cp4vasp.Chgcar_nx_get
    if _newclass:nx = property(_cp4vasp.Chgcar_nx_get, _cp4vasp.Chgcar_nx_set)
    __swig_setmethods__["ny"] = _cp4vasp.Chgcar_ny_set
    __swig_getmethods__["ny"] = _cp4vasp.Chgcar_ny_get
    if _newclass:ny = property(_cp4vasp.Chgcar_ny_get, _cp4vasp.Chgcar_ny_set)
    __swig_setmethods__["nz"] = _cp4vasp.Chgcar_nz_set
    __swig_getmethods__["nz"] = _cp4vasp.Chgcar_nz_get
    if _newclass:nz = property(_cp4vasp.Chgcar_nz_get, _cp4vasp.Chgcar_nz_set)
    __swig_setmethods__["data"] = _cp4vasp.Chgcar_data_set
    __swig_getmethods__["data"] = _cp4vasp.Chgcar_data_get
    if _newclass:data = property(_cp4vasp.Chgcar_data_get, _cp4vasp.Chgcar_data_set)
    def __init__(self, *args):
        _swig_setattr(self, Chgcar, 'this', _cp4vasp.new_Chgcar(*args))
        _swig_setattr(self, Chgcar, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_Chgcar):
        try:
            if self.thisown: destroy(self)
        except: pass

    def subtractChgcar(*args): return _cp4vasp.Chgcar_subtractChgcar(*args)
    def calculateStatistics(*args): return _cp4vasp.Chgcar_calculateStatistics(*args)
    def getMinimum(*args): return _cp4vasp.Chgcar_getMinimum(*args)
    def getMaximum(*args): return _cp4vasp.Chgcar_getMaximum(*args)
    def getAverage(*args): return _cp4vasp.Chgcar_getAverage(*args)
    def getVariance(*args): return _cp4vasp.Chgcar_getVariance(*args)
    def getSigma(*args): return _cp4vasp.Chgcar_getSigma(*args)
    def clean(*args): return _cp4vasp.Chgcar_clean(*args)
    def read(*args): return _cp4vasp.Chgcar_read(*args)
    def createReadProcess(*args): return _cp4vasp.Chgcar_createReadProcess(*args)
    def write(*args): return _cp4vasp.Chgcar_write(*args)
    def get(*args): return _cp4vasp.Chgcar_get(*args)
    def getRaw(*args): return _cp4vasp.Chgcar_getRaw(*args)
    def set(*args): return _cp4vasp.Chgcar_set(*args)
    def setRaw(*args): return _cp4vasp.Chgcar_setRaw(*args)
    def getDirGrad(*args): return _cp4vasp.Chgcar_getDirGrad(*args)
    def getGrad(*args): return _cp4vasp.Chgcar_getGrad(*args)
    def downSampleByFactors(*args): return _cp4vasp.Chgcar_downSampleByFactors(*args)
    def sumElectrons(*args): return _cp4vasp.Chgcar_sumElectrons(*args)
    def gaussianSmearingX(*args): return _cp4vasp.Chgcar_gaussianSmearingX(*args)
    def gaussianSmearingY(*args): return _cp4vasp.Chgcar_gaussianSmearingY(*args)
    def gaussianSmearingZ(*args): return _cp4vasp.Chgcar_gaussianSmearingZ(*args)
    def clone(*args): return _cp4vasp.Chgcar_clone(*args)
    def setChgcar(*args): return _cp4vasp.Chgcar_setChgcar(*args)
    def calculatePlaneStatisticsX(*args): return _cp4vasp.Chgcar_calculatePlaneStatisticsX(*args)
    def calculatePlaneStatisticsY(*args): return _cp4vasp.Chgcar_calculatePlaneStatisticsY(*args)
    def calculatePlaneStatisticsZ(*args): return _cp4vasp.Chgcar_calculatePlaneStatisticsZ(*args)
    def searchMinPlaneX(*args): return _cp4vasp.Chgcar_searchMinPlaneX(*args)
    def searchMinPlaneY(*args): return _cp4vasp.Chgcar_searchMinPlaneY(*args)
    def searchMinPlaneZ(*args): return _cp4vasp.Chgcar_searchMinPlaneZ(*args)
    def getPlaneX(*args): return _cp4vasp.Chgcar_getPlaneX(*args)
    def getPlaneY(*args): return _cp4vasp.Chgcar_getPlaneY(*args)
    def getPlaneZ(*args): return _cp4vasp.Chgcar_getPlaneZ(*args)
    def createSmoothPlaneProcessX(*args): return _cp4vasp.Chgcar_createSmoothPlaneProcessX(*args)
    def createSmoothPlaneProcessY(*args): return _cp4vasp.Chgcar_createSmoothPlaneProcessY(*args)
    def createSmoothPlaneProcessZ(*args): return _cp4vasp.Chgcar_createSmoothPlaneProcessZ(*args)
    def createCCPlaneX(*args): return _cp4vasp.Chgcar_createCCPlaneX(*args)
    def createCCPlaneY(*args): return _cp4vasp.Chgcar_createCCPlaneY(*args)
    def createCCPlaneZ(*args): return _cp4vasp.Chgcar_createCCPlaneZ(*args)
    def createCCPlaneCubicX(*args): return _cp4vasp.Chgcar_createCCPlaneCubicX(*args)
    def createCCPlaneCubicY(*args): return _cp4vasp.Chgcar_createCCPlaneCubicY(*args)
    def createCCPlaneCubicZ(*args): return _cp4vasp.Chgcar_createCCPlaneCubicZ(*args)

class ChgcarPtr(Chgcar):
    def __init__(self, this):
        _swig_setattr(self, Chgcar, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Chgcar, 'thisown', 0)
        _swig_setattr(self, Chgcar,self.__class__,Chgcar)
_cp4vasp.Chgcar_swigregister(ChgcarPtr)

class ChgcarSmear(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ChgcarSmear, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ChgcarSmear, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ChgcarSmear instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ChgcarSmear, 'this', _cp4vasp.new_ChgcarSmear(*args))
        _swig_setattr(self, ChgcarSmear, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.ChgcarSmear_getClassName(*args)
    def setChgcar(*args): return _cp4vasp.ChgcarSmear_setChgcar(*args)
    def get(*args): return _cp4vasp.ChgcarSmear_get(*args)
    def __del__(self, destroy=_cp4vasp.delete_ChgcarSmear):
        try:
            if self.thisown: destroy(self)
        except: pass


class ChgcarSmearPtr(ChgcarSmear):
    def __init__(self, this):
        _swig_setattr(self, ChgcarSmear, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ChgcarSmear, 'thisown', 0)
        _swig_setattr(self, ChgcarSmear,self.__class__,ChgcarSmear)
_cp4vasp.ChgcarSmear_swigregister(ChgcarSmearPtr)

class GaussianChgcarSmear(ChgcarSmear):
    __swig_setmethods__ = {}
    for _s in [ChgcarSmear]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, GaussianChgcarSmear, name, value)
    __swig_getmethods__ = {}
    for _s in [ChgcarSmear]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, GaussianChgcarSmear, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ GaussianChgcarSmear instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["lx"] = _cp4vasp.GaussianChgcarSmear_lx_set
    __swig_getmethods__["lx"] = _cp4vasp.GaussianChgcarSmear_lx_get
    if _newclass:lx = property(_cp4vasp.GaussianChgcarSmear_lx_get, _cp4vasp.GaussianChgcarSmear_lx_set)
    __swig_setmethods__["ly"] = _cp4vasp.GaussianChgcarSmear_ly_set
    __swig_getmethods__["ly"] = _cp4vasp.GaussianChgcarSmear_ly_get
    if _newclass:ly = property(_cp4vasp.GaussianChgcarSmear_ly_get, _cp4vasp.GaussianChgcarSmear_ly_set)
    __swig_setmethods__["lz"] = _cp4vasp.GaussianChgcarSmear_lz_set
    __swig_getmethods__["lz"] = _cp4vasp.GaussianChgcarSmear_lz_get
    if _newclass:lz = property(_cp4vasp.GaussianChgcarSmear_lz_get, _cp4vasp.GaussianChgcarSmear_lz_set)
    __swig_setmethods__["dir"] = _cp4vasp.GaussianChgcarSmear_dir_set
    __swig_getmethods__["dir"] = _cp4vasp.GaussianChgcarSmear_dir_get
    if _newclass:dir = property(_cp4vasp.GaussianChgcarSmear_dir_get, _cp4vasp.GaussianChgcarSmear_dir_set)
    __swig_setmethods__["horizontal_sigma"] = _cp4vasp.GaussianChgcarSmear_horizontal_sigma_set
    __swig_getmethods__["horizontal_sigma"] = _cp4vasp.GaussianChgcarSmear_horizontal_sigma_get
    if _newclass:horizontal_sigma = property(_cp4vasp.GaussianChgcarSmear_horizontal_sigma_get, _cp4vasp.GaussianChgcarSmear_horizontal_sigma_set)
    __swig_setmethods__["vertical_sigma"] = _cp4vasp.GaussianChgcarSmear_vertical_sigma_set
    __swig_getmethods__["vertical_sigma"] = _cp4vasp.GaussianChgcarSmear_vertical_sigma_get
    if _newclass:vertical_sigma = property(_cp4vasp.GaussianChgcarSmear_vertical_sigma_get, _cp4vasp.GaussianChgcarSmear_vertical_sigma_set)
    def __init__(self, *args):
        _swig_setattr(self, GaussianChgcarSmear, 'this', _cp4vasp.new_GaussianChgcarSmear(*args))
        _swig_setattr(self, GaussianChgcarSmear, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.GaussianChgcarSmear_getClassName(*args)
    def setChgcar(*args): return _cp4vasp.GaussianChgcarSmear_setChgcar(*args)
    def get(*args): return _cp4vasp.GaussianChgcarSmear_get(*args)
    def __del__(self, destroy=_cp4vasp.delete_GaussianChgcarSmear):
        try:
            if self.thisown: destroy(self)
        except: pass


class GaussianChgcarSmearPtr(GaussianChgcarSmear):
    def __init__(self, this):
        _swig_setattr(self, GaussianChgcarSmear, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, GaussianChgcarSmear, 'thisown', 0)
        _swig_setattr(self, GaussianChgcarSmear,self.__class__,GaussianChgcarSmear)
_cp4vasp.GaussianChgcarSmear_swigregister(GaussianChgcarSmearPtr)

class ChgcarSmearProcess(Process):
    __swig_setmethods__ = {}
    for _s in [Process]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ChgcarSmearProcess, name, value)
    __swig_getmethods__ = {}
    for _s in [Process]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ChgcarSmearProcess, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ChgcarSmearProcess instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ChgcarSmearProcess, 'this', _cp4vasp.new_ChgcarSmearProcess(*args))
        _swig_setattr(self, ChgcarSmearProcess, 'thisown', 1)
    def next(*args): return _cp4vasp.ChgcarSmearProcess_next(*args)
    def get(*args): return _cp4vasp.ChgcarSmearProcess_get(*args)
    def getClassName(*args): return _cp4vasp.ChgcarSmearProcess_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_ChgcarSmearProcess):
        try:
            if self.thisown: destroy(self)
        except: pass


class ChgcarSmearProcessPtr(ChgcarSmearProcess):
    def __init__(self, this):
        _swig_setattr(self, ChgcarSmearProcess, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ChgcarSmearProcess, 'thisown', 0)
        _swig_setattr(self, ChgcarSmearProcess,self.__class__,ChgcarSmearProcess)
_cp4vasp.ChgcarSmearProcess_swigregister(ChgcarSmearProcessPtr)

class ChgcarSmearPlaneProcess(Process):
    __swig_setmethods__ = {}
    for _s in [Process]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ChgcarSmearPlaneProcess, name, value)
    __swig_getmethods__ = {}
    for _s in [Process]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ChgcarSmearPlaneProcess, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ChgcarSmearPlaneProcess instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ChgcarSmearPlaneProcess, 'this', _cp4vasp.new_ChgcarSmearPlaneProcess(*args))
        _swig_setattr(self, ChgcarSmearPlaneProcess, 'thisown', 1)
    def planeName(*args): return _cp4vasp.ChgcarSmearPlaneProcess_planeName(*args)
    def next(*args): return _cp4vasp.ChgcarSmearPlaneProcess_next(*args)
    def __del__(self, destroy=_cp4vasp.delete_ChgcarSmearPlaneProcess):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getClassName(*args): return _cp4vasp.ChgcarSmearPlaneProcess_getClassName(*args)
    def getPlane(*args): return _cp4vasp.ChgcarSmearPlaneProcess_getPlane(*args)

class ChgcarSmearPlaneProcessPtr(ChgcarSmearPlaneProcess):
    def __init__(self, this):
        _swig_setattr(self, ChgcarSmearPlaneProcess, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ChgcarSmearPlaneProcess, 'thisown', 0)
        _swig_setattr(self, ChgcarSmearPlaneProcess,self.__class__,ChgcarSmearPlaneProcess)
_cp4vasp.ChgcarSmearPlaneProcess_swigregister(ChgcarSmearPlaneProcessPtr)

class STMSearchProcess(Process):
    __swig_setmethods__ = {}
    for _s in [Process]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, STMSearchProcess, name, value)
    __swig_getmethods__ = {}
    for _s in [Process]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, STMSearchProcess, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ STMSearchProcess instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["mode"] = _cp4vasp.STMSearchProcess_mode_set
    __swig_getmethods__["mode"] = _cp4vasp.STMSearchProcess_mode_get
    if _newclass:mode = property(_cp4vasp.STMSearchProcess_mode_get, _cp4vasp.STMSearchProcess_mode_set)
    __swig_setmethods__["pstep"] = _cp4vasp.STMSearchProcess_pstep_set
    __swig_getmethods__["pstep"] = _cp4vasp.STMSearchProcess_pstep_get
    if _newclass:pstep = property(_cp4vasp.STMSearchProcess_pstep_get, _cp4vasp.STMSearchProcess_pstep_set)
    __swig_setmethods__["delta"] = _cp4vasp.STMSearchProcess_delta_set
    __swig_getmethods__["delta"] = _cp4vasp.STMSearchProcess_delta_get
    if _newclass:delta = property(_cp4vasp.STMSearchProcess_delta_get, _cp4vasp.STMSearchProcess_delta_set)
    __swig_setmethods__["n0"] = _cp4vasp.STMSearchProcess_n0_set
    __swig_getmethods__["n0"] = _cp4vasp.STMSearchProcess_n0_get
    if _newclass:n0 = property(_cp4vasp.STMSearchProcess_n0_get, _cp4vasp.STMSearchProcess_n0_set)
    __swig_setmethods__["autoplane"] = _cp4vasp.STMSearchProcess_autoplane_set
    __swig_getmethods__["autoplane"] = _cp4vasp.STMSearchProcess_autoplane_get
    if _newclass:autoplane = property(_cp4vasp.STMSearchProcess_autoplane_get, _cp4vasp.STMSearchProcess_autoplane_set)
    __swig_setmethods__["value"] = _cp4vasp.STMSearchProcess_value_set
    __swig_getmethods__["value"] = _cp4vasp.STMSearchProcess_value_get
    if _newclass:value = property(_cp4vasp.STMSearchProcess_value_get, _cp4vasp.STMSearchProcess_value_set)
    def __init__(self, *args):
        _swig_setattr(self, STMSearchProcess, 'this', _cp4vasp.new_STMSearchProcess(*args))
        _swig_setattr(self, STMSearchProcess, 'thisown', 1)
    def update(*args): return _cp4vasp.STMSearchProcess_update(*args)
    def getDir(*args): return _cp4vasp.STMSearchProcess_getDir(*args)
    def setDir(*args): return _cp4vasp.STMSearchProcess_setDir(*args)
    def setChgcar(*args): return _cp4vasp.STMSearchProcess_setChgcar(*args)
    def setSmear(*args): return _cp4vasp.STMSearchProcess_setSmear(*args)
    def searchFast(*args): return _cp4vasp.STMSearchProcess_searchFast(*args)
    def searchSlow(*args): return _cp4vasp.STMSearchProcess_searchSlow(*args)
    def getHeightFast(*args): return _cp4vasp.STMSearchProcess_getHeightFast(*args)
    def getHeightSlow(*args): return _cp4vasp.STMSearchProcess_getHeightSlow(*args)
    def getHeightFastCubic(*args): return _cp4vasp.STMSearchProcess_getHeightFastCubic(*args)
    def getHeightSlowCubic(*args): return _cp4vasp.STMSearchProcess_getHeightSlowCubic(*args)
    def getClassName(*args): return _cp4vasp.STMSearchProcess_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_STMSearchProcess):
        try:
            if self.thisown: destroy(self)
        except: pass

    def next(*args): return _cp4vasp.STMSearchProcess_next(*args)
    def processAll(*args): return _cp4vasp.STMSearchProcess_processAll(*args)
    def getPlane(*args): return _cp4vasp.STMSearchProcess_getPlane(*args)

class STMSearchProcessPtr(STMSearchProcess):
    def __init__(self, this):
        _swig_setattr(self, STMSearchProcess, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, STMSearchProcess, 'thisown', 0)
        _swig_setattr(self, STMSearchProcess,self.__class__,STMSearchProcess)
_cp4vasp.STMSearchProcess_swigregister(STMSearchProcessPtr)


VisInit = _cp4vasp.VisInit

VisMainLoop = _cp4vasp.VisMainLoop

VisMainLoopInThread = _cp4vasp.VisMainLoopInThread

VisSync = _cp4vasp.VisSync

VisCheck = _cp4vasp.VisCheck

checkThreadsSupport = _cp4vasp.checkThreadsSupport
class VisBackEvent(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisBackEvent, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VisBackEvent, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisBackEvent instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_getmethods__["create"] = lambda x: _cp4vasp.VisBackEvent_create
    if _newclass:create = staticmethod(_cp4vasp.VisBackEvent_create)
    __swig_setmethods__["type"] = _cp4vasp.VisBackEvent_type_set
    __swig_getmethods__["type"] = _cp4vasp.VisBackEvent_type_get
    if _newclass:type = property(_cp4vasp.VisBackEvent_type_get, _cp4vasp.VisBackEvent_type_set)
    __swig_setmethods__["index"] = _cp4vasp.VisBackEvent_index_set
    __swig_getmethods__["index"] = _cp4vasp.VisBackEvent_index_get
    if _newclass:index = property(_cp4vasp.VisBackEvent_index_get, _cp4vasp.VisBackEvent_index_set)
    __swig_setmethods__["nx"] = _cp4vasp.VisBackEvent_nx_set
    __swig_getmethods__["nx"] = _cp4vasp.VisBackEvent_nx_get
    if _newclass:nx = property(_cp4vasp.VisBackEvent_nx_get, _cp4vasp.VisBackEvent_nx_set)
    __swig_setmethods__["ny"] = _cp4vasp.VisBackEvent_ny_set
    __swig_getmethods__["ny"] = _cp4vasp.VisBackEvent_ny_get
    if _newclass:ny = property(_cp4vasp.VisBackEvent_ny_get, _cp4vasp.VisBackEvent_ny_set)
    __swig_setmethods__["nz"] = _cp4vasp.VisBackEvent_nz_set
    __swig_getmethods__["nz"] = _cp4vasp.VisBackEvent_nz_get
    if _newclass:nz = property(_cp4vasp.VisBackEvent_nz_get, _cp4vasp.VisBackEvent_nz_set)
    def getStructureDrawer(*args): return _cp4vasp.VisBackEvent_getStructureDrawer(*args)
    def getWindow(*args): return _cp4vasp.VisBackEvent_getWindow(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisBackEvent):
        try:
            if self.thisown: destroy(self)
        except: pass


class VisBackEventPtr(VisBackEvent):
    def __init__(self, this):
        _swig_setattr(self, VisBackEvent, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisBackEvent, 'thisown', 0)
        _swig_setattr(self, VisBackEvent,self.__class__,VisBackEvent)
_cp4vasp.VisBackEvent_swigregister(VisBackEventPtr)
cvar = _cp4vasp.cvar
BE_NONE = cvar.BE_NONE
BE_SELECTED = cvar.BE_SELECTED
BE_DESELECTED = cvar.BE_DESELECTED
BE_WIN_ACTIVATE = cvar.BE_WIN_ACTIVATE
BE_WIN_DEACTIVATE = cvar.BE_WIN_DEACTIVATE
BE_WIN_SHOW = cvar.BE_WIN_SHOW
BE_WIN_HIDE = cvar.BE_WIN_HIDE
BE_WIN_CLOSE = cvar.BE_WIN_CLOSE

VisBackEvent_create = _cp4vasp.VisBackEvent_create

class VisBackEventQueue(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisBackEventQueue, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VisBackEventQueue, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisBackEventQueue instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_getmethods__["get"] = lambda x: _cp4vasp.VisBackEventQueue_get
    if _newclass:get = staticmethod(_cp4vasp.VisBackEventQueue_get)
    def __init__(self, *args):
        _swig_setattr(self, VisBackEventQueue, 'this', _cp4vasp.new_VisBackEventQueue(*args))
        _swig_setattr(self, VisBackEventQueue, 'thisown', 1)
    def current(*args): return _cp4vasp.VisBackEventQueue_current(*args)
    def last(*args): return _cp4vasp.VisBackEventQueue_last(*args)
    def pop(*args): return _cp4vasp.VisBackEventQueue_pop(*args)
    def append(*args): return _cp4vasp.VisBackEventQueue_append(*args)
    def prepend(*args): return _cp4vasp.VisBackEventQueue_prepend(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisBackEventQueue):
        try:
            if self.thisown: destroy(self)
        except: pass


class VisBackEventQueuePtr(VisBackEventQueue):
    def __init__(self, this):
        _swig_setattr(self, VisBackEventQueue, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisBackEventQueue, 'thisown', 0)
        _swig_setattr(self, VisBackEventQueue,self.__class__,VisBackEventQueue)
_cp4vasp.VisBackEventQueue_swigregister(VisBackEventQueuePtr)

VisBackEventQueue_get = _cp4vasp.VisBackEventQueue_get

class VisWindow(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisWindow, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisWindow, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisWindow instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["mouse_x"] = _cp4vasp.VisWindow_mouse_x_set
    __swig_getmethods__["mouse_x"] = _cp4vasp.VisWindow_mouse_x_get
    if _newclass:mouse_x = property(_cp4vasp.VisWindow_mouse_x_get, _cp4vasp.VisWindow_mouse_x_set)
    __swig_setmethods__["mouse_y"] = _cp4vasp.VisWindow_mouse_y_set
    __swig_getmethods__["mouse_y"] = _cp4vasp.VisWindow_mouse_y_get
    if _newclass:mouse_y = property(_cp4vasp.VisWindow_mouse_y_get, _cp4vasp.VisWindow_mouse_y_set)
    __swig_setmethods__["mouse_button1"] = _cp4vasp.VisWindow_mouse_button1_set
    __swig_getmethods__["mouse_button1"] = _cp4vasp.VisWindow_mouse_button1_get
    if _newclass:mouse_button1 = property(_cp4vasp.VisWindow_mouse_button1_get, _cp4vasp.VisWindow_mouse_button1_set)
    __swig_setmethods__["mouse_button2"] = _cp4vasp.VisWindow_mouse_button2_set
    __swig_getmethods__["mouse_button2"] = _cp4vasp.VisWindow_mouse_button2_get
    if _newclass:mouse_button2 = property(_cp4vasp.VisWindow_mouse_button2_get, _cp4vasp.VisWindow_mouse_button2_set)
    __swig_setmethods__["mouse_button3"] = _cp4vasp.VisWindow_mouse_button3_set
    __swig_getmethods__["mouse_button3"] = _cp4vasp.VisWindow_mouse_button3_get
    if _newclass:mouse_button3 = property(_cp4vasp.VisWindow_mouse_button3_get, _cp4vasp.VisWindow_mouse_button3_set)
    __swig_setmethods__["mouse_button"] = _cp4vasp.VisWindow_mouse_button_set
    __swig_getmethods__["mouse_button"] = _cp4vasp.VisWindow_mouse_button_get
    if _newclass:mouse_button = property(_cp4vasp.VisWindow_mouse_button_get, _cp4vasp.VisWindow_mouse_button_set)
    __swig_setmethods__["key"] = _cp4vasp.VisWindow_key_set
    __swig_getmethods__["key"] = _cp4vasp.VisWindow_key_get
    if _newclass:key = property(_cp4vasp.VisWindow_key_get, _cp4vasp.VisWindow_key_set)
    def getClassName(*args): return _cp4vasp.VisWindow_getClassName(*args)
    __swig_getmethods__["x"] = _cp4vasp.VisWindow_x_get
    if _newclass:x = property(_cp4vasp.VisWindow_x_get)
    __swig_getmethods__["y"] = _cp4vasp.VisWindow_y_get
    if _newclass:y = property(_cp4vasp.VisWindow_y_get)
    __swig_getmethods__["w"] = _cp4vasp.VisWindow_w_get
    if _newclass:w = property(_cp4vasp.VisWindow_w_get)
    __swig_getmethods__["h"] = _cp4vasp.VisWindow_h_get
    if _newclass:h = property(_cp4vasp.VisWindow_h_get)
    def __init__(self, *args):
        _swig_setattr(self, VisWindow, 'this', _cp4vasp.new_VisWindow(*args))
        _swig_setattr(self, VisWindow, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_VisWindow):
        try:
            if self.thisown: destroy(self)
        except: pass

    __swig_getmethods__["getFirstWindow"] = lambda x: _cp4vasp.VisWindow_getFirstWindow
    if _newclass:getFirstWindow = staticmethod(_cp4vasp.VisWindow_getFirstWindow)
    __swig_getmethods__["getLastWindow"] = lambda x: _cp4vasp.VisWindow_getLastWindow
    if _newclass:getLastWindow = staticmethod(_cp4vasp.VisWindow_getLastWindow)
    __swig_getmethods__["getWindow"] = lambda x: _cp4vasp.VisWindow_getWindow
    if _newclass:getWindow = staticmethod(_cp4vasp.VisWindow_getWindow)
    __swig_getmethods__["getWindowByOutput"] = lambda x: _cp4vasp.VisWindow_getWindowByOutput
    if _newclass:getWindowByOutput = staticmethod(_cp4vasp.VisWindow_getWindowByOutput)
    __swig_getmethods__["windowsCount"] = lambda x: _cp4vasp.VisWindow_windowsCount
    if _newclass:windowsCount = staticmethod(_cp4vasp.VisWindow_windowsCount)
    __swig_getmethods__["getWindowIndex"] = lambda x: _cp4vasp.VisWindow_getWindowIndex
    if _newclass:getWindowIndex = staticmethod(_cp4vasp.VisWindow_getWindowIndex)
    __swig_getmethods__["deleteAllWindows"] = lambda x: _cp4vasp.VisWindow_deleteAllWindows
    if _newclass:deleteAllWindows = staticmethod(_cp4vasp.VisWindow_deleteAllWindows)
    def getNextWindow(*args): return _cp4vasp.VisWindow_getNextWindow(*args)
    def getPreviousWindow(*args): return _cp4vasp.VisWindow_getPreviousWindow(*args)
    __swig_getmethods__["deleteWindow"] = lambda x: _cp4vasp.VisWindow_deleteWindow
    if _newclass:deleteWindow = staticmethod(_cp4vasp.VisWindow_deleteWindow)
    def getTitle(*args): return _cp4vasp.VisWindow_getTitle(*args)
    def setTitle(*args): return _cp4vasp.VisWindow_setTitle(*args)
    def position(*args): return _cp4vasp.VisWindow_position(*args)
    def size(*args): return _cp4vasp.VisWindow_size(*args)
    def resize(*args): return _cp4vasp.VisWindow_resize(*args)
    def show(*args): return _cp4vasp.VisWindow_show(*args)
    def hide(*args): return _cp4vasp.VisWindow_hide(*args)
    def redraw(*args): return _cp4vasp.VisWindow_redraw(*args)
    def setDrawer(*args): return _cp4vasp.VisWindow_setDrawer(*args)
    def getDrawer(*args): return _cp4vasp.VisWindow_getDrawer(*args)
    def saveScreenshot(*args): return _cp4vasp.VisWindow_saveScreenshot(*args)

class VisWindowPtr(VisWindow):
    def __init__(self, this):
        _swig_setattr(self, VisWindow, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisWindow, 'thisown', 0)
        _swig_setattr(self, VisWindow,self.__class__,VisWindow)
_cp4vasp.VisWindow_swigregister(VisWindowPtr)

VisWindow_getFirstWindow = _cp4vasp.VisWindow_getFirstWindow

VisWindow_getLastWindow = _cp4vasp.VisWindow_getLastWindow

VisWindow_getWindow = _cp4vasp.VisWindow_getWindow

VisWindow_getWindowByOutput = _cp4vasp.VisWindow_getWindowByOutput

VisWindow_windowsCount = _cp4vasp.VisWindow_windowsCount

VisWindow_getWindowIndex = _cp4vasp.VisWindow_getWindowIndex

VisWindow_deleteAllWindows = _cp4vasp.VisWindow_deleteAllWindows

VisWindow_deleteWindow = _cp4vasp.VisWindow_deleteWindow

class VisDrawer(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisDrawer, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisDrawer, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisDrawer instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, VisDrawer, 'this', _cp4vasp.new_VisDrawer(*args))
        _swig_setattr(self, VisDrawer, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.VisDrawer_getClassName(*args)
    def getWindow(*args): return _cp4vasp.VisDrawer_getWindow(*args)
    def getPrevious(*args): return _cp4vasp.VisDrawer_getPrevious(*args)
    def getNext(*args): return _cp4vasp.VisDrawer_getNext(*args)
    def getFirst(*args): return _cp4vasp.VisDrawer_getFirst(*args)
    def getLast(*args): return _cp4vasp.VisDrawer_getLast(*args)
    def countBefore(*args): return _cp4vasp.VisDrawer_countBefore(*args)
    def countAfter(*args): return _cp4vasp.VisDrawer_countAfter(*args)
    def count(*args): return _cp4vasp.VisDrawer_count(*args)
    def insertAfter(*args): return _cp4vasp.VisDrawer_insertAfter(*args)
    def insertBefore(*args): return _cp4vasp.VisDrawer_insertBefore(*args)
    def insertSequenceAfter(*args): return _cp4vasp.VisDrawer_insertSequenceAfter(*args)
    def insertSequenceBefore(*args): return _cp4vasp.VisDrawer_insertSequenceBefore(*args)
    def append(*args): return _cp4vasp.VisDrawer_append(*args)
    def appendSequence(*args): return _cp4vasp.VisDrawer_appendSequence(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisDrawer):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getMouseX(*args): return _cp4vasp.VisDrawer_getMouseX(*args)
    def getMouseY(*args): return _cp4vasp.VisDrawer_getMouseY(*args)
    def getMouseButton(*args): return _cp4vasp.VisDrawer_getMouseButton(*args)
    def getMouseButton1(*args): return _cp4vasp.VisDrawer_getMouseButton1(*args)
    def getMouseButton2(*args): return _cp4vasp.VisDrawer_getMouseButton2(*args)
    def getMouseButton3(*args): return _cp4vasp.VisDrawer_getMouseButton3(*args)
    def getKey(*args): return _cp4vasp.VisDrawer_getKey(*args)
    def getWidth(*args): return _cp4vasp.VisDrawer_getWidth(*args)
    def getHeight(*args): return _cp4vasp.VisDrawer_getHeight(*args)
    def redraw(*args): return _cp4vasp.VisDrawer_redraw(*args)

class VisDrawerPtr(VisDrawer):
    def __init__(self, this):
        _swig_setattr(self, VisDrawer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisDrawer, 'thisown', 0)
        _swig_setattr(self, VisDrawer,self.__class__,VisDrawer)
_cp4vasp.VisDrawer_swigregister(VisDrawerPtr)

class VisNavDrawer(VisDrawer):
    __swig_setmethods__ = {}
    for _s in [VisDrawer]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisNavDrawer, name, value)
    __swig_getmethods__ = {}
    for _s in [VisDrawer]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisNavDrawer, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisNavDrawer instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_getmethods__["bg_red"] = _cp4vasp.VisNavDrawer_bg_red_get
    if _newclass:bg_red = property(_cp4vasp.VisNavDrawer_bg_red_get)
    __swig_getmethods__["bg_green"] = _cp4vasp.VisNavDrawer_bg_green_get
    if _newclass:bg_green = property(_cp4vasp.VisNavDrawer_bg_green_get)
    __swig_getmethods__["bg_blue"] = _cp4vasp.VisNavDrawer_bg_blue_get
    if _newclass:bg_blue = property(_cp4vasp.VisNavDrawer_bg_blue_get)
    def __init__(self, *args):
        _swig_setattr(self, VisNavDrawer, 'this', _cp4vasp.new_VisNavDrawer(*args))
        _swig_setattr(self, VisNavDrawer, 'thisown', 1)
    def setBackground(*args): return _cp4vasp.VisNavDrawer_setBackground(*args)
    def setAntialiasing(*args): return _cp4vasp.VisNavDrawer_setAntialiasing(*args)
    def getAntialiasing(*args): return _cp4vasp.VisNavDrawer_getAntialiasing(*args)
    def getClassName(*args): return _cp4vasp.VisNavDrawer_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisNavDrawer):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setHome(*args): return _cp4vasp.VisNavDrawer_setHome(*args)
    def setFrontView(*args): return _cp4vasp.VisNavDrawer_setFrontView(*args)
    def setBackView(*args): return _cp4vasp.VisNavDrawer_setBackView(*args)
    def setLeftView(*args): return _cp4vasp.VisNavDrawer_setLeftView(*args)
    def setRightView(*args): return _cp4vasp.VisNavDrawer_setRightView(*args)
    def setTopView(*args): return _cp4vasp.VisNavDrawer_setTopView(*args)
    def setBottomView(*args): return _cp4vasp.VisNavDrawer_setBottomView(*args)
    def getPerspective(*args): return _cp4vasp.VisNavDrawer_getPerspective(*args)
    def setPerspective(*args): return _cp4vasp.VisNavDrawer_setPerspective(*args)
    def getZoom(*args): return _cp4vasp.VisNavDrawer_getZoom(*args)
    def setZoom(*args): return _cp4vasp.VisNavDrawer_setZoom(*args)
    def mulZoom(*args): return _cp4vasp.VisNavDrawer_mulZoom(*args)

class VisNavDrawerPtr(VisNavDrawer):
    def __init__(self, this):
        _swig_setattr(self, VisNavDrawer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisNavDrawer, 'thisown', 0)
        _swig_setattr(self, VisNavDrawer,self.__class__,VisNavDrawer)
_cp4vasp.VisNavDrawer_swigregister(VisNavDrawerPtr)


getDefaultPrimitivesResolution = _cp4vasp.getDefaultPrimitivesResolution

setDefaultPrimitivesResolution = _cp4vasp.setDefaultPrimitivesResolution
class VisPrimitiveDrawer(VisDrawer):
    __swig_setmethods__ = {}
    for _s in [VisDrawer]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisPrimitiveDrawer, name, value)
    __swig_getmethods__ = {}
    for _s in [VisDrawer]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisPrimitiveDrawer, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisPrimitiveDrawer instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["arrow_radius"] = _cp4vasp.VisPrimitiveDrawer_arrow_radius_set
    __swig_getmethods__["arrow_radius"] = _cp4vasp.VisPrimitiveDrawer_arrow_radius_get
    if _newclass:arrow_radius = property(_cp4vasp.VisPrimitiveDrawer_arrow_radius_get, _cp4vasp.VisPrimitiveDrawer_arrow_radius_set)
    __swig_setmethods__["arrowhead_radius"] = _cp4vasp.VisPrimitiveDrawer_arrowhead_radius_set
    __swig_getmethods__["arrowhead_radius"] = _cp4vasp.VisPrimitiveDrawer_arrowhead_radius_get
    if _newclass:arrowhead_radius = property(_cp4vasp.VisPrimitiveDrawer_arrowhead_radius_get, _cp4vasp.VisPrimitiveDrawer_arrowhead_radius_set)
    __swig_setmethods__["arrowhead_length"] = _cp4vasp.VisPrimitiveDrawer_arrowhead_length_set
    __swig_getmethods__["arrowhead_length"] = _cp4vasp.VisPrimitiveDrawer_arrowhead_length_get
    if _newclass:arrowhead_length = property(_cp4vasp.VisPrimitiveDrawer_arrowhead_length_get, _cp4vasp.VisPrimitiveDrawer_arrowhead_length_set)
    def __init__(self, *args):
        _swig_setattr(self, VisPrimitiveDrawer, 'this', _cp4vasp.new_VisPrimitiveDrawer(*args))
        _swig_setattr(self, VisPrimitiveDrawer, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.VisPrimitiveDrawer_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisPrimitiveDrawer):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setPrimitivesResolution(*args): return _cp4vasp.VisPrimitiveDrawer_setPrimitivesResolution(*args)
    def color(*args): return _cp4vasp.VisPrimitiveDrawer_color(*args)
    def sphere(*args): return _cp4vasp.VisPrimitiveDrawer_sphere(*args)
    def cylinder(*args): return _cp4vasp.VisPrimitiveDrawer_cylinder(*args)
    def cone(*args): return _cp4vasp.VisPrimitiveDrawer_cone(*args)
    def line(*args): return _cp4vasp.VisPrimitiveDrawer_line(*args)
    def arrow(*args): return _cp4vasp.VisPrimitiveDrawer_arrow(*args)

class VisPrimitiveDrawerPtr(VisPrimitiveDrawer):
    def __init__(self, this):
        _swig_setattr(self, VisPrimitiveDrawer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisPrimitiveDrawer, 'thisown', 0)
        _swig_setattr(self, VisPrimitiveDrawer,self.__class__,VisPrimitiveDrawer)
_cp4vasp.VisPrimitiveDrawer_swigregister(VisPrimitiveDrawerPtr)

class AtomId(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, AtomId, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, AtomId, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ AtomId instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, AtomId, 'this', _cp4vasp.new_AtomId(*args))
        _swig_setattr(self, AtomId, 'thisown', 1)
    __swig_setmethods__["atom"] = _cp4vasp.AtomId_atom_set
    __swig_getmethods__["atom"] = _cp4vasp.AtomId_atom_get
    if _newclass:atom = property(_cp4vasp.AtomId_atom_get, _cp4vasp.AtomId_atom_set)
    __swig_setmethods__["nx"] = _cp4vasp.AtomId_nx_set
    __swig_getmethods__["nx"] = _cp4vasp.AtomId_nx_get
    if _newclass:nx = property(_cp4vasp.AtomId_nx_get, _cp4vasp.AtomId_nx_set)
    __swig_setmethods__["ny"] = _cp4vasp.AtomId_ny_set
    __swig_getmethods__["ny"] = _cp4vasp.AtomId_ny_get
    if _newclass:ny = property(_cp4vasp.AtomId_ny_get, _cp4vasp.AtomId_ny_set)
    __swig_setmethods__["nz"] = _cp4vasp.AtomId_nz_set
    __swig_getmethods__["nz"] = _cp4vasp.AtomId_nz_get
    if _newclass:nz = property(_cp4vasp.AtomId_nz_get, _cp4vasp.AtomId_nz_set)
    def __del__(self, destroy=_cp4vasp.delete_AtomId):
        try:
            if self.thisown: destroy(self)
        except: pass


class AtomIdPtr(AtomId):
    def __init__(self, this):
        _swig_setattr(self, AtomId, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, AtomId, 'thisown', 0)
        _swig_setattr(self, AtomId,self.__class__,AtomId)
_cp4vasp.AtomId_swigregister(AtomIdPtr)

class VisStructureDrawer(VisPrimitiveDrawer):
    __swig_setmethods__ = {}
    for _s in [VisPrimitiveDrawer]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisStructureDrawer, name, value)
    __swig_getmethods__ = {}
    for _s in [VisPrimitiveDrawer]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisStructureDrawer, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisStructureDrawer instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["cell_red"] = _cp4vasp.VisStructureDrawer_cell_red_set
    __swig_getmethods__["cell_red"] = _cp4vasp.VisStructureDrawer_cell_red_get
    if _newclass:cell_red = property(_cp4vasp.VisStructureDrawer_cell_red_get, _cp4vasp.VisStructureDrawer_cell_red_set)
    __swig_setmethods__["cell_green"] = _cp4vasp.VisStructureDrawer_cell_green_set
    __swig_getmethods__["cell_green"] = _cp4vasp.VisStructureDrawer_cell_green_get
    if _newclass:cell_green = property(_cp4vasp.VisStructureDrawer_cell_green_get, _cp4vasp.VisStructureDrawer_cell_green_set)
    __swig_setmethods__["cell_blue"] = _cp4vasp.VisStructureDrawer_cell_blue_set
    __swig_getmethods__["cell_blue"] = _cp4vasp.VisStructureDrawer_cell_blue_get
    if _newclass:cell_blue = property(_cp4vasp.VisStructureDrawer_cell_blue_get, _cp4vasp.VisStructureDrawer_cell_blue_set)
    __swig_setmethods__["bond_red"] = _cp4vasp.VisStructureDrawer_bond_red_set
    __swig_getmethods__["bond_red"] = _cp4vasp.VisStructureDrawer_bond_red_get
    if _newclass:bond_red = property(_cp4vasp.VisStructureDrawer_bond_red_get, _cp4vasp.VisStructureDrawer_bond_red_set)
    __swig_setmethods__["bond_green"] = _cp4vasp.VisStructureDrawer_bond_green_set
    __swig_getmethods__["bond_green"] = _cp4vasp.VisStructureDrawer_bond_green_get
    if _newclass:bond_green = property(_cp4vasp.VisStructureDrawer_bond_green_get, _cp4vasp.VisStructureDrawer_bond_green_set)
    __swig_setmethods__["bond_blue"] = _cp4vasp.VisStructureDrawer_bond_blue_set
    __swig_getmethods__["bond_blue"] = _cp4vasp.VisStructureDrawer_bond_blue_get
    if _newclass:bond_blue = property(_cp4vasp.VisStructureDrawer_bond_blue_get, _cp4vasp.VisStructureDrawer_bond_blue_set)
    __swig_getmethods__["info"] = _cp4vasp.VisStructureDrawer_info_get
    if _newclass:info = property(_cp4vasp.VisStructureDrawer_info_get)
    __swig_getmethods__["cell_line_width"] = _cp4vasp.VisStructureDrawer_cell_line_width_get
    if _newclass:cell_line_width = property(_cp4vasp.VisStructureDrawer_cell_line_width_get)
    __swig_getmethods__["showcellflag"] = _cp4vasp.VisStructureDrawer_showcellflag_get
    if _newclass:showcellflag = property(_cp4vasp.VisStructureDrawer_showcellflag_get)
    def updateStructure(*args): return _cp4vasp.VisStructureDrawer_updateStructure(*args)
    def setMultiple(*args): return _cp4vasp.VisStructureDrawer_setMultiple(*args)
    def setMultiple1(*args): return _cp4vasp.VisStructureDrawer_setMultiple1(*args)
    def setMultiple2(*args): return _cp4vasp.VisStructureDrawer_setMultiple2(*args)
    def setMultiple3(*args): return _cp4vasp.VisStructureDrawer_setMultiple3(*args)
    def getMultiple1(*args): return _cp4vasp.VisStructureDrawer_getMultiple1(*args)
    def getMultiple2(*args): return _cp4vasp.VisStructureDrawer_getMultiple2(*args)
    def getMultiple3(*args): return _cp4vasp.VisStructureDrawer_getMultiple3(*args)
    def __init__(self, *args):
        _swig_setattr(self, VisStructureDrawer, 'this', _cp4vasp.new_VisStructureDrawer(*args))
        _swig_setattr(self, VisStructureDrawer, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.VisStructureDrawer_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisStructureDrawer):
        try:
            if self.thisown: destroy(self)
        except: pass

    def switchSelectionByPick(*args): return _cp4vasp.VisStructureDrawer_switchSelectionByPick(*args)
    def setPrimitivesResolution(*args): return _cp4vasp.VisStructureDrawer_setPrimitivesResolution(*args)
    def selectObject(*args): return _cp4vasp.VisStructureDrawer_selectObject(*args)
    def fillInfo(*args): return _cp4vasp.VisStructureDrawer_fillInfo(*args)
    def setStructure(*args): return _cp4vasp.VisStructureDrawer_setStructure(*args)
    def getStructure(*args): return _cp4vasp.VisStructureDrawer_getStructure(*args)
    def getRadiusFactor(*args): return _cp4vasp.VisStructureDrawer_getRadiusFactor(*args)
    def setRadiusFactor(*args): return _cp4vasp.VisStructureDrawer_setRadiusFactor(*args)
    def getBondRadius(*args): return _cp4vasp.VisStructureDrawer_getBondRadius(*args)
    def setBondRadius(*args): return _cp4vasp.VisStructureDrawer_setBondRadius(*args)
    def getBondFactor(*args): return _cp4vasp.VisStructureDrawer_getBondFactor(*args)
    def setBondFactor(*args): return _cp4vasp.VisStructureDrawer_setBondFactor(*args)
    def showCell(*args): return _cp4vasp.VisStructureDrawer_showCell(*args)
    def getCellLineWidth(*args): return _cp4vasp.VisStructureDrawer_getCellLineWidth(*args)
    def setCellLineWidth(*args): return _cp4vasp.VisStructureDrawer_setCellLineWidth(*args)
    def setCellColor(*args): return _cp4vasp.VisStructureDrawer_setCellColor(*args)
    def setBondColor(*args): return _cp4vasp.VisStructureDrawer_setBondColor(*args)
    def getSelected(*args): return _cp4vasp.VisStructureDrawer_getSelected(*args)
    def getSelectedCount(*args): return _cp4vasp.VisStructureDrawer_getSelectedCount(*args)
    def appendSelected(*args): return _cp4vasp.VisStructureDrawer_appendSelected(*args)
    def notifySelected(*args): return _cp4vasp.VisStructureDrawer_notifySelected(*args)
    def notifyDeselected(*args): return _cp4vasp.VisStructureDrawer_notifyDeselected(*args)
    def removeSelectedAll(*args): return _cp4vasp.VisStructureDrawer_removeSelectedAll(*args)
    def removeSelectedItem(*args): return _cp4vasp.VisStructureDrawer_removeSelectedItem(*args)
    def findSelectedAtom(*args): return _cp4vasp.VisStructureDrawer_findSelectedAtom(*args)
    def selectAtom(*args): return _cp4vasp.VisStructureDrawer_selectAtom(*args)
    def deselectAtom(*args): return _cp4vasp.VisStructureDrawer_deselectAtom(*args)
    def switchAtomSelection(*args): return _cp4vasp.VisStructureDrawer_switchAtomSelection(*args)

class VisStructureDrawerPtr(VisStructureDrawer):
    def __init__(self, this):
        _swig_setattr(self, VisStructureDrawer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisStructureDrawer, 'thisown', 0)
        _swig_setattr(self, VisStructureDrawer,self.__class__,VisStructureDrawer)
_cp4vasp.VisStructureDrawer_swigregister(VisStructureDrawerPtr)

class Clamp(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Clamp, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Clamp, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ Clamp instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getClassName(*args): return _cp4vasp.Clamp_getClassName(*args)
    def f(*args): return _cp4vasp.Clamp_f(*args)
    def __del__(self, destroy=_cp4vasp.delete_Clamp):
        try:
            if self.thisown: destroy(self)
        except: pass


class ClampPtr(Clamp):
    def __init__(self, this):
        _swig_setattr(self, Clamp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Clamp, 'thisown', 0)
        _swig_setattr(self, Clamp,self.__class__,Clamp)
_cp4vasp.Clamp_swigregister(ClampPtr)

class ThresholdClamp(Clamp):
    __swig_setmethods__ = {}
    for _s in [Clamp]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ThresholdClamp, name, value)
    __swig_getmethods__ = {}
    for _s in [Clamp]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ThresholdClamp, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ThresholdClamp instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ThresholdClamp, 'this', _cp4vasp.new_ThresholdClamp(*args))
        _swig_setattr(self, ThresholdClamp, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ThresholdClamp):
        try:
            if self.thisown: destroy(self)
        except: pass


class ThresholdClampPtr(ThresholdClamp):
    def __init__(self, this):
        _swig_setattr(self, ThresholdClamp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ThresholdClamp, 'thisown', 0)
        _swig_setattr(self, ThresholdClamp,self.__class__,ThresholdClamp)
_cp4vasp.ThresholdClamp_swigregister(ThresholdClampPtr)

class SawtoothClamp(Clamp):
    __swig_setmethods__ = {}
    for _s in [Clamp]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, SawtoothClamp, name, value)
    __swig_getmethods__ = {}
    for _s in [Clamp]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, SawtoothClamp, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ SawtoothClamp instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, SawtoothClamp, 'this', _cp4vasp.new_SawtoothClamp(*args))
        _swig_setattr(self, SawtoothClamp, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_SawtoothClamp):
        try:
            if self.thisown: destroy(self)
        except: pass


class SawtoothClampPtr(SawtoothClamp):
    def __init__(self, this):
        _swig_setattr(self, SawtoothClamp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, SawtoothClamp, 'thisown', 0)
        _swig_setattr(self, SawtoothClamp,self.__class__,SawtoothClamp)
_cp4vasp.SawtoothClamp_swigregister(SawtoothClampPtr)

class CosClamp(Clamp):
    __swig_setmethods__ = {}
    for _s in [Clamp]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, CosClamp, name, value)
    __swig_getmethods__ = {}
    for _s in [Clamp]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, CosClamp, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ CosClamp instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, CosClamp, 'this', _cp4vasp.new_CosClamp(*args))
        _swig_setattr(self, CosClamp, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_CosClamp):
        try:
            if self.thisown: destroy(self)
        except: pass


class CosClampPtr(CosClamp):
    def __init__(self, this):
        _swig_setattr(self, CosClamp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, CosClamp, 'thisown', 0)
        _swig_setattr(self, CosClamp,self.__class__,CosClamp)
_cp4vasp.CosClamp_swigregister(CosClampPtr)

class WaveClamp(Clamp):
    __swig_setmethods__ = {}
    for _s in [Clamp]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, WaveClamp, name, value)
    __swig_getmethods__ = {}
    for _s in [Clamp]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, WaveClamp, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ WaveClamp instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, WaveClamp, 'this', _cp4vasp.new_WaveClamp(*args))
        _swig_setattr(self, WaveClamp, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_WaveClamp):
        try:
            if self.thisown: destroy(self)
        except: pass


class WaveClampPtr(WaveClamp):
    def __init__(self, this):
        _swig_setattr(self, WaveClamp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, WaveClamp, 'thisown', 0)
        _swig_setattr(self, WaveClamp,self.__class__,WaveClamp)
_cp4vasp.WaveClamp_swigregister(WaveClampPtr)

class AtanClamp(Clamp):
    __swig_setmethods__ = {}
    for _s in [Clamp]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, AtanClamp, name, value)
    __swig_getmethods__ = {}
    for _s in [Clamp]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, AtanClamp, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ AtanClamp instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, AtanClamp, 'this', _cp4vasp.new_AtanClamp(*args))
        _swig_setattr(self, AtanClamp, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_AtanClamp):
        try:
            if self.thisown: destroy(self)
        except: pass


class AtanClampPtr(AtanClamp):
    def __init__(self, this):
        _swig_setattr(self, AtanClamp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, AtanClamp, 'thisown', 0)
        _swig_setattr(self, AtanClamp,self.__class__,AtanClamp)
_cp4vasp.AtanClamp_swigregister(AtanClampPtr)

class FermiClamp(Clamp):
    __swig_setmethods__ = {}
    for _s in [Clamp]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FermiClamp, name, value)
    __swig_getmethods__ = {}
    for _s in [Clamp]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FermiClamp, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ FermiClamp instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, FermiClamp, 'this', _cp4vasp.new_FermiClamp(*args))
        _swig_setattr(self, FermiClamp, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_FermiClamp):
        try:
            if self.thisown: destroy(self)
        except: pass


class FermiClampPtr(FermiClamp):
    def __init__(self, this):
        _swig_setattr(self, FermiClamp, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FermiClamp, 'thisown', 0)
        _swig_setattr(self, FermiClamp,self.__class__,FermiClamp)
_cp4vasp.FermiClamp_swigregister(FermiClampPtr)

class ColorGradient(ClassInterface):
    __swig_setmethods__ = {}
    for _s in [ClassInterface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ColorGradient, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassInterface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ColorGradient, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ ColorGradient instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["color"] = _cp4vasp.ColorGradient_color_set
    __swig_getmethods__["color"] = _cp4vasp.ColorGradient_color_get
    if _newclass:color = property(_cp4vasp.ColorGradient_color_get, _cp4vasp.ColorGradient_color_set)
    def getClassName(*args): return _cp4vasp.ColorGradient_getClassName(*args)
    def f(*args): return _cp4vasp.ColorGradient_f(*args)
    def glcolor(*args): return _cp4vasp.ColorGradient_glcolor(*args)
    def __del__(self, destroy=_cp4vasp.delete_ColorGradient):
        try:
            if self.thisown: destroy(self)
        except: pass


class ColorGradientPtr(ColorGradient):
    def __init__(self, this):
        _swig_setattr(self, ColorGradient, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ColorGradient, 'thisown', 0)
        _swig_setattr(self, ColorGradient,self.__class__,ColorGradient)
_cp4vasp.ColorGradient_swigregister(ColorGradientPtr)

class GrayColorGradient(ColorGradient):
    __swig_setmethods__ = {}
    for _s in [ColorGradient]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, GrayColorGradient, name, value)
    __swig_getmethods__ = {}
    for _s in [ColorGradient]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, GrayColorGradient, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ GrayColorGradient instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getClassName(*args): return _cp4vasp.GrayColorGradient_getClassName(*args)
    def f(*args): return _cp4vasp.GrayColorGradient_f(*args)
    def __init__(self, *args):
        _swig_setattr(self, GrayColorGradient, 'this', _cp4vasp.new_GrayColorGradient(*args))
        _swig_setattr(self, GrayColorGradient, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_GrayColorGradient):
        try:
            if self.thisown: destroy(self)
        except: pass


class GrayColorGradientPtr(GrayColorGradient):
    def __init__(self, this):
        _swig_setattr(self, GrayColorGradient, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, GrayColorGradient, 'thisown', 0)
        _swig_setattr(self, GrayColorGradient,self.__class__,GrayColorGradient)
_cp4vasp.GrayColorGradient_swigregister(GrayColorGradientPtr)

class RainbowColorGradient(ColorGradient):
    __swig_setmethods__ = {}
    for _s in [ColorGradient]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, RainbowColorGradient, name, value)
    __swig_getmethods__ = {}
    for _s in [ColorGradient]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, RainbowColorGradient, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ RainbowColorGradient instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["saturation"] = _cp4vasp.RainbowColorGradient_saturation_set
    __swig_getmethods__["saturation"] = _cp4vasp.RainbowColorGradient_saturation_get
    if _newclass:saturation = property(_cp4vasp.RainbowColorGradient_saturation_get, _cp4vasp.RainbowColorGradient_saturation_set)
    __swig_setmethods__["value"] = _cp4vasp.RainbowColorGradient_value_set
    __swig_getmethods__["value"] = _cp4vasp.RainbowColorGradient_value_get
    if _newclass:value = property(_cp4vasp.RainbowColorGradient_value_get, _cp4vasp.RainbowColorGradient_value_set)
    def __init__(self, *args):
        _swig_setattr(self, RainbowColorGradient, 'this', _cp4vasp.new_RainbowColorGradient(*args))
        _swig_setattr(self, RainbowColorGradient, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.RainbowColorGradient_getClassName(*args)
    def f(*args): return _cp4vasp.RainbowColorGradient_f(*args)
    def __del__(self, destroy=_cp4vasp.delete_RainbowColorGradient):
        try:
            if self.thisown: destroy(self)
        except: pass


class RainbowColorGradientPtr(RainbowColorGradient):
    def __init__(self, this):
        _swig_setattr(self, RainbowColorGradient, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, RainbowColorGradient, 'thisown', 0)
        _swig_setattr(self, RainbowColorGradient,self.__class__,RainbowColorGradient)
_cp4vasp.RainbowColorGradient_swigregister(RainbowColorGradientPtr)

class VisSlideDrawer(VisDrawer):
    __swig_setmethods__ = {}
    for _s in [VisDrawer]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisSlideDrawer, name, value)
    __swig_getmethods__ = {}
    for _s in [VisDrawer]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisSlideDrawer, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisSlideDrawer instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["n1"] = _cp4vasp.VisSlideDrawer_n1_set
    __swig_getmethods__["n1"] = _cp4vasp.VisSlideDrawer_n1_get
    if _newclass:n1 = property(_cp4vasp.VisSlideDrawer_n1_get, _cp4vasp.VisSlideDrawer_n1_set)
    __swig_setmethods__["n2"] = _cp4vasp.VisSlideDrawer_n2_set
    __swig_getmethods__["n2"] = _cp4vasp.VisSlideDrawer_n2_get
    if _newclass:n2 = property(_cp4vasp.VisSlideDrawer_n2_get, _cp4vasp.VisSlideDrawer_n2_set)
    __swig_setmethods__["lo"] = _cp4vasp.VisSlideDrawer_lo_set
    __swig_getmethods__["lo"] = _cp4vasp.VisSlideDrawer_lo_get
    if _newclass:lo = property(_cp4vasp.VisSlideDrawer_lo_get, _cp4vasp.VisSlideDrawer_lo_set)
    __swig_setmethods__["hi"] = _cp4vasp.VisSlideDrawer_hi_set
    __swig_getmethods__["hi"] = _cp4vasp.VisSlideDrawer_hi_get
    if _newclass:hi = property(_cp4vasp.VisSlideDrawer_hi_get, _cp4vasp.VisSlideDrawer_hi_set)
    __swig_setmethods__["scale"] = _cp4vasp.VisSlideDrawer_scale_set
    __swig_getmethods__["scale"] = _cp4vasp.VisSlideDrawer_scale_get
    if _newclass:scale = property(_cp4vasp.VisSlideDrawer_scale_get, _cp4vasp.VisSlideDrawer_scale_set)
    def __init__(self, *args):
        _swig_setattr(self, VisSlideDrawer, 'this', _cp4vasp.new_VisSlideDrawer(*args))
        _swig_setattr(self, VisSlideDrawer, 'thisown', 1)
    def setShadow(*args): return _cp4vasp.VisSlideDrawer_setShadow(*args)
    def getShadow(*args): return _cp4vasp.VisSlideDrawer_getShadow(*args)
    def setB1(*args): return _cp4vasp.VisSlideDrawer_setB1(*args)
    def setB2(*args): return _cp4vasp.VisSlideDrawer_setB2(*args)
    def setOrigin(*args): return _cp4vasp.VisSlideDrawer_setOrigin(*args)
    def setFArray(*args): return _cp4vasp.VisSlideDrawer_setFArray(*args)
    def setGradient(*args): return _cp4vasp.VisSlideDrawer_setGradient(*args)
    def setClamp(*args): return _cp4vasp.VisSlideDrawer_setClamp(*args)
    def assureClampAndGradient(*args): return _cp4vasp.VisSlideDrawer_assureClampAndGradient(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisSlideDrawer):
        try:
            if self.thisown: destroy(self)
        except: pass


class VisSlideDrawerPtr(VisSlideDrawer):
    def __init__(self, this):
        _swig_setattr(self, VisSlideDrawer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisSlideDrawer, 'thisown', 0)
        _swig_setattr(self, VisSlideDrawer,self.__class__,VisSlideDrawer)
_cp4vasp.VisSlideDrawer_swigregister(VisSlideDrawerPtr)

class VisStructureArrowsDrawer(VisDrawer):
    __swig_setmethods__ = {}
    for _s in [VisDrawer]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisStructureArrowsDrawer, name, value)
    __swig_getmethods__ = {}
    for _s in [VisDrawer]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisStructureArrowsDrawer, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisStructureArrowsDrawer instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["arrow_radius"] = _cp4vasp.VisStructureArrowsDrawer_arrow_radius_set
    __swig_getmethods__["arrow_radius"] = _cp4vasp.VisStructureArrowsDrawer_arrow_radius_get
    if _newclass:arrow_radius = property(_cp4vasp.VisStructureArrowsDrawer_arrow_radius_get, _cp4vasp.VisStructureArrowsDrawer_arrow_radius_set)
    __swig_setmethods__["arrowhead_radius"] = _cp4vasp.VisStructureArrowsDrawer_arrowhead_radius_set
    __swig_getmethods__["arrowhead_radius"] = _cp4vasp.VisStructureArrowsDrawer_arrowhead_radius_get
    if _newclass:arrowhead_radius = property(_cp4vasp.VisStructureArrowsDrawer_arrowhead_radius_get, _cp4vasp.VisStructureArrowsDrawer_arrowhead_radius_set)
    __swig_setmethods__["arrowhead_length"] = _cp4vasp.VisStructureArrowsDrawer_arrowhead_length_set
    __swig_getmethods__["arrowhead_length"] = _cp4vasp.VisStructureArrowsDrawer_arrowhead_length_get
    if _newclass:arrowhead_length = property(_cp4vasp.VisStructureArrowsDrawer_arrowhead_length_get, _cp4vasp.VisStructureArrowsDrawer_arrowhead_length_set)
    __swig_setmethods__["red"] = _cp4vasp.VisStructureArrowsDrawer_red_set
    __swig_getmethods__["red"] = _cp4vasp.VisStructureArrowsDrawer_red_get
    if _newclass:red = property(_cp4vasp.VisStructureArrowsDrawer_red_get, _cp4vasp.VisStructureArrowsDrawer_red_set)
    __swig_setmethods__["green"] = _cp4vasp.VisStructureArrowsDrawer_green_set
    __swig_getmethods__["green"] = _cp4vasp.VisStructureArrowsDrawer_green_get
    if _newclass:green = property(_cp4vasp.VisStructureArrowsDrawer_green_get, _cp4vasp.VisStructureArrowsDrawer_green_set)
    __swig_setmethods__["blue"] = _cp4vasp.VisStructureArrowsDrawer_blue_set
    __swig_getmethods__["blue"] = _cp4vasp.VisStructureArrowsDrawer_blue_get
    if _newclass:blue = property(_cp4vasp.VisStructureArrowsDrawer_blue_get, _cp4vasp.VisStructureArrowsDrawer_blue_set)
    __swig_setmethods__["arrows_scale"] = _cp4vasp.VisStructureArrowsDrawer_arrows_scale_set
    __swig_getmethods__["arrows_scale"] = _cp4vasp.VisStructureArrowsDrawer_arrows_scale_get
    if _newclass:arrows_scale = property(_cp4vasp.VisStructureArrowsDrawer_arrows_scale_get, _cp4vasp.VisStructureArrowsDrawer_arrows_scale_set)
    def __init__(self, *args):
        _swig_setattr(self, VisStructureArrowsDrawer, 'this', _cp4vasp.new_VisStructureArrowsDrawer(*args))
        _swig_setattr(self, VisStructureArrowsDrawer, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.VisStructureArrowsDrawer_getClassName(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisStructureArrowsDrawer):
        try:
            if self.thisown: destroy(self)
        except: pass

    def updateStructure(*args): return _cp4vasp.VisStructureArrowsDrawer_updateStructure(*args)
    def len(*args): return _cp4vasp.VisStructureArrowsDrawer_len(*args)
    def getArrow(*args): return _cp4vasp.VisStructureArrowsDrawer_getArrow(*args)
    def setArrow(*args): return _cp4vasp.VisStructureArrowsDrawer_setArrow(*args)
    def setScale(*args): return _cp4vasp.VisStructureArrowsDrawer_setScale(*args)
    def getScale(*args): return _cp4vasp.VisStructureArrowsDrawer_getScale(*args)

class VisStructureArrowsDrawerPtr(VisStructureArrowsDrawer):
    def __init__(self, this):
        _swig_setattr(self, VisStructureArrowsDrawer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisStructureArrowsDrawer, 'thisown', 0)
        _swig_setattr(self, VisStructureArrowsDrawer,self.__class__,VisStructureArrowsDrawer)
_cp4vasp.VisStructureArrowsDrawer_swigregister(VisStructureArrowsDrawerPtr)

class VisIsosurfaceDrawer(VisDrawer):
    __swig_setmethods__ = {}
    for _s in [VisDrawer]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, VisIsosurfaceDrawer, name, value)
    __swig_getmethods__ = {}
    for _s in [VisDrawer]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, VisIsosurfaceDrawer, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ VisIsosurfaceDrawer instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["list"] = _cp4vasp.VisIsosurfaceDrawer_list_set
    __swig_getmethods__["list"] = _cp4vasp.VisIsosurfaceDrawer_list_get
    if _newclass:list = property(_cp4vasp.VisIsosurfaceDrawer_list_get, _cp4vasp.VisIsosurfaceDrawer_list_set)
    __swig_setmethods__["list_update_required"] = _cp4vasp.VisIsosurfaceDrawer_list_update_required_set
    __swig_getmethods__["list_update_required"] = _cp4vasp.VisIsosurfaceDrawer_list_update_required_get
    if _newclass:list_update_required = property(_cp4vasp.VisIsosurfaceDrawer_list_update_required_get, _cp4vasp.VisIsosurfaceDrawer_list_update_required_set)
    __swig_setmethods__["level"] = _cp4vasp.VisIsosurfaceDrawer_level_set
    __swig_getmethods__["level"] = _cp4vasp.VisIsosurfaceDrawer_level_get
    if _newclass:level = property(_cp4vasp.VisIsosurfaceDrawer_level_get, _cp4vasp.VisIsosurfaceDrawer_level_set)
    __swig_setmethods__["draw_as_points"] = _cp4vasp.VisIsosurfaceDrawer_draw_as_points_set
    __swig_getmethods__["draw_as_points"] = _cp4vasp.VisIsosurfaceDrawer_draw_as_points_get
    if _newclass:draw_as_points = property(_cp4vasp.VisIsosurfaceDrawer_draw_as_points_get, _cp4vasp.VisIsosurfaceDrawer_draw_as_points_set)
    __swig_setmethods__["mx"] = _cp4vasp.VisIsosurfaceDrawer_mx_set
    __swig_getmethods__["mx"] = _cp4vasp.VisIsosurfaceDrawer_mx_get
    if _newclass:mx = property(_cp4vasp.VisIsosurfaceDrawer_mx_get, _cp4vasp.VisIsosurfaceDrawer_mx_set)
    __swig_setmethods__["my"] = _cp4vasp.VisIsosurfaceDrawer_my_set
    __swig_getmethods__["my"] = _cp4vasp.VisIsosurfaceDrawer_my_get
    if _newclass:my = property(_cp4vasp.VisIsosurfaceDrawer_my_get, _cp4vasp.VisIsosurfaceDrawer_my_set)
    __swig_setmethods__["mz"] = _cp4vasp.VisIsosurfaceDrawer_mz_set
    __swig_getmethods__["mz"] = _cp4vasp.VisIsosurfaceDrawer_mz_get
    if _newclass:mz = property(_cp4vasp.VisIsosurfaceDrawer_mz_get, _cp4vasp.VisIsosurfaceDrawer_mz_set)
    __swig_setmethods__["chgcar"] = _cp4vasp.VisIsosurfaceDrawer_chgcar_set
    __swig_getmethods__["chgcar"] = _cp4vasp.VisIsosurfaceDrawer_chgcar_get
    if _newclass:chgcar = property(_cp4vasp.VisIsosurfaceDrawer_chgcar_get, _cp4vasp.VisIsosurfaceDrawer_chgcar_set)
    def paint_isosurface(*args): return _cp4vasp.VisIsosurfaceDrawer_paint_isosurface(*args)
    def handle_tetrahedron(*args): return _cp4vasp.VisIsosurfaceDrawer_handle_tetrahedron(*args)
    def handle_type1(*args): return _cp4vasp.VisIsosurfaceDrawer_handle_type1(*args)
    def handle_type2(*args): return _cp4vasp.VisIsosurfaceDrawer_handle_type2(*args)
    __swig_setmethods__["red"] = _cp4vasp.VisIsosurfaceDrawer_red_set
    __swig_getmethods__["red"] = _cp4vasp.VisIsosurfaceDrawer_red_get
    if _newclass:red = property(_cp4vasp.VisIsosurfaceDrawer_red_get, _cp4vasp.VisIsosurfaceDrawer_red_set)
    __swig_setmethods__["green"] = _cp4vasp.VisIsosurfaceDrawer_green_set
    __swig_getmethods__["green"] = _cp4vasp.VisIsosurfaceDrawer_green_get
    if _newclass:green = property(_cp4vasp.VisIsosurfaceDrawer_green_get, _cp4vasp.VisIsosurfaceDrawer_green_set)
    __swig_setmethods__["blue"] = _cp4vasp.VisIsosurfaceDrawer_blue_set
    __swig_getmethods__["blue"] = _cp4vasp.VisIsosurfaceDrawer_blue_get
    if _newclass:blue = property(_cp4vasp.VisIsosurfaceDrawer_blue_get, _cp4vasp.VisIsosurfaceDrawer_blue_set)
    def __init__(self, *args):
        _swig_setattr(self, VisIsosurfaceDrawer, 'this', _cp4vasp.new_VisIsosurfaceDrawer(*args))
        _swig_setattr(self, VisIsosurfaceDrawer, 'thisown', 1)
    def getClassName(*args): return _cp4vasp.VisIsosurfaceDrawer_getClassName(*args)
    def setMultiple(*args): return _cp4vasp.VisIsosurfaceDrawer_setMultiple(*args)
    def setMultiple1(*args): return _cp4vasp.VisIsosurfaceDrawer_setMultiple1(*args)
    def setMultiple2(*args): return _cp4vasp.VisIsosurfaceDrawer_setMultiple2(*args)
    def setMultiple3(*args): return _cp4vasp.VisIsosurfaceDrawer_setMultiple3(*args)
    def getMultiple1(*args): return _cp4vasp.VisIsosurfaceDrawer_getMultiple1(*args)
    def getMultiple2(*args): return _cp4vasp.VisIsosurfaceDrawer_getMultiple2(*args)
    def getMultiple3(*args): return _cp4vasp.VisIsosurfaceDrawer_getMultiple3(*args)
    def updateIsosurface(*args): return _cp4vasp.VisIsosurfaceDrawer_updateIsosurface(*args)
    def setLevel(*args): return _cp4vasp.VisIsosurfaceDrawer_setLevel(*args)
    def setChgcar(*args): return _cp4vasp.VisIsosurfaceDrawer_setChgcar(*args)
    def getLevel(*args): return _cp4vasp.VisIsosurfaceDrawer_getLevel(*args)
    def getDrawAsPoints(*args): return _cp4vasp.VisIsosurfaceDrawer_getDrawAsPoints(*args)
    def setDrawAsPoints(*args): return _cp4vasp.VisIsosurfaceDrawer_setDrawAsPoints(*args)
    def __del__(self, destroy=_cp4vasp.delete_VisIsosurfaceDrawer):
        try:
            if self.thisown: destroy(self)
        except: pass


class VisIsosurfaceDrawerPtr(VisIsosurfaceDrawer):
    def __init__(self, this):
        _swig_setattr(self, VisIsosurfaceDrawer, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, VisIsosurfaceDrawer, 'thisown', 0)
        _swig_setattr(self, VisIsosurfaceDrawer,self.__class__,VisIsosurfaceDrawer)
_cp4vasp.VisIsosurfaceDrawer_swigregister(VisIsosurfaceDrawerPtr)


ODP_parseString = _cp4vasp.ODP_parseString

ODP_parseFile = _cp4vasp.ODP_parseFile
class ODPNode(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPNode, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ODPNode, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPNode instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["doc"] = _cp4vasp.ODPNode_doc_set
    __swig_getmethods__["doc"] = _cp4vasp.ODPNode_doc_get
    if _newclass:doc = property(_cp4vasp.ODPNode_doc_get, _cp4vasp.ODPNode_doc_set)
    __swig_setmethods__["pos"] = _cp4vasp.ODPNode_pos_set
    __swig_getmethods__["pos"] = _cp4vasp.ODPNode_pos_get
    if _newclass:pos = property(_cp4vasp.ODPNode_pos_get, _cp4vasp.ODPNode_pos_set)
    def __init__(self, *args):
        _swig_setattr(self, ODPNode, 'this', _cp4vasp.new_ODPNode(*args))
        _swig_setattr(self, ODPNode, 'thisown', 1)
    ELEMENT_NODE = _cp4vasp.ODPNode_ELEMENT_NODE
    ATTRIBUTE_NODE = _cp4vasp.ODPNode_ATTRIBUTE_NODE
    TEXT_NODE = _cp4vasp.ODPNode_TEXT_NODE
    CDATA_SECTION_NODE = _cp4vasp.ODPNode_CDATA_SECTION_NODE
    ENTITY_REFERENCE_NODE = _cp4vasp.ODPNode_ENTITY_REFERENCE_NODE
    ENTITY_NODE = _cp4vasp.ODPNode_ENTITY_NODE
    PROCESSING_INSTRUCTION_NODE = _cp4vasp.ODPNode_PROCESSING_INSTRUCTION_NODE
    COMMENT_NODE = _cp4vasp.ODPNode_COMMENT_NODE
    DOCUMENT_NODE = _cp4vasp.ODPNode_DOCUMENT_NODE
    DOCUMENT_TYPE_NODE = _cp4vasp.ODPNode_DOCUMENT_TYPE_NODE
    DOCUMENT_FRAGMENT_NODE = _cp4vasp.ODPNode_DOCUMENT_FRAGMENT_NODE
    NOTATION_NODE = _cp4vasp.ODPNode_NOTATION_NODE
    def getNodeName(*args): return _cp4vasp.ODPNode_getNodeName(*args)
    def getNodeValue(*args): return _cp4vasp.ODPNode_getNodeValue(*args)
    def setNodeValue(*args): return _cp4vasp.ODPNode_setNodeValue(*args)
    def getNodeType(*args): return _cp4vasp.ODPNode_getNodeType(*args)
    def getParentNode(*args): return _cp4vasp.ODPNode_getParentNode(*args)
    def getChildNodes(*args): return _cp4vasp.ODPNode_getChildNodes(*args)
    def getFirstChild(*args): return _cp4vasp.ODPNode_getFirstChild(*args)
    def getLastChild(*args): return _cp4vasp.ODPNode_getLastChild(*args)
    def getPreviousSibling(*args): return _cp4vasp.ODPNode_getPreviousSibling(*args)
    def getNextSibling(*args): return _cp4vasp.ODPNode_getNextSibling(*args)
    def getAttributes(*args): return _cp4vasp.ODPNode_getAttributes(*args)
    def getOwnerDocument(*args): return _cp4vasp.ODPNode_getOwnerDocument(*args)
    def hasChildNodes(*args): return _cp4vasp.ODPNode_hasChildNodes(*args)
    def up(*args): return _cp4vasp.ODPNode_up(*args)
    def down(*args): return _cp4vasp.ODPNode_down(*args)
    def next(*args): return _cp4vasp.ODPNode_next(*args)
    def previous(*args): return _cp4vasp.ODPNode_previous(*args)
    def poschar(*args): return _cp4vasp.ODPNode_poschar(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPNode):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPNodePtr(ODPNode):
    def __init__(self, this):
        _swig_setattr(self, ODPNode, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPNode, 'thisown', 0)
        _swig_setattr(self, ODPNode,self.__class__,ODPNode)
_cp4vasp.ODPNode_swigregister(ODPNodePtr)

class ODPDOMImplementation(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPDOMImplementation, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ODPDOMImplementation, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPDOMImplementation instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def hasFeature(*args): return _cp4vasp.ODPDOMImplementation_hasFeature(*args)
    def __init__(self, *args):
        _swig_setattr(self, ODPDOMImplementation, 'this', _cp4vasp.new_ODPDOMImplementation(*args))
        _swig_setattr(self, ODPDOMImplementation, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ODPDOMImplementation):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPDOMImplementationPtr(ODPDOMImplementation):
    def __init__(self, this):
        _swig_setattr(self, ODPDOMImplementation, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPDOMImplementation, 'thisown', 0)
        _swig_setattr(self, ODPDOMImplementation,self.__class__,ODPDOMImplementation)
_cp4vasp.ODPDOMImplementation_swigregister(ODPDOMImplementationPtr)

class ODPDocument(ODPNode):
    __swig_setmethods__ = {}
    for _s in [ODPNode]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPDocument, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNode]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPDocument, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPDocument instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getDoctype(*args): return _cp4vasp.ODPDocument_getDoctype(*args)
    def getImplementation(*args): return _cp4vasp.ODPDocument_getImplementation(*args)
    def getDocumentElement(*args): return _cp4vasp.ODPDocument_getDocumentElement(*args)
    def createElement(*args): return _cp4vasp.ODPDocument_createElement(*args)
    def createDocumentFragment(*args): return _cp4vasp.ODPDocument_createDocumentFragment(*args)
    def createTextNode(*args): return _cp4vasp.ODPDocument_createTextNode(*args)
    def createComment(*args): return _cp4vasp.ODPDocument_createComment(*args)
    def createCDATASection(*args): return _cp4vasp.ODPDocument_createCDATASection(*args)
    def createProcessingInstruction(*args): return _cp4vasp.ODPDocument_createProcessingInstruction(*args)
    def createAttribute(*args): return _cp4vasp.ODPDocument_createAttribute(*args)
    def createEntityReference(*args): return _cp4vasp.ODPDocument_createEntityReference(*args)
    def getElementsByTagName(*args): return _cp4vasp.ODPDocument_getElementsByTagName(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPDocument):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPDocumentPtr(ODPDocument):
    def __init__(self, this):
        _swig_setattr(self, ODPDocument, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPDocument, 'thisown', 0)
        _swig_setattr(self, ODPDocument,self.__class__,ODPDocument)
_cp4vasp.ODPDocument_swigregister(ODPDocumentPtr)

class ODPDocumentParent(ODPDocument):
    __swig_setmethods__ = {}
    for _s in [ODPDocument]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPDocumentParent, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPDocument]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPDocumentParent, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPDocumentParent instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_cp4vasp.delete_ODPDocumentParent):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPDocumentParentPtr(ODPDocumentParent):
    def __init__(self, this):
        _swig_setattr(self, ODPDocumentParent, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPDocumentParent, 'thisown', 0)
        _swig_setattr(self, ODPDocumentParent,self.__class__,ODPDocumentParent)
_cp4vasp.ODPDocumentParent_swigregister(ODPDocumentParentPtr)

class ODPElement(ODPNode):
    __swig_setmethods__ = {}
    for _s in [ODPNode]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPElement, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNode]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPElement, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPElement instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPElement, 'this', _cp4vasp.new_ODPElement(*args))
        _swig_setattr(self, ODPElement, 'thisown', 1)
    def getTagName(*args): return _cp4vasp.ODPElement_getTagName(*args)
    def getAttribute(*args): return _cp4vasp.ODPElement_getAttribute(*args)
    def setAttribute(*args): return _cp4vasp.ODPElement_setAttribute(*args)
    def removeAttribute(*args): return _cp4vasp.ODPElement_removeAttribute(*args)
    def getAttributeNode(*args): return _cp4vasp.ODPElement_getAttributeNode(*args)
    def setAttributeNode(*args): return _cp4vasp.ODPElement_setAttributeNode(*args)
    def removeAttributeNode(*args): return _cp4vasp.ODPElement_removeAttributeNode(*args)
    def getElementsByTagName(*args): return _cp4vasp.ODPElement_getElementsByTagName(*args)
    def normalize(*args): return _cp4vasp.ODPElement_normalize(*args)
    def refreshAttr(*args): return _cp4vasp.ODPElement_refreshAttr(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPElement):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPElementPtr(ODPElement):
    def __init__(self, this):
        _swig_setattr(self, ODPElement, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPElement, 'thisown', 0)
        _swig_setattr(self, ODPElement,self.__class__,ODPElement)
_cp4vasp.ODPElement_swigregister(ODPElementPtr)

class ODPCharacterData(ODPNode):
    __swig_setmethods__ = {}
    for _s in [ODPNode]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPCharacterData, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNode]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPCharacterData, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPCharacterData instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getData(*args): return _cp4vasp.ODPCharacterData_getData(*args)
    def setData(*args): return _cp4vasp.ODPCharacterData_setData(*args)
    def __init__(self, *args):
        _swig_setattr(self, ODPCharacterData, 'this', _cp4vasp.new_ODPCharacterData(*args))
        _swig_setattr(self, ODPCharacterData, 'thisown', 1)
    def getLength(*args): return _cp4vasp.ODPCharacterData_getLength(*args)
    def substringData(*args): return _cp4vasp.ODPCharacterData_substringData(*args)
    def appendData(*args): return _cp4vasp.ODPCharacterData_appendData(*args)
    def insertData(*args): return _cp4vasp.ODPCharacterData_insertData(*args)
    def deleteData(*args): return _cp4vasp.ODPCharacterData_deleteData(*args)
    def replaceData(*args): return _cp4vasp.ODPCharacterData_replaceData(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPCharacterData):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPCharacterDataPtr(ODPCharacterData):
    def __init__(self, this):
        _swig_setattr(self, ODPCharacterData, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPCharacterData, 'thisown', 0)
        _swig_setattr(self, ODPCharacterData,self.__class__,ODPCharacterData)
_cp4vasp.ODPCharacterData_swigregister(ODPCharacterDataPtr)

class ODPText(ODPCharacterData):
    __swig_setmethods__ = {}
    for _s in [ODPCharacterData]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPText, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPCharacterData]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPText, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPText instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPText, 'this', _cp4vasp.new_ODPText(*args))
        _swig_setattr(self, ODPText, 'thisown', 1)
    def splitText(*args): return _cp4vasp.ODPText_splitText(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPText):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPTextPtr(ODPText):
    def __init__(self, this):
        _swig_setattr(self, ODPText, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPText, 'thisown', 0)
        _swig_setattr(self, ODPText,self.__class__,ODPText)
_cp4vasp.ODPText_swigregister(ODPTextPtr)

class ODPComment(ODPCharacterData):
    __swig_setmethods__ = {}
    for _s in [ODPCharacterData]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPComment, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPCharacterData]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPComment, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPComment instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPComment, 'this', _cp4vasp.new_ODPComment(*args))
        _swig_setattr(self, ODPComment, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ODPComment):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPCommentPtr(ODPComment):
    def __init__(self, this):
        _swig_setattr(self, ODPComment, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPComment, 'thisown', 0)
        _swig_setattr(self, ODPComment,self.__class__,ODPComment)
_cp4vasp.ODPComment_swigregister(ODPCommentPtr)

class ODPCDATASection(ODPText):
    __swig_setmethods__ = {}
    for _s in [ODPText]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPCDATASection, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPText]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPCDATASection, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPCDATASection instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPCDATASection, 'this', _cp4vasp.new_ODPCDATASection(*args))
        _swig_setattr(self, ODPCDATASection, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ODPCDATASection):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPCDATASectionPtr(ODPCDATASection):
    def __init__(self, this):
        _swig_setattr(self, ODPCDATASection, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPCDATASection, 'thisown', 0)
        _swig_setattr(self, ODPCDATASection,self.__class__,ODPCDATASection)
_cp4vasp.ODPCDATASection_swigregister(ODPCDATASectionPtr)

class ODPNodeList(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPNodeList, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ODPNodeList, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPNodeList instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def item(*args): return _cp4vasp.ODPNodeList_item(*args)
    def getLength(*args): return _cp4vasp.ODPNodeList_getLength(*args)
    def __init__(self, *args):
        _swig_setattr(self, ODPNodeList, 'this', _cp4vasp.new_ODPNodeList(*args))
        _swig_setattr(self, ODPNodeList, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ODPNodeList):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPNodeListPtr(ODPNodeList):
    def __init__(self, this):
        _swig_setattr(self, ODPNodeList, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPNodeList, 'thisown', 0)
        _swig_setattr(self, ODPNodeList,self.__class__,ODPNodeList)
_cp4vasp.ODPNodeList_swigregister(ODPNodeListPtr)

class ODPChildList(ODPNodeList):
    __swig_setmethods__ = {}
    for _s in [ODPNodeList]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPChildList, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNodeList]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPChildList, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPChildList instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPChildList, 'this', _cp4vasp.new_ODPChildList(*args))
        _swig_setattr(self, ODPChildList, 'thisown', 1)
    def item(*args): return _cp4vasp.ODPChildList_item(*args)
    def getLength(*args): return _cp4vasp.ODPChildList_getLength(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPChildList):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPChildListPtr(ODPChildList):
    def __init__(self, this):
        _swig_setattr(self, ODPChildList, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPChildList, 'thisown', 0)
        _swig_setattr(self, ODPChildList,self.__class__,ODPChildList)
_cp4vasp.ODPChildList_swigregister(ODPChildListPtr)

class ODPElementsByTagNameList(ODPNodeList):
    __swig_setmethods__ = {}
    for _s in [ODPNodeList]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPElementsByTagNameList, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNodeList]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPElementsByTagNameList, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPElementsByTagNameList instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPElementsByTagNameList, 'this', _cp4vasp.new_ODPElementsByTagNameList(*args))
        _swig_setattr(self, ODPElementsByTagNameList, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ODPElementsByTagNameList):
        try:
            if self.thisown: destroy(self)
        except: pass

    def item(*args): return _cp4vasp.ODPElementsByTagNameList_item(*args)
    def getLength(*args): return _cp4vasp.ODPElementsByTagNameList_getLength(*args)

class ODPElementsByTagNameListPtr(ODPElementsByTagNameList):
    def __init__(self, this):
        _swig_setattr(self, ODPElementsByTagNameList, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPElementsByTagNameList, 'thisown', 0)
        _swig_setattr(self, ODPElementsByTagNameList,self.__class__,ODPElementsByTagNameList)
_cp4vasp.ODPElementsByTagNameList_swigregister(ODPElementsByTagNameListPtr)

class ODPChildrenByTagNameList(ODPNodeList):
    __swig_setmethods__ = {}
    for _s in [ODPNodeList]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPChildrenByTagNameList, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNodeList]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPChildrenByTagNameList, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPChildrenByTagNameList instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPChildrenByTagNameList, 'this', _cp4vasp.new_ODPChildrenByTagNameList(*args))
        _swig_setattr(self, ODPChildrenByTagNameList, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ODPChildrenByTagNameList):
        try:
            if self.thisown: destroy(self)
        except: pass

    def item(*args): return _cp4vasp.ODPChildrenByTagNameList_item(*args)
    def getLength(*args): return _cp4vasp.ODPChildrenByTagNameList_getLength(*args)

class ODPChildrenByTagNameListPtr(ODPChildrenByTagNameList):
    def __init__(self, this):
        _swig_setattr(self, ODPChildrenByTagNameList, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPChildrenByTagNameList, 'thisown', 0)
        _swig_setattr(self, ODPChildrenByTagNameList,self.__class__,ODPChildrenByTagNameList)
_cp4vasp.ODPChildrenByTagNameList_swigregister(ODPChildrenByTagNameListPtr)

class ODPNamedNodeMap(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPNamedNodeMap, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ODPNamedNodeMap, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPNamedNodeMap instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def getNamedItem(*args): return _cp4vasp.ODPNamedNodeMap_getNamedItem(*args)
    def setNamedItem(*args): return _cp4vasp.ODPNamedNodeMap_setNamedItem(*args)
    def removeNamedItem(*args): return _cp4vasp.ODPNamedNodeMap_removeNamedItem(*args)
    def item(*args): return _cp4vasp.ODPNamedNodeMap_item(*args)
    def getLength(*args): return _cp4vasp.ODPNamedNodeMap_getLength(*args)
    def __init__(self, *args):
        _swig_setattr(self, ODPNamedNodeMap, 'this', _cp4vasp.new_ODPNamedNodeMap(*args))
        _swig_setattr(self, ODPNamedNodeMap, 'thisown', 1)
    def __del__(self, destroy=_cp4vasp.delete_ODPNamedNodeMap):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPNamedNodeMapPtr(ODPNamedNodeMap):
    def __init__(self, this):
        _swig_setattr(self, ODPNamedNodeMap, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPNamedNodeMap, 'thisown', 0)
        _swig_setattr(self, ODPNamedNodeMap,self.__class__,ODPNamedNodeMap)
_cp4vasp.ODPNamedNodeMap_swigregister(ODPNamedNodeMapPtr)

class ODPAttributeMap(ODPNamedNodeMap):
    __swig_setmethods__ = {}
    for _s in [ODPNamedNodeMap]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPAttributeMap, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNamedNodeMap]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPAttributeMap, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPAttributeMap instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPAttributeMap, 'this', _cp4vasp.new_ODPAttributeMap(*args))
        _swig_setattr(self, ODPAttributeMap, 'thisown', 1)
    def setNode(*args): return _cp4vasp.ODPAttributeMap_setNode(*args)
    def getAttribute(*args): return _cp4vasp.ODPAttributeMap_getAttribute(*args)
    def getNamedItem(*args): return _cp4vasp.ODPAttributeMap_getNamedItem(*args)
    def item(*args): return _cp4vasp.ODPAttributeMap_item(*args)
    def getLength(*args): return _cp4vasp.ODPAttributeMap_getLength(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPAttributeMap):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPAttributeMapPtr(ODPAttributeMap):
    def __init__(self, this):
        _swig_setattr(self, ODPAttributeMap, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPAttributeMap, 'thisown', 0)
        _swig_setattr(self, ODPAttributeMap,self.__class__,ODPAttributeMap)
_cp4vasp.ODPAttributeMap_swigregister(ODPAttributeMapPtr)

class ODPAttr(ODPNode):
    __swig_setmethods__ = {}
    for _s in [ODPNode]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODPAttr, name, value)
    __swig_getmethods__ = {}
    for _s in [ODPNode]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ODPAttr, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ODPAttr instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ODPAttr, 'this', _cp4vasp.new_ODPAttr(*args))
        _swig_setattr(self, ODPAttr, 'thisown', 1)
    def getName(*args): return _cp4vasp.ODPAttr_getName(*args)
    def getSpecified(*args): return _cp4vasp.ODPAttr_getSpecified(*args)
    def getValue(*args): return _cp4vasp.ODPAttr_getValue(*args)
    def setValue(*args): return _cp4vasp.ODPAttr_setValue(*args)
    def __del__(self, destroy=_cp4vasp.delete_ODPAttr):
        try:
            if self.thisown: destroy(self)
        except: pass


class ODPAttrPtr(ODPAttr):
    def __init__(self, this):
        _swig_setattr(self, ODPAttr, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODPAttr, 'thisown', 0)
        _swig_setattr(self, ODPAttr,self.__class__,ODPAttr)
_cp4vasp.ODPAttr_swigregister(ODPAttrPtr)


