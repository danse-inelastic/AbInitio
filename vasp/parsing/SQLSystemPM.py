#!/usr/bin/python
#
# HappyDoc:docStringFormat='ClassicStructuredText'
#
#  p4vasp is a GUI-program and a library for processing outputs of the
#  Vienna Ab-inition Simulation Package (VASP)
#  (see http://cms.mpi.univie.ac.at/vasp/Welcome.html)
#  
#  Copyright (C) 2003  Orest Dubay <odubay@users.sourceforge.net>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



"""
SQLSystemPM - access to the data stored in a SQL database
"""

from __future__ import generators

from vasp.parsing import *
from vasp.parsing.util import *
import vasp.parsing.Dictionary
from vasp.parsing.store import *
from vasp.parsing.Array import *
from vasp.parsing.Property import *
import os.path
from string import *
from vasp.parsing.util import *
import traceback
import os
import cp4vasp
from vasp.parsing.Structure import *
import vasp.parsing.cStructure
import vasp.parsing.repository as repository
from vasp.parsing.SystemPM import *


class SQLStructuresLL(LateList):
  def __init__(self,plist,dbi):
    LateList.__init__(self,plist)
    self.dbi=dbi    
  def parse(self,x):
    print "parse",x  
    return self.dbi.readStructure(x)
    
class SQLSystemPM(SystemPM):
  def __init__(self,dbi,Id):
    SystemPM.__init__(self)
    self.dbi=dbi
    self.Id=Id
    self.add("NAME",                read = self.getName)
    self.add("INITIAL_STRUCTURE",   read = self.getInitialStructure)
    self.add("FINAL_STRUCTURE",     read = self.getFinalStructure)
    self.add("STRUCTURE_SEQUENCE",  read = self.getStructureSequence)
    self.add("STRUCTURE_SEQUENCE_L",read = self.getStructureSequenceL)
    self.add("FREE_ENERGY_SEQUENCE",read = self.getFreeEnergySequence)
    self.add("FREE_ENERGY",         read = lambda x,i=self.Id,d=self.dbi:d.fetchvalue(
                                           "SELECT energy FROM #CALC WHERE id=%d"%i))
  def getName(self,x):
    n=self.dbi.fetchvalue("SELECT name FROM #CALC WHERE id=%d"%s.Id)
    if n is None:
      n=""
    print "getName",n
    return n  
  def getInitialStructure(self,x):
    sid=self.dbi.fetchvalue("SELECT id FROM #STRUCT WHERE calc_id=%d AND step=-1"%self.Id)
    return self.dbi.readStructure(sid)
  def getFinalStructure(self,x):
    sid=self.dbi.fetchvalue("SELECT id FROM #STRUCT WHERE calc_id=%d AND step=-2"%self.Id)
    return self.dbi.readStructure(sid)
  def getStructureSequence(self,x):
    sids=self.dbi.fetchvalues("SELECT id FROM #STRUCT WHERE calc_id=%d AND step>=0 ORDER BY step"%self.Id)
    l=[]
    for i in sids:
      l.append(self.dbi.readStructure(i))
    return l

  def getFreeEnergySequence(self,x):
    return self.dbi.fetchvalues("SELECT energy FROM #ENERGY WHERE calc_id=%d AND step>=0 ORDER BY step"%self.Id)
  def getStructureSequenceL(self,x):
    sids=self.dbi.fetchvalues("SELECT id FROM #STRUCT WHERE calc_id=%d AND step>=0 ORDER BY step"%self.Id)
    return SQLStructuresLL(sids,self.dbi)
    
