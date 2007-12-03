#!/usr/bin/python

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
 setupstore serves for creation of setup StoreProfile
"""


from vasp.parsing.store import *
from string import *
from types import *

import vasp.parsing.applet.Applet 
import vasp.parsing.applet.EnergyConvergenceApplet
import vasp.parsing.applet.ForcesConvergenceApplet
import vasp.parsing.applet.SelForcesConvergenceApplet
import vasp.parsing.applet.StructureApplet

import vasp.parsing.applet
import vasp.parsing.Property
import vasp.parsing.SystemPM


          
class SetupProfile(Profile):
#  frame=None
  def __init__(self):
    Profile.__init__(self)
    self.addClass(vasp.parsing.applet.applets().getStoreProfile())
    self.addClass(vasp.parsing.SystemPM.systemlist().store_profile)

class SetupProfile_old(Profile):
#  frame=None
  def __init__(self):
    Profile.__init__(self)

    sp=vasp.parsing.applet.applets().store_profile
    
    for x in [
              vasp.parsing.applet.Applet.Applet,
              vasp.parsing.applet.EnergyConvergenceApplet.EnergyConvergenceApplet,
              vasp.parsing.applet.ForcesConvergenceApplet.ForcesConvergenceApplet,
              vasp.parsing.applet.SelForcesConvergenceApplet.SelForcesConvergenceApplet,
              vasp.parsing.applet.StructureApplet.StructureApplet,
             ]:
      #print "subprofile:",x
      sp.addClass(x.store_profile)


#if __name__=="__main__":
#  import LDOSApplet
#  from sys import *
#  f=stdout
#  sp=SetupProfile()
#  obj=LDOSApplet.LDOSApplet()
#  obj.lines.append(LDOSApplet.Line())
#  sp.dump(obj,f)
  
