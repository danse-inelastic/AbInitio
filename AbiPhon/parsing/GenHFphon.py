#!/usr/bin/env python

# this script extracts the atomic forces from the vasprun.xml files
# corresponding to different atomic displacements 
# and creates a 

from p4vasp.SystemPM import *
from p4vasp.Structure import Structure
#run=XMLSystemPM('vasprun.xml')
#struct=run.FINAL_STRUCTURE

struct = Structure("POSCAR")

dispfile=open("DISP",'r')
displist=[]
forcefile=open("FORCES", 'w')

for line in dispfile:
   str=eval(line)
   displist.append(str)

dispfile.close()
print "Number of displacements:", len(displist)
print "Displacements:"
for disp in displist:
   print disp.strip()

npos=len(displist)

forcefile.write(repr(npos)+"\n")

posnum=range(npos)

for i in posnum:
   run=XMLSystemPM('vasprun.xml_'+repr(i+1))
   forces=run.FORCES_SEQUENCE
   forcefile.write(displist[i])
   forcefile.write("\n")
   for atomforce in forces[0]:
      fx = "%8f" % atomforce[0]
      fy = "%8f" % atomforce[1]
      fz = "%8f" % atomforce[2]
      forcefile.write(fx+' '+fy+' '+fz+'\n')     

forcefile.close()
