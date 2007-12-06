import UnitCell as UC
from Atom import *

uc = UnitCell( )
at1=Atom(symbol='Fe', mass=57) ; pos1=(0.0,0.0,0.0)
at2=Atom(symbol='Al') ; pos2=(0.5,0.5,0.5)
    
site1 = Site(pos1, at1)
site2 = Site(pos2, at2)

uc.addAtom( at1, pos1, "Fe1" )
uc.addAtom( at2, pos2, "Al1" )

print uc
