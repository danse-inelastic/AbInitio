#!/usr/bin/env python

import pyphon, numpy

a = pyphon._pyphon.hpsort([1.,2,3,4,5])
print a
b = pyphon._pyphon.hpsort([1.,2,3,4,5,1])
c = pyphon._pyphon.hpsort([1.,2,3,4,5,1])
d = pyphon._pyphon.hpsort([1.,2,3,4,5,1,6])

print d

