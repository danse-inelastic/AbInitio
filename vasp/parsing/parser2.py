# Olivier Delaire

"""This is a new implementation for a parser for VASP INCAR file.
The goal of this implementation is to suppress the duplication of
keyword entries in the INCAR file when keyword values are modified."""

import os
from OrderedDict import OrderedDict

from pyparsing import *

# numbers: 1, 30.0, 1e-5, -99
number = Combine( Optional('-') + ( '0' | Word('123456789',nums) ) + \
                  Optional( '.' + Word(nums) ) + \
                  Optional( Word('eE',exact=1) + Word(nums+'+-',nums) ) )

digits = "0123456789"
integer = Word(digits)

def convertNumbers(s,l,toks):
    n = toks[0]
    try: return int(n)
    except ValueError, ve: return float(n)
    raise

number.setParseAction( convertNumbers )
integer.setParseAction( convertNumbers )

lineCont = Suppress("\\")
tag = Word(alphas)
numValue = OneOrMore(Or([number, lineCont]))
#vaspOption = OneOrMore(Or([Word(alphas+"."), lineCont]))
# this does not work as the keyword from the next line
# may be interpreted as an option in continuing list of options.
# It should work for numerical values, however, as INCAR file lines
# may not begin with a number, so 'isolated' numbers must come from
# a continuing line.
vaspOption = Word(alphas+".")

numericalAssignment = tag + Suppress('=') + numValue
optionAssignment = tag + Suppress('=') + vaspOption

assignment = Or([numericalAssignment, optionAssignment])

class INPUT2(OrderedDict):
    '''
    We store all the key=vals in a dictionary to prevent
    having multiple entries.
    The INCAR file is scanned using the pyparsing package.
    The dictionnary is derived from OrderedDict to preserve the
    order in which the entries are added.
    '''
    def __init__(self,filename='INCAR'):
        OrderedDict.__init__(self)

        self.filename = filename

        # if the INCAR file does not exist, we create it:
        if not os.path.exists(filename):
               f = open(filename, 'w')
               f.close()

        f = open(filename,'r')
        incarString = f.read()
        f.close()

        incarlist = []

        for data, dataStart, dataEnd in assignment.scanString(incarString):
            incarlist.append(data)
            continue

        # The OrteredDict is populated from the incarlist created above.
        for item in incarlist:
            key = item[0]
            val = item[1:]
            # We call the OrderedDict __setitem__ to avoid updating the file at this point
            OrderedDict.__setitem__(self, key, val)
            continue

    def __setitem__(self,name,value):
        """modify the builtin setitem behavior so that everytime
        a parameter is added or changed the INCAR file is written out
        so it is always up to date.
        """
        OrderedDict.__setitem__(self,name,value)

        if __debug__:
            print 'Updating after setting an item!'
        self.Write()

    def Write(self):
        """
        write each key=val pair back out.
        This should handle gracefully list items (eg: MAGMOM for all atoms),
        by checking the type and writing the values in the list one by one.
        [ This is improved by using an OrderedDict class,
        that keeps track of the order in which items are
        added to the input dictionary.]
        """
        
        f = open(self.filename,'w')
        for key in self:
            # if we have a list of values, print those one by one
            if type(self[key]) == type([]):
                f.write(key + ' = ')
                for item in self[key]:
                    f.write(str(item) + ' ')
                    continue
                f.write('\n')
            else:
                f.write('%s = %s\n' % (key,self[key]))
            pass
        f.close()
                    


if __name__ == '__main__':

    c = INPUT2('tests/CO/INCAR')
    print c

    c2 = INPUT2('tests/INCARS/incar2')
    print c2

    c3 = INPUT2('tests/INCARS/incar3')
    print c3
    c3.Write('testincar3')
