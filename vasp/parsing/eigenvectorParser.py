from pyparsing import *
import pickle
import numpy as np
#import scipy.io

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

vector3 = Group(number + number + number)
two_vectors = Group(Suppress(vector3) + vector3)

mode = OneOrMore(Suppress(Literal("X         Y         Z           dx          dy          dz")) + OneOrMore(two_vectors))

eandmode = number + Suppress(Literal('meV')) + mode

def fuse(listOfStrings):
    bigString=''
    for line in listOfStrings:
        bigString+=line
    return bigString

def getDynamicalMatrixOutput(filename='OUTCAR'):
    try:
        f = open(filename, 'r')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    lines=f.readlines()
    for i in range(len(lines)):
        if lines[i]==' Eigenvectors and eigenvalues of the dynamical matrix\n':
            start=i; continue
        if lines[i]==' Eigenvectors after division by SQRT(mass)\n':
            stop=i; break
    return fuse(lines[start:stop])
    
def parseModes():
    """Parses the result eigenvectors from a VASP calculation,
    The VASP output file is passed as argument."""
    inputString=getDynamicalMatrixOutput()
    rt = []
    dataSource=eandmode.scanString(inputString)
    #print dataSource
    for data, dataStart, dataEnd in dataSource:
        rt.append(data.asList())    
    return rt 

def parseEsModes():
    """Parses the result energies and eigenvectors from a VASP calculation,
    The VASP output file is passed as argument."""
    inputString=getDynamicalMatrixOutput()
    rt = []
    dataSource=eandmode.scanString(inputString)
    #print dataSource
    for data, dataStart, dataEnd in dataSource:
        rt.append(data.asList())    
    return rt 
    
