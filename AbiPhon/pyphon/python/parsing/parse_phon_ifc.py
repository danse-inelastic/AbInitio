#!/usr/bin/env python

from pyparsing import *
import pickle

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
matrix33 = Group(vector3 + vector3 + vector3)

rcutoff = number
nbonds = number

header = rcutoff + Literal("cutoff radius") + nbonds + Literal("number of vectors")

kappa = number
kappa_ = number
nrr = number
weight = number

#bond = Suppress(Literal("vector:")) + Suppress(kappa) + Suppress(kappa_) + Suppress(nrr) + Suppress(weight)
bond = Suppress(Literal("vector:")) + kappa + kappa_ + nrr + weight + vector3
FCmatrix = Suppress(Literal("force constant matrix:")) + matrix33
bondFC = bond + FCmatrix

def parseIFC(filename):
    """Parses the interatomic force-constant tensor from the Phon output.
    """
    try:
        f = open(filename, 'r')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    inputString = f.read()

    info = header.parseString(inputString)

    rt = []
    for data, dataStart, dataEnd in bondFC.scanString(inputString):
        rt.append(data.asList())    
    return rt 

def test():
    results = parse("HARMONIC")

    try:
        outfile = open("harmonic.pkl", 'w')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    
    pickle.dump(results, outfile)
    outfile.close()
    return


def main():
    test()
    return


if __name__ == "__main__": main()
    
