#!/usr/bin/env python


from pyparsing import *
import pickle

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
vector2 = Group(number + number)
vector6 = Group(number + number + number + number + number + number)
three_vectors = Group(vector3 + vector3 + vector3)

lattice_vectors = Suppress(Literal("Lattice vectors:")) + three_vectors
reciprocal_vectors = Suppress(Literal("Reciprocal vectors:")) + three_vectors

eigenvalue = Suppress(Literal("Eigenvalue")) + Suppress(number) + number
atomeigenvector = Group(vector2 + vector2 + vector2)
atomeigenvecoutput = Suppress(Literal("atom")) + Suppress(number) + atomeigenvector

multiAtomEigenvec = OneOrMore(atomeigenvecoutput)
mode = Group(eigenvalue + Suppress(Literal("Eigenvector")) + Suppress(number) + multiAtomEigenvec)
modesAtPoint = OneOrMore(mode)


def parse(filename):
    """Parses the results of Phon calculation for phonon eigenenergies and polarization vectors.
    The Phon output file is passed as argument."""

    try:
        f = open(filename, 'r')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    inputString = f.read()
    rt = []
    for data, dataStart, dataEnd in modesAtPoint.scanString(inputString):
        rt.append(data.asList())    
    return rt 

def test():
    results = parse("phon.out_dos_noIBZ")

    try:
        outfile = open("phonoutput.pkl", 'w')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    
    pickle.dump(results, outfile)
    outfile.close()
    return


def main():
    test()
    return


if __name__ == "__main__": main()
