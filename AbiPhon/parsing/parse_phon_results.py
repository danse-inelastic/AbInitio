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


def parsePhon2IDF(inputfilename='phon.out',
                  polarizationsfile='polarizations.idf',
                  omega2sfile='energies.idf',
                  D=3):
    """Parses the Phon output and writes a IDF-format file for polarization vectors.
    D is the dimensionality for the crystal; it should be equal to 3 for Phon outputs."""

    import numpy
    from idf.Polarizations import write as writePols
    from idf.Omega2 import write as writeOmega2s
    
    try:
        infile = open(inputfilename, 'r')
        #polfile = open(polarizationsfile, 'w')
        #om2file = open(omega2sfile, 'w')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    inputString = infile.read()
    res = []
    for data, dataStart, dataEnd in modesAtPoint.scanString(inputString):
        rt.append(data.asList())    

    nkpts = len(res)
    if nkpts == 0: raise ValueError, 'No phonon mode was parsed.'
    nmodes = len(res[0])
    natoms = nmodes / D
    pols = numpy.zeros((nkpts, nmodes, natoms, D), dtype='float')
    omega2s = numpy.zeros((nkpts, nmodes), dtype='float')
    for kIndex in range(nkpts):
        for modeIndex in range(nmodes):
            omega2s[kIndex][modeIndex] = res[kIndex][modeIndex]
            for atomIndex in range(natoms):
                polvec = numpy.array([x+1j*y for (x,y) in res[kIndex][modeIndex][atomIndex]])
                pols[kIndex,modeIndex,atomIndex] = polvec # this will crash if D's don't match

    writePols(pols, filename=polarizationsfile,
              comment='Parsed from Phon output file '+inputfilename)
    writeOmega2s(omega2s, filename=omega2sfile,
                 comment='Parsed from Phon output file '+inputfilename, D=3)
    return


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
