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


def parseModes(filename='phon.out'):
    """Parses the results of Phon calculation,
    for phonon eigenenergies and polarization vectors.
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


def parseQpoints(filename='QPOINTS',
                 dumpQpts=1,
                 dumpWeights=1,
                 outQfile='qpoints.pkl',
                 outWfile='weights.pkl'):
    """Parses the list of points in reciprocal space used by Phon.
    Returns two Numpy arrays:
    one for Qpoints and one for associated weights, eg:
    >>> qs, ws = parseQpoints()
    The points are obtained from the file QPOINTS by default.
    If dumpQpts is set to one, the Q-points are dumped in 'qponits.pkl' pickle file.
    If dumpWeights is set to one, the weights are dumped in a 'weights.pkl' pickle file.
    """
    try:
        f = open(filename, 'r')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)

    try:
        numpoints = int(f.readline().strip())
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    except ValueError:
        print "Could not convert data to an integer."
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise
    
    print "\n Number of Q-points: %s" % numpoints

    data = [map(float,line.split()) for line in f]

    data_array = np.array(data)
    qpoints = data_array[:,0:3]
    weights = data_array[:,3]

    if dumpQpts == 1:
        try:
            qptfile = open(outQfile, 'w')
        except IOError, (errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        qpoints.dump(qptfile)
        qptfile.close()

    if dumpWeights == 1:
        try:
            wghtfile = open(outWfile, 'w')
        except IOError, (errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        weights.dump(wghtfile)
        wghtfile.close()
    
    f.close()
    return qpoints, weights

def plotQpoints():
    """Plots the Qpoints used for phonon calculation."""
    import pylab as pl
    import matplotlib.axes3d as p3

    qpoints, weights = parseQpoints()

    fig = pl.figure()
    ax = p3.Axes3D(fig)
    ax.scatter3D(qpoints[:,0],qpoints[:,1],qpoints[:,2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    pl.show()


def test():
    results = parseModes("phon.out_dos_noIBZ")

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
