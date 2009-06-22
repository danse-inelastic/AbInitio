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
    If dumpQpts is set to one, the Q-points are dumped in 'qpoints.pkl' pickle file.
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

def plotQpoints(qpoints):
    """Plots the q-points used for phonon calculation.
    The q-points can be parsed from the Phon output with:
    >>> qpoints, weights = parseQpoints(filename='QPOINTS')"""
    import pylab as pl
    import matplotlib.axes3d as p3

    fig = pl.figure()
    ax = p3.Axes3D(fig)
    ax.scatter3D(qpoints[:,0],qpoints[:,1],qpoints[:,2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    pl.show()



def parsePhon2IDF(inputfilename='phon.out',
                  polarizationsfile='polarizations.idf',
                  omega2sfile='energies.idf',
                  padding = False, qgridsize = None,
                  D=3):
    """Parses the Phon output and writes a IDF-format file for polarization vectors.
    D is the dimensionality for the crystal; it should be equal to 3 for Phon outputs."""

    from idf.Polarizations import write as writePols
    from idf.Omega2 import write as writeOmega2s

    pols, omega2s = readPols_and_Omega2s(inputfilename=inputfilename, D=D)
    if padding:
        if not qgridsize: raise ValueError, 'qgridsize not specified'
        pols, omega2s = applyPadding(pols, omega2s, qgridsize)
    
    writePols(pols, filename=polarizationsfile,
              comment='Parsed from Phon output file '+inputfilename)
    writeOmega2s(omega2s, filename=omega2sfile,
                 comment='Parsed from Phon output file '+inputfilename, D=3)
    return


def applyPadding(pols, omega2s, qgridsize, D=3):
    '''add paddings to the ends of q grids
    The Monkhorst Pack grid, when going through gamma point,
    covers 0, 1/n, 2/n, ..., (n-1)/n
    This padding makes the grid to cover 0, 1/n, 2/n, ..., (n-1)/n, 1
    by simply applying periodic condiction.
    '''
    import operator
    nqs = reduce(operator.__mul__, qgridsize)
    
    nqs1, nmodes, natoms, D = pols.shape
    assert nqs == nqs1 and nmodes == natoms*D

    nqs2, nmodes2 = omega2s.shape
    assert nqs2 == nqs1 and nmodes2 == nmodes2

    import numpy
    nqx,nqy,nqz = qgridsize
    newnqx = nqx+1
    newnqy = nqy+1
    newnqz = nqz+1
    newpols = numpy.zeros((newnqx, newnqy, newnqz, nmodes, natoms,D), dtype='complex')
    newomega2s = numpy.zeros((newnqx, newnqy, newnqz, nmodes), float)

    pols.shape = nqx,nqy,nqz,nmodes, natoms,D
    omega2s.shape = nqx,nqy,nqz,nmodes

    def pad(arr, newarr):
        newarr[:-1, :-1, :-1] = arr
        newarr[-1, :-1, :-1] = arr[0, :, :]
        newarr[:-1, -1, :-1] = arr[:, 0, :]
        newarr[:-1, :-1, -1] = arr[:, :, 0]
        newarr[-1, -1, :-1] = arr[0, 0, :]
        newarr[-1, :-1, -1] = arr[0, :, 0]
        newarr[:-1, -1, -1] = arr[:, 0, 0]
        newarr[-1, -1, -1] = arr[0, 0, 0]

    pad(pols, newpols)
    pad(omega2s, newomega2s)

    newpols.shape = -1, nmodes, natoms, D
    newomega2s.shape = -1, nmodes

    return newpols, newomega2s


def readPols_and_Omega2s(
    inputfilename='phon.out',
    D=3):
    """Parses the Phon output and read polarization vectors and omega2s (angular frequency sqaures)
    D is the dimensionality for the crystal; it should be equal to 3 for Phon outputs."""

    if D != 3: raise NotImplementedError
    
    import numpy
    
    try:
        infile = open(inputfilename, 'r')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    inputString = infile.read()
    res = []
    for data, dataStart, dataEnd in modesAtPoint.scanString(inputString):
        res.append(data.asList())    

    nkpts = len(res)
    if nkpts == 0: raise ValueError, 'No phonon mode was parsed.'
    nmodes = len(res[0])
    natoms = nmodes / D
    pols = numpy.zeros((nkpts, nmodes, natoms, D), dtype='complex')
    omega2s = numpy.zeros((nkpts, nmodes), dtype='float')
    twopi2 = (2*numpy.pi)**2
    for kIndex in range(nkpts):
        for modeIndex in range(nmodes):
            # we convert the omega^2 vaules from Phon (frequencies in Thz^2)
            # into angular frequncies in Hz^2
            omega2s[kIndex][modeIndex] = res[kIndex][modeIndex][0] * 1e24 * twopi2
            for atomIndex in range(natoms):
                polvec = numpy.array([x+1j*y for (x,y) in res[kIndex][modeIndex][atomIndex+1]])
                pols[kIndex,modeIndex,atomIndex] = polvec # this will crash if D's don't match

    return pols, omega2s


def parseFastPhon2IDF(inputfilename='phon.out',
                  polarizationsfile='polarizations.idf',
                  omega2sfile='energies.idf',
                  D=3):
    """Parses the Phon output and writes IDF-format files for the polarization vectors
    and for the eigenvalues (frequencies squared).
    D is the dimensionality for the crystal; it should be equal to 3 for Phon outputs.
    example:
    >>> parseFastPhon2IDF(inputfilename='phon_FeAl222.out',
                          polarizationsfile='polarizations.idf',
                          omega2sfile='energies.idf', D=3)
    """
    from idf.Polarizations import write as writePols
    from idf.Omega2 import write as writeOmega2s
    try:
        infile = open(inputfilename, 'r')
        #polfile = open(polarizationsfile, 'w')
        #om2file = open(omega2sfile, 'w')
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
    numkpts = 0
    res = []
    line=''
    # store eigenvals and eigenvecs into lists
    # these will be cast into numpy arrays,
    # once we know the dimensions of k-point grid and number of atoms
    eigvals = []
    eigvecs = []
    # !!! this does not work if Phon read the k-point grid from file !!!
    # Instead the number of k-points has to be passed explicitely
    # to parseFastPhon2IDF by the user.
    ## we want to determine the dimensions of the calculation:
    ## seek to the line defining the k-point grid:
    ##while line[0:10]!= 'Generating':
    ##    line=infile.readline()
    ##line = line.lstrip('Generating IBZ points ....')
    ##line = line.strip()
    ##griddims = [int(x) for x in line.split()]
    ##print 'Found a %s x %s x %s k-point grid.\n' % (griddims[0], griddims[1], griddims[2])

    # seek to first dynamical matrix (skip symmetry outputs):
    while line.strip() != 'Integrating frequencies...':
        line=infile.readline()
    infile.readline()  # skip blank line
    #infile.readline()  # skip 'Using ...'
    line = infile.readline() # read 'Using...' number of k-points
    stuff = line.split()
    numkpts = int(stuff[1])
    print 'Reading phonon modes at %s k-points.' % numkpts
    infile.readline()  # skip blank line
    # we are now at the first dynamical matrix.
    infile.readline()  # skip 'Dynamical matrix'
    numatoms = 1
    linecount = 0
    while line.strip() != 'Calling zheev.':
        linecount +=1
        line=infile.readline()
    linecount -= 1  # we read the extra 'Calling zheev' line.
    numatoms = int(np.sqrt(linecount/4)) # dyn mat is 3N * 3N, and each block is labelled
    print 'Found %s atom(s) \n' % numatoms
    # now we can read all eigenvalues and eigenvectors for all k-points:
    for ikpt in range(numkpts):
        # we just read in 'Calling zheev.'
        # loop over number of modes:
        modevals = []
        modevecs = []
        for modeindex in range(3*numatoms):
            infile.readline() # skip 'Eigenvalue  N'
            line = infile.readline()  # read eigenvalue
            modevals.append(float(line.strip()))
            infile.readline() # skip 'Eigenvector'
            vec = []
            for atomindex in range(numatoms):
                infile.readline() # skip 'Atom'
                atomvec = []
                for x in range(3):
                    vxstring = infile.readline()
                    vxlist = [float(x) for x in vxstring.strip().split()]
                    vx = complex(vxlist[0] + 1j * vxlist[1])
                    atomvec.append(vx)
                vec.append(atomvec)
            modevecs.append(vec)
        # we finished reading eigenvals and eigenvecs at current k-point
        eigvals.append(modevals)
        eigvecs.append(modevecs)
        #print "eigen-values:", eigvals
        # now skip next dynamical matrix:
        while ((line.strip() != 'Calling zheev.') and (line != '\n')):
            line=infile.readline()
    # write IDF files:
    omega2s = np.array(eigvals)
    pols = np.array(eigvecs)
    writeOmega2s(omega2s,
                 filename=omega2sfile,
                 comment='Parsed from'+inputfilename,
                 D=3)
    writePols(pols,
              filename=polarizationsfile,
              comment='Parsed from'+inputfilename)
    return
                        
                       
def parseDOS_meV2IDF(infile='DOS.meV', outfile='DOS_meV.idf'):
    f = open(infile)
    lines = f.readlines()
    es = []; Is = []
    for line in lines:
        line = line.strip()
        if not line: continue
        tokens = line.split()
        if len(tokens)!=2:
            #
            continue
        e, I = tokens
        e = float(e); es.append(e)
        I = float(I); Is.append(I)
        continue
    import numpy
    E = numpy.array(es)
    DOS = numpy.array(Is)
    from idf.DOS import write
    write(E, DOS, filename=outfile)
    return


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
