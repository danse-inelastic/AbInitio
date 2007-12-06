#!/usr/bin/env python

__doc__ = """Parses the list of points in reciprocal space used by Phon to a Numpy array.
The points are obtained from the file QPOINTS."""

import numpy as np
import scipy.io
import pylab as pl
import matplotlib.axes3d as p3
import pickle

def parse_qpoints():

    try:
        f = open("QPOINTS", 'r')
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

    try:
        qptfile = open("qpoints.pkl", 'w')
        wghtfile = open("weights.pkl", "w")
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)

    qpoints.dump(qptfile)
    weights.dump(wghtfile)

    qptfile.close()
    wghtfile.close()
    

    fig = pl.figure()
    ax = p3.Axes3D(fig)
    ax.scatter3D(qpoints[:,0],qpoints[:,1],qpoints[:,2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    pl.show()

    f.close()
    return qpoints, weights




if __name__ == "__main__": parse_qpoints()
