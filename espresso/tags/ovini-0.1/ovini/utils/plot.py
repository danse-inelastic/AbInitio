#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

def createPHPlot(infile,  imagefile):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    from ovini.utils import utils

    # Populate 'x' and 'y' lists from file
    (e,  x) = utils.parsePHFile(infile)

    pyplot.plot(e, x, 'g')
    pyplot.xlabel('Energy')
    pyplot.ylabel('DOS')
    pyplot.xlim(0, 300)

    pyplot.savefig(imagefile)

def createPWPlot(infile,  imagefile):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    from ovini.utils import utils

    # Populate 'x', 'y' and 'z' (cumulative) lists from file
    (e,  x,  y,  z) = utils.parseFile(infile)

    pyplot.plot(e, x, 'b', e, y, 'g')
    pyplot.xlabel('Energy')
    pyplot.ylabel('DOS')
    pyplot.xlim(5, 25)

    pyplot.savefig(imagefile)

__date__ = "$Jul 30, 2009 12:09:43 PM$"


