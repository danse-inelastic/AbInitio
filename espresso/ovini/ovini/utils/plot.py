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

def create_ph_plot(infile,  imagefile):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    # Populate 'x' and 'y' lists from file
    (e,  x) = parse_ph_file(infile)

    pyplot.plot(e, x, 'r')
    pyplot.xlabel('Energy')
    pyplot.ylabel('DOS')
    pyplot.xlim(0, 300)

    pyplot.savefig(imagefile)

def create_pw_plot(infile,  imagefile):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    # Populate 'x' and 'y' lists from file
    (e,  x,  y,  z) = parse_file(infile)

    pyplot.plot(e, x, 'b', e, y, 'g', e, z, 'r' )
    pyplot.xlabel('Energy')
    pyplot.ylabel('DOS')
    pyplot.xlim(5, 25)

    pyplot.savefig(imagefile)

__date__ = "$Jul 30, 2009 12:09:43 PM$"


