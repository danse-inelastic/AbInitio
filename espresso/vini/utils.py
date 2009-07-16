import sys
import os
from subprocess import Popen,  PIPE

def run_pw_simulation(infile, outfile):
    f = open(infile)
    buf = f.read() # buf is 'string'
    f.close()
    proc = Popen("bin/pw.x", stdin=PIPE, stdout=PIPE,  stderr=PIPE,  shell="/bin/bash") 
    (stdout,   stderr) = proc.communicate(buf) # stdout is 'string'
    out = open(outfile,  "w")
    out.write(stdout)
    out.close()
    
def run_pw_dos(infile):
    f = open(infile)
    buf = f.read() # buf is 'string'
    f.close()
    proc = Popen("bin/dos.x", stdin=PIPE, stdout=PIPE,  stderr=PIPE,  shell="/bin/bash") 
    (stdout,   stderr) = proc.communicate(buf) # stdout is 'string'

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
    
def parse_file(filename):
    e = []
    x = []
    y = []
    z = []
    f = open(filename,  "r")
    line = f.readline() # Skip the first line with header 
    line = f.readline() 
    while line:
        list = line.split()
        #print list
        e.append(list[0])
        x.append(list[1])
        y.append(list[2])
        z.append(list[3])
        line = f.readline()
    f.close()
    return (e,  x,  y,  z)
    
    
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

def parse_ph_file(filename):
    e = []
    x = []
    f = open(filename,  "r")
    line = f.readline() 
    while line:
        list = line.split()
        #print list
        e.append(list[0])
        x.append(list[1])
        line = f.readline()
    f.close()
    return (e,  x)
