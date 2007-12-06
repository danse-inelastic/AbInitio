#!/usr/bin/env python

import pyphon, numpy, os, sys

def read_input(filename = 'INPHON'):
    f = open(filename,'r').readlines()
    return f

def set_lsuper(lines, lsuper = '.T.'):
    def helper(ln):
        if ln.strip().upper()[0:6] == 'LSUPER':
            ln = ln.replace('.T.', lsuper)
            ln = ln.replace('.F.', lsuper)
            return ln
        else:
            return ln
    return [helper(ln) for ln in lines]

def run_phon(inputfile, phonexe = 'phon.exe'):
    f = read_input(inputfile)
    if os.path.exists('INPHON'):
        print "File INPHON already exists, copying to INPHON.bkup"
        import shutil
        shutil.copyfile('INPHON', 'INPHON.bkup')

    # first pass
    open('INPHON','w').write(''.join(set_lsuper(f, '.T.')))
    open('pass1.log',"w").write(''.join(os.popen(phonexe).readlines()))

    # second pass
    open('INPHON','w').write(''.join(set_lsuper(f, '.F.')))
    open('pass2.log',"w").write(''.join(os.popen(phonexe).readlines()))
        

if __name__ == '__main__':
    run_phon('INPHON')
    print "done"
