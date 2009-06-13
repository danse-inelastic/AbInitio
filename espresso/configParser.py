#!/usr/bin/env python

# Parses configuration file (for pw.x) and stores it in dictionary.
# It also allows to dump the dictionary to create the configuration file
# The main use case is to parse the existing configuration file, change 
# some values and save it back to the configuration file.

# Input parameters are defined in INPUT_PW.html 

ref = {}


import sys

# Dictionary for Quantum Espresso configuration file. 
qe = QEConfig = {}

# Parses the configuration file and stores the values in qe dictionary
def parse(filename=None):
    if filename is None: 
        if len(sys.argv) != 2:
            print "Usage: configParser.py <config_file>"
            return
        else:
            filename = sys.argv[1]

    try:
        f = open(filename)
    except IOError:
        print "I/O error"
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise


    """
    - Parsing goes line by line. 
    - Namelist starts with '&' and ends with '/' on a separate line
    - Card starts with card title on a separate line and values between card titles.
    - Prints both Namelists and Cards in capital 
    1. Regex: ()
    """
    f.close()

# Adds parameter to qe
def add(card, parameter, value):
    print "stub: add"

# Removes the parameter from qe
def remove(card, parameter):
    print "stub: remove"

# Edits parameter in qe
def edit(card, parameter, value):
    print "stub: edit"

# Saves the qe dictionary in the configuration file
def save(filename=None):
    print "stub: save"

def test():
    parse()
    add('control', 'title', 'Ni')
    remove('control', 'verbosity')
    edit('control', 'calculation', 'nscf')
    save()
     

if __name__ == "__main__":
    test()
