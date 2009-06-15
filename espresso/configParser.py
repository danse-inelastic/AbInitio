#!/usr/bin/env python

# Parses configuration file (for pw.x) and stores it in dictionary.
# It also allows to dump the dictionary to create the configuration file
# The main use case is to parse the existing configuration file, change 
# some values and save it back to the configuration file.

# Input parameters are defined in INPUT_PW.html 

"""
Stability issues:
- Parsing goes line by line. 
- Namelist starts with '&' and ends with '/' on a separate line
- Card starts with card title on a separate line and values between card titles.
- Prints both Namelists and Cards in capital 
- Would like to use ordered dictionary?
- Refactoring?  Introduce class relation: Namelist(Block), Card(Block)
1. Regex: ()
"""


import configPW

ref = {'control': configPW.namelist_control}

namelistsPW = ('control',
             'system',
             'electrons',
             'ions',
             'cell',
             'phonon')

cardsPW = ('atomic_species', 
           'atomic_positions',
           'k_points automatic',
           'cell_parameters',
           'climbing_images',
           'constraints',
           'collective_vars',
           'occupations')

import sys

# Dictionary for Quantum Espresso configuration file.
# First record in the qe['name'] refers to 'type' of the block which can be either
#  'namelist' or 'card'. E.g. qe = {'control': {'type': 'namelist', ...}}
 
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

    lines       = f.readlines()
    lines       = clearLines(lines)
    marks       = getMarks(lines)
    
    for m in marks:
        block           = getBlock(m[0], lines[m[1]:m[2]]) # (type, slice)
        qe[block[0]]    = block[1]
                
    print qe
    
    f.close()

def clearLines(lines):
    cl = []     # Lines without white spaces and empty lines
    for l in lines:
        l = l.strip().strip(',') # Remove both lead and trailing whitespace, including '\n' and comma
        if l == '':
            continue
        
        cl.append(l)

    return cl

# Returns tuple (namelistName/cardName, parametersDictionary)
def getBlock(type, slice):
    # Improve calls?
    if type == 'namelist':
        return getNamelist(slice)
    elif type == 'card':
        return getCard(slice)
    
    return
    
# Returns (namelistName, parametersDictionary)
def getNamelist(slice):
    size    = len(slice)    # 8, numParams = 6
    name    = slice[0].strip('&')
    block   = {}
    block['type'] = 'namelist'
    
    for s in slice[1:]:
        p           = getParam(s)
        block[p[0]] = p[1]

    return (name, block)

# Returns (cardName, parametersDictionary)
def getCard(slice):
    name    = slice[0].lower()
    block   = {}
    block['type'] = 'card'
    vals    = []
    for s in slice[1:]:
        vals.append(s)
    block['values'] = vals

    return (name, block) 
    
    #Example return: ('atomic_species', {'type': 'card', 'values': ('Ni  26.98  Ni.pbe-nd-rrkjus.UPF', 'Other line', 'Another line')})


def getMarks(lines):
    """
    Determines start and end of namelist and card blocks: [type, start, end]
    E.g ['namelist', 0, 7] for CONTROL namelist
    Iterate over number of lines. Empty lines are included
    Not tested very well yet
    """
    blocklist   = []
    isNamelist  = False
    isCard      = False
    size        = len(lines)

    for i in range(size):
        l = lines[i]
        # We suppose that namelists and card do not intersect
        # Namelist part

        # Namelist end
        if l[0] == '/' and isNamelist:
            isNamelist  = False
            block.append(i)
            blocklist.append(block)

        # Namelist start
        if l[0] == '&' and not isNamelist:
            name = l[1:].lower()
            
            if not name in namelistsPW:
                continue             # namelist is not recognizable
            
            block       = []
            isNamelist  = True
            block.append('namelist')
            block.append(i)

        # Card part
        line    = l.lower()
        # Card end
        if line in cardsPW and isCard:
            #print "End: %s, line: %d" % (line, i-1)
            isCard  = False
            block.append(i)
            blocklist.append(block)

        if i == size-1 and isCard:
            isCard  = False
            block.append(i+1)
            blocklist.append(block)
                
        # Card start
        if line in cardsPW and not isCard:
            #print "Start: %s, line: %d" % (line, i)
            block   = []
            isCard  = True
            block.append('card')
            block.append(i)
                    
    return blocklist

    # Example return: [['namelist', 0, 7], ['namelist', 8, 20]]

# Takes string like 'a = 2' and returns tuple ('a', 2)
# Checks if the value is int, float or string
def getParam(s):
    ss = s.split('=')
    for i in range(len(ss)):
        ss[i] = ss[i].strip()
        
    # Strip comma in case if present
    val = ss[1].strip(',')
    
    # Assume that there are two values only: (variable, value) pair
    assert len(ss) == 2

    return (ss[0], val)

# Add parameter to namelist
def addNamelistParam(namelist, param, value):
    if namelist in namelistsPW:
        qe[namelist][param] = value
        
    return  # Do nothing

# Add record to card
def addCardParam(card, record):
    if card in cardsPW:
        qe[card]['values'].append(record)
        
    return  # Do nothing

# Remove parameter from namelist
def removeNamelistParam(namelist, param):
    if namelist in namelistsPW:
        try:
            del(qe[namelist][param])
        except KeyError:    # parameter is not present
            return    

# Remove parameter from card
def removeCard(card):
    try:
        qe[card]
        del(qe[card])
    except KeyError:    # parameter is not present
        return

# Edit namelist parameter in qe
def editNamelistParam(namelist, param, value):
    if namelist in namelistsPW:
        try:
            qe[namelist][param]     # Make sure that parameter exists
            qe[namelist][param] = value
        except KeyError:    # parameter is not present
            return
        
    return  # Do nothing

# Edit card parameter in qe. 'record' is list of values
def editCardParam(card, record):
    if card in cardsPW:
        qe[card]['values'] = record
        
    return  # Do nothing

# Saves the qe dictionary in the configuration file
def save(filename="config.saved"):
    nind    = "    "
    cind    = " "
    br      = "\n"
    s = ''
    f = open(filename, "w")
    namelists   = []
    cards       = []
    for e in qe.keys():
        if qe[e]['type'] == 'namelist': # namelist
            namelists.append(e)
        elif qe[e]['type'] == 'card':   # card
            cards.append(e)
    
    # 1. Dump namelists
    for n in namelists:
        s += "&%s%s" % (n.upper(), br)
        for np in qe[n].keys():
            if np == 'type':
                continue
            s += "%s%s = %s%s" % (nind, np, qe[n][np], br)
        s += "/%s" % br

    # 2. Dump cards
    for c in cards:
        s += "%s%s" % (c.upper(), br)
        
        for cp in qe[c]['values']:
            s += "%s%s%s" % (cind, cp, br)

    f.write(s)
    f.close()

def test():
    parse("vini/ni.scf.in")
    addNamelistParam('control', 'title', 'Ni')
    removeNamelistParam('control', 'title') #'verbosity')
    removeCard('atomic_species', )
    editNamelistParam('control', 'calculation', "'nscf'")
    editCardParam('atomic_positions', ['blah'])
    save("ni.scf.in.saved")

if __name__ == "__main__":
    test()
