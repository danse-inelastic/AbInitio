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

from vinil.utils.orderedDict import OrderedDict

# Available packages
PACKAGES    = ("Quantum Espresso",)  #, "VASP", "GULP"]  # Packages

# Type of configuration files
TYPES       = ("PW", "PH", "PP")  # "BANDS", "CPPP", "D3", "DOS", "DYNMAT", "INITIAL_STATE", "GIPAW", "D1", "MATDYN", "PROJWFC", "PWCOND", "Q2R" 

# Steps of job creation
STEPS       = ("Create Simulation",
               "Create Configuration",
               "Set Simulation Parameters",
               "Review Simulation")

# Types of simulations
SIMULATIONS = OrderedDict()
SIMULATIONS["Total Energy"]             = ("PW",)
SIMULATIONS["Electron DOS"]             = ("PW", "DOS")
SIMULATIONS["Electron Dispersion"]      = ("PW", "DOS")
SIMULATIONS["Geometry Optimization"]    = ("PW",)
SIMULATIONS["Single-Phonon"]            = ("PW", "PH", "DYNMAT")
SIMULATIONS["Multi-Phonon DOS"]         = ("PW", "PH", "Q2R", "MATDYN")
SIMULATIONS["Multi-Phonon Dispersion"]  = ("PW", "PH", "Q2R", "MATDYN")


# Available servers
SERVERS     = ("foxtrot.danse.us",)
               #"octopod.danse.us",
               #"upgrayedd.danse.us",
               #"teragrid"

# States of a job
STATES = {
        'C': 'finished',
        'R': 'running',
        'Q': 'queued',
        'E': 'exiting', #after having run
        'H': 'onhold',
        'W': 'waiting',
        'S': 'suspend',
        }



__date__ = "$Nov 3, 2009 3:12:34 PM$"


